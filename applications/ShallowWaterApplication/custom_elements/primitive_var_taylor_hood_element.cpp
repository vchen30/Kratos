//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#include "primitive_var_taylor_hood_element.hpp"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/checks.h"
#include "shallow_water_application.h"

namespace Kratos {

Element::Pointer PrimitiveVarTaylorHoodElement::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new PrimitiveVarTaylorHoodElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties) );
}

int PrimitiveVarTaylorHoodElement::Check(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int result = Element::Check(rCurrentProcessInfo);
    if(result != 0) return result;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(HEIGHT)
    KRATOS_CHECK_VARIABLE_KEY(PROJECTED_VECTOR1)
    KRATOS_CHECK_VARIABLE_KEY(PROJECTED_SCALAR1)
    KRATOS_CHECK_VARIABLE_KEY(BATHYMETRY)
    KRATOS_CHECK_VARIABLE_KEY(GRAVITY)
    KRATOS_CHECK_VARIABLE_KEY(MANNING)
    KRATOS_CHECK_VARIABLE_KEY(RAIN)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for( unsigned int i = 0; i < this->GetGeometry().size(); ++i )
    {
        Node<3> &rNode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, rNode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, rNode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_VECTOR1, rNode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_SCALAR1, rNode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BATHYMETRY, rNode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RAIN, rNode)
    }

    // If this is a 2D problem, check that nodes are in XY plane
    //~ if (this->GetGeometry().WorkingSpaceDimension() == 2)
    //~ {
        //~ for (unsigned int i=0; i<this->GetGeometry().size(); ++i)
        //~ {
            //~ if (this->GetGeometry()[i].Z() != 0.0)
                //~ KRATOS_ERROR << "Node with non-zero Z coordinate found. Id: " << this->GetGeometry()[i].Id() << std:endl;
        //~ }
    //~ }

    return result;

    KRATOS_CATCH("")
}

void PrimitiveVarTaylorHoodElement::Initialize()
{
    KRATOS_TRY;

    const GeometryType& rGeom = this->GetGeometry();
    const SizeType Dim = rGeom.WorkingSpaceDimension();
    const SizeType NumVNodes = rGeom.PointsNumber();

    // Define a geometry container for pressure nodes
    switch (NumVNodes)
    {
    case 3: // 2D P1P1, not div-stable !!
        mpPressureGeometry = this->pGetGeometry();
        break;
    case 4: // 2D Q1Q1, not div-stable !!
        mpPressureGeometry = this->pGetGeometry();
        break;
    case 6: // 2D P2P1
        mpPressureGeometry = GeometryType::Pointer( new Triangle2D3< Node<3> >(rGeom(0), rGeom(1), rGeom(2)) );
        break;
    case 9: // 2D Q2Q1
        mpPressureGeometry = GeometryType::Pointer( new Quadrilateral2D4< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
        break;
    case 10: // 3D P2P1
        mpPressureGeometry = GeometryType::Pointer( new Tetrahedra3D4< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3)) );
        break;
    case 27: // 3D Q2Q1
        mpPressureGeometry = GeometryType::Pointer( new Hexahedra3D8< Node<3> >(rGeom(0), rGeom(1), rGeom(2), rGeom(3), rGeom(4), rGeom(5), rGeom(6), rGeom(7)) );
        break;
    default:
        KRATOS_ERROR << "Unexpected geometry type for Fluid Taylor-Hood elements" << std::endl;
    }

    if (NumVNodes > 4)
        this->mIntegrationMethod = GeometryData::GI_GAUSS_4; // Quadratic velocities
    else
        this->mIntegrationMethod = GeometryData::GI_GAUSS_2; // Linear velocities

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints( this->mIntegrationMethod );

    // Initialize member variables
    mDNv_DX.resize( IntegrationPoints.size() ); // Shape function derivatives container
    mDetJ.resize( IntegrationPoints.size() ); // determinant of Jacobian at each integration point

    GeometryType::JacobiansType J;
    J = GetGeometry().Jacobian( J, mIntegrationMethod );

    const GeometryType::ShapeFunctionsGradientsType& DNv_De = rGeom.ShapeFunctionsLocalGradients( this->mIntegrationMethod );

    // Temporary container for inverse of J
    Matrix InvJ;

    //calculating the inverse J
    for ( SizeType g = 0; g < IntegrationPoints.size(); g++ )
    {
        //calculating and storing inverse of the jacobian and the parameters needed
        MathUtils<double>::InvertMatrix( J[g], InvJ, mDetJ[g] );

        //calculating the shape function derivatives in global coordinates
        mDNv_DX[g].resize(NumVNodes,Dim);
        noalias( mDNv_DX[g] ) = prod( DNv_De[g], InvJ );
    }

    KRATOS_CATCH( "" )
}

void PrimitiveVarTaylorHoodElement::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix, VectorType &rRightHandSideVector, ProcessInfo &rCurrentProcessInfo)
{
    // Obtain required constants
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumVNodes = this->GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType NumGauss = this->GetGeometry().IntegrationPoints(this->mIntegrationMethod).size();

    const SizeType LocalSize = Dim * NumVNodes + NumPNodes;

    const Matrix NvContainer = this->GetGeometry().ShapeFunctionsValues(this->mIntegrationMethod);
    const Matrix NpContainer = mpPressureGeometry->ShapeFunctionsValues(this->mIntegrationMethod);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints( this->mIntegrationMethod );

    // Initialize local contribution
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    rRightHandSideVector = ZeroVector(LocalSize);

    // Loop on integration points
    for (SizeType g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& Nv = row(NvContainer,g);
        const ShapeFunctionsType& Np = row(NpContainer,g);
        const ShapeDerivativesType& DNv_DX = mDNv_DX[g];
        const double GaussWeight = mDetJ[g] * IntegrationPoints[g].Weight();

        double Density;
        double Viscosity;
        array_1d<double,3> BodyForce(3,0.0);
        array_1d<double,3> Velocity(3,0.0);
        array_1d<double,3> MeshVelocity(3,0.0);

        // Interpolation using pressure is linear
        this->EvaluateInPoint(Density,DENSITY,Np,*mpPressureGeometry);
        this->EvaluateInPoint(Viscosity,VISCOSITY,Np,*mpPressureGeometry);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,Np,*mpPressureGeometry);

        this->EvaluateInPoint(Velocity,VELOCITY,Nv,this->GetGeometry());
        this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,Nv,this->GetGeometry());

        // For ALE: convective velocity
        array_1d<double,3> ConvVel = Velocity - MeshVelocity;

        // Evaluate convection operator Velocity * Grad(N)
        Vector UGradN(NumVNodes);

        this->EvaluateConvection(UGradN,ConvVel,DNv_DX);

        // Add velocity terms in momentum equation
        this->AddMomentumTerms(rLeftHandSideMatrix,rRightHandSideVector,UGradN,Density,Viscosity,BodyForce,Nv,DNv_DX,GaussWeight);

        // Add velocity-pressure terms
        this->AddContinuityTerms(rLeftHandSideMatrix,Np,DNv_DX,GaussWeight);
    }

    // Add residual of previous iteration to RHS
    VectorType LastValues = ZeroVector(LocalSize);
    this->GetFirstDerivativesVector(LastValues);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,LastValues);
}

void PrimitiveVarTaylorHoodElement::MassMatrix(MatrixType &rMassMatrix,
                                 ProcessInfo &rCurrentProcessInfo)
{
    // Obtain required constants
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumVNodes = this->GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();
    const SizeType NumGauss = this->GetGeometry().IntegrationPoints(this->mIntegrationMethod).size();

    const SizeType LocalSize = Dim * NumVNodes + NumPNodes;

    const Matrix NvContainer = this->GetGeometry().ShapeFunctionsValues(this->mIntegrationMethod);
    const Matrix NpContainer = mpPressureGeometry->ShapeFunctionsValues(this->mIntegrationMethod);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = this->GetGeometry().IntegrationPoints( this->mIntegrationMethod );

    // Initialize local contribution
    if (rMassMatrix.size1() != LocalSize)
        rMassMatrix.resize(LocalSize, LocalSize, false);

    rMassMatrix = ZeroMatrix(LocalSize,LocalSize);

//    double Density = this->GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
//    const double A = 0.5 * mElementSize;
//    const double W = Density * A / 180.0;

//    rMassMatrix(0,0) = 6.0*W; rMassMatrix(0,2) = -1.0*W; rMassMatrix(0,4) = -1.0*W; rMassMatrix(0,8) = -4.0*W;
//    rMassMatrix(1,1) = 6.0*W; rMassMatrix(1,3) = -1.0*W; rMassMatrix(1,5) = -1.0*W; rMassMatrix(1,9) = -4.0*W;
//    rMassMatrix(2,0) = -1.0*W; rMassMatrix(2,2) = 6.0*W; rMassMatrix(2,4) = -1.0*W; rMassMatrix(2,10) = -4.0*W;
//    rMassMatrix(3,1) = -1.0*W; rMassMatrix(3,3) = 6.0*W; rMassMatrix(3,5) = -1.0*W; rMassMatrix(3,11) = -4.0*W;
//    rMassMatrix(4,0) = -1.0*W; rMassMatrix(4,2) = -1.0*W; rMassMatrix(4,4) = -6.0*W; rMassMatrix(4,6) = -4.0*W;
//    rMassMatrix(5,1) = -1.0*W; rMassMatrix(5,3) = -1.0*W; rMassMatrix(5,5) = -6.0*W; rMassMatrix(5,7) = -4.0*W;
//    rMassMatrix(6,4) = -4.0*W; rMassMatrix(6,6) = 32.0 * W; rMassMatrix(6,8) = 16.0 * W; rMassMatrix(6,10) = 16.0 * W;
//    rMassMatrix(7,5) = -4.0*W; rMassMatrix(7,7) = 32.0 * W; rMassMatrix(7,9) = 16.0 * W; rMassMatrix(7,11) = 16.0 * W;
//    rMassMatrix(8,0) = -4.0*W; rMassMatrix(8,6) = 16.0 * W; rMassMatrix(8,8) = 32.0 * W; rMassMatrix(8,10) = 16.0 * W;
//    rMassMatrix(9,1) = -4.0*W; rMassMatrix(9,7) = 16.0 * W; rMassMatrix(9,9) = 32.0 * W; rMassMatrix(9,11) = 16.0 * W;
//    rMassMatrix(10,2) = -4.0*W; rMassMatrix(10,6) = 16.0 * W; rMassMatrix(10,8) = 16.0 * W; rMassMatrix(10,10) = 32.0 * W;
//    rMassMatrix(11,3) = -4.0*W; rMassMatrix(11,7) = 16.0 * W; rMassMatrix(11,9) = 16.0 * W; rMassMatrix(11,11) = 32.0 * W;

//    double Density = this->GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
//    const double A = 0.5 * mElementSize;
//    const double W = Density * A / 57.0;
//    rMassMatrix(0,0) = 3.0*W;
//    rMassMatrix(1,1) = 3.0*W;
//    rMassMatrix(2,2) = 3.0*W;
//    rMassMatrix(3,3) = 3.0*W;
//    rMassMatrix(4,4) = 3.0*W;
//    rMassMatrix(5,5) = 3.0*W;
//    rMassMatrix(6,6) = 16.0*W;
//    rMassMatrix(7,7) = 16.0*W;
//    rMassMatrix(8,8) = 16.0*W;
//    rMassMatrix(9,9) = 16.0*W;
//    rMassMatrix(10,10) = 16.0*W;
//    rMassMatrix(11,11) = 16.0*W;

    // Loop on integration points
    for (SizeType g = 0; g < NumGauss; g++)
    {
        const ShapeFunctionsType& Nv = row(NvContainer,g);
        const ShapeFunctionsType& Np = row(NpContainer,g);

        double Density;
        // Interpolation using pressure is linear
        this->EvaluateInPoint(Density,DENSITY,Np,*mpPressureGeometry);

        const double Weight = Density * mDetJ[g] * IntegrationPoints[g].Weight();

        this->AddMassTerm(rMassMatrix,Nv,Weight);
    }
}


void PrimitiveVarTaylorHoodElement::GetDofList(DofsVectorType &rElementalDofList,
                                 ProcessInfo &rCurrentProcessInfo)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumVNodes = this->GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    const SizeType LocalSize = NumVNodes * Dim + NumPNodes;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumVNodes; i++)
    {
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof(VELOCITY_X);
        rElementalDofList[Index++] = GetGeometry()[i].pGetDof(VELOCITY_Y);
        if(Dim > 2) rElementalDofList[Index++] = GetGeometry()[i].pGetDof(VELOCITY_Z);
    }

    for (SizeType i = 0; i < NumPNodes; i++)
        rElementalDofList[Index++] = mpPressureGeometry->operator[](i).pGetDof(PRESSURE);
}


void PrimitiveVarTaylorHoodElement::EquationIdVector(Element::EquationIdVectorType &rResult, ProcessInfo &rCurrentProcessInfo)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumVNodes = this->GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    const SizeType LocalSize = NumVNodes * Dim + NumPNodes;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumVNodes; i++)
    {
        rResult[Index++] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
        rResult[Index++] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        if(Dim > 2) rResult[Index++] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
    }

    for (SizeType i = 0; i < NumPNodes; i++)
        rResult[Index++] = mpPressureGeometry->operator[](i).GetDof(PRESSURE).EquationId();
}

//void PrimitiveVarTaylorHoodElement::GetValuesVector(Vector &rValues, int Step)
//{
//    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
//    const SizeType NumVNodes = this->GetGeometry().PointsNumber();
//    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

//    const SizeType LocalSize = NumVNodes * Dim + NumPNodes;

//    if (rValues.size() != LocalSize)
//        rValues.resize(LocalSize);

//    SizeType Index = 0;

//    for (SizeType i = 0; i < NumVNodes; i++)
//    {
//        rValues[Index++] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X,Step);
//        rValues[Index++] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
//        if(Dim > 2) rValues[Index++] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Z,Step);
//    }

//    for (SizeType i = 0; i < NumPNodes; i++)
//        rValues[Index++] = mpPressureGeometry->operator[](i).FastGetSolutionStepValue(PRESSURE,Step);
//}

void PrimitiveVarTaylorHoodElement::GetFirstDerivativesVector(Vector &rValues, int Step)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumVNodes = this->GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    const SizeType LocalSize = NumVNodes * Dim + NumPNodes;

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumVNodes; i++)
    {
        rValues[Index++] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        rValues[Index++] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
        if(Dim > 2) rValues[Index++] = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Z,Step);
    }

    for (SizeType i = 0; i < NumPNodes; i++)
        rValues[Index++] = mpPressureGeometry->operator[](i).FastGetSolutionStepValue(PRESSURE,Step);
}

void PrimitiveVarTaylorHoodElement::GetSecondDerivativesVector(Vector &rValues, int Step)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumVNodes = this->GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    const SizeType LocalSize = NumVNodes * Dim + NumPNodes;

    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumVNodes; i++)
    {
        rValues[Index++] = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_X,Step);
        rValues[Index++] = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_Y,Step);
        if(Dim > 2) rValues[Index++] = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION_Z,Step);
    }

    for (SizeType i = 0; i < NumPNodes; i++)
        rValues[Index++] = 0.0;
}

void PrimitiveVarTaylorHoodElement::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType NumVNodes = this->GetGeometry().PointsNumber();
    GeometryType& rGeom = this->GetGeometry();
    switch (NumVNodes)
    {
    case 3: // P1P1, both geometries have the same nodes, do nothing
    {
        break;
    }
    case 4: // Q1Q1, both geometries have the same nodes, do nothing
    {
        break;
    }
    case 6: // triangle
    {
        const double p0 = rGeom[0].FastGetSolutionStepValue(PRESSURE);
        const double p1 = rGeom[1].FastGetSolutionStepValue(PRESSURE);
        const double p2 = rGeom[2].FastGetSolutionStepValue(PRESSURE);
        ThreadSafeNodeWrite(rGeom[3],PRESSURE, 0.5 * (p0 + p1) );
        ThreadSafeNodeWrite(rGeom[4],PRESSURE, 0.5 * (p1 + p2) );
        ThreadSafeNodeWrite(rGeom[5],PRESSURE, 0.5 * (p2 + p0) );
        break;
    }
    case 9: // quadrilateral
    {
        const double p0 = rGeom[0].FastGetSolutionStepValue(PRESSURE);
        const double p1 = rGeom[1].FastGetSolutionStepValue(PRESSURE);
        const double p2 = rGeom[2].FastGetSolutionStepValue(PRESSURE);
        const double p3 = rGeom[3].FastGetSolutionStepValue(PRESSURE);
        ThreadSafeNodeWrite(rGeom[4],PRESSURE, 0.5 * (p0 + p1) );
        ThreadSafeNodeWrite(rGeom[5],PRESSURE, 0.5 * (p1 + p2) );
        ThreadSafeNodeWrite(rGeom[6],PRESSURE, 0.5 * (p2 + p3) );
        ThreadSafeNodeWrite(rGeom[7],PRESSURE, 0.5 * (p3 + p0) );
        ThreadSafeNodeWrite(rGeom[8],PRESSURE, 0.25 * (p0 + p1 + p2 + p3) );
        break;
    }
    case 10: // tetrahedron
    {
        const double p0 = rGeom[0].FastGetSolutionStepValue(PRESSURE);
        const double p1 = rGeom[1].FastGetSolutionStepValue(PRESSURE);
        const double p2 = rGeom[2].FastGetSolutionStepValue(PRESSURE);
        const double p3 = rGeom[3].FastGetSolutionStepValue(PRESSURE);
        ThreadSafeNodeWrite(rGeom[4],PRESSURE, 0.5 * (p0 + p1) );
        ThreadSafeNodeWrite(rGeom[5],PRESSURE, 0.5 * (p1 + p2) );
        ThreadSafeNodeWrite(rGeom[6],PRESSURE, 0.5 * (p2 + p0) );
        ThreadSafeNodeWrite(rGeom[7],PRESSURE, 0.5 * (p0 + p3) );
        ThreadSafeNodeWrite(rGeom[8],PRESSURE, 0.5 * (p1 + p3) );
        ThreadSafeNodeWrite(rGeom[9],PRESSURE, 0.5 * (p2 + p3) );
        break;
    }
    case 27: // hexahedron
    {
        const double p0 = rGeom[0].FastGetSolutionStepValue(PRESSURE);
        const double p1 = rGeom[1].FastGetSolutionStepValue(PRESSURE);
        const double p2 = rGeom[2].FastGetSolutionStepValue(PRESSURE);
        const double p3 = rGeom[3].FastGetSolutionStepValue(PRESSURE);
        const double p4 = rGeom[4].FastGetSolutionStepValue(PRESSURE);
        const double p5 = rGeom[5].FastGetSolutionStepValue(PRESSURE);
        const double p6 = rGeom[6].FastGetSolutionStepValue(PRESSURE);
        const double p7 = rGeom[7].FastGetSolutionStepValue(PRESSURE);
        // edges -- bottom
        ThreadSafeNodeWrite(rGeom[8],PRESSURE, 0.5 * (p0 + p1) );
        ThreadSafeNodeWrite(rGeom[9],PRESSURE, 0.5 * (p1 + p2) );
        ThreadSafeNodeWrite(rGeom[10],PRESSURE, 0.5 * (p2 + p3) );
        ThreadSafeNodeWrite(rGeom[11],PRESSURE, 0.5 * (p3 + p0) );
        // edges -- middle
        ThreadSafeNodeWrite(rGeom[12],PRESSURE, 0.5 * (p4 + p0) );
        ThreadSafeNodeWrite(rGeom[13],PRESSURE, 0.5 * (p5 + p1) );
        ThreadSafeNodeWrite(rGeom[14],PRESSURE, 0.5 * (p6 + p2) );
        ThreadSafeNodeWrite(rGeom[15],PRESSURE, 0.5 * (p7 + p3) );
        // edges -- top
        ThreadSafeNodeWrite(rGeom[16],PRESSURE, 0.5 * (p4 + p5) );
        ThreadSafeNodeWrite(rGeom[17],PRESSURE, 0.5 * (p5 + p6) );
        ThreadSafeNodeWrite(rGeom[18],PRESSURE, 0.5 * (p6 + p7) );
        ThreadSafeNodeWrite(rGeom[19],PRESSURE, 0.5 * (p7 + p0) );
        // face centers
        ThreadSafeNodeWrite(rGeom[20],PRESSURE, 0.25 * (p0 + p1 + p2 + p3) );
        ThreadSafeNodeWrite(rGeom[21],PRESSURE, 0.25 * (p0 + p1 + p4 + p5) );
        ThreadSafeNodeWrite(rGeom[22],PRESSURE, 0.25 * (p1 + p2 + p5 + p6) );
        ThreadSafeNodeWrite(rGeom[23],PRESSURE, 0.25 * (p2 + p3 + p6 + p7) );
        ThreadSafeNodeWrite(rGeom[24],PRESSURE, 0.25 * (p3 + p0 + p7 + p4) );
        ThreadSafeNodeWrite(rGeom[25],PRESSURE, 0.25 * (p4 + p5 + p6 + p7) );
        // element center
        ThreadSafeNodeWrite(rGeom[26],PRESSURE, 0.125 * (p0+p1+p2+p3+p4+p5+p6+p7) );
        break;
    }
    default:
        KRATOS_ERROR << "Unexpected geometry type for Fluid Taylor-Hood elements" << std::endl;
    }

    KRATOS_CATCH("");
}

void PrimitiveVarTaylorHoodElement::AddMassTerm(MatrixType &rMassMatrix,
                                  const ShapeFunctionsType &Nv,
                                  const double Weight)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumVNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;
    SizeType FirstCol = 0;
    double Term_ij = 0.0;

    for (unsigned int i = 0; i < NumVNodes; ++i)
    {
        for (unsigned int j = 0; j < NumVNodes; ++j)
        {
            Term_ij = Weight * Nv[i] * Nv[j];

            for (unsigned int d = 0; d < Dim; ++d)
                rMassMatrix(FirstRow+d,FirstCol+d) += Term_ij;
            FirstCol += Dim;
        }
        FirstRow += Dim;
        FirstCol = 0;
    }
}

void PrimitiveVarTaylorHoodElement::AddMomentumTerms(MatrixType &rLHS,
                                       VectorType &rRHS,
                                       const Vector &UGradN,
                                       const double Density,
                                       const double Viscosity,
                                       const array_1d<double, 3> &BodyForce,
                                       const ShapeFunctionsType &Nv,
                                       const ShapeDerivativesType &DNv_DX,
                                       const double Weigth)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumVNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;
    SizeType FirstCol = 0;
    double Term_ij = 0.0;

    for(SizeType i = 0; i < NumVNodes; ++i)
    {
        // Body force
        for(SizeType d = 0; d < Dim; ++d)
            rRHS[FirstRow + d] += Weigth * Nv[i] * Density * BodyForce[d];

        for(SizeType j = 0; j < NumVNodes; ++j)
        {
            Term_ij = 0.0;

            // Viscous term
            for(SizeType d = 0; d < Dim; ++d)
                Term_ij += DNv_DX(i,d) * DNv_DX(j,d);
            Term_ij *= Viscosity;

            // Convection
            Term_ij += Nv[i] * UGradN[j];

            Term_ij *= Density * Weigth;

            for(SizeType d = 0; d < Dim; ++d)
                rLHS(FirstRow+d,FirstCol+d) += Term_ij;

            // Update column index
            FirstCol += Dim;
        }
        // Update matrix indices
        FirstRow += Dim;
        FirstCol = 0;
    }
}

void PrimitiveVarTaylorHoodElement::AddContinuityTerms(MatrixType &rLHS,
                                         const ShapeFunctionsType &Np,
                                         const ShapeDerivativesType &DNv_DX,
                                         const double Weight)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumVNodes = this->GetGeometry().PointsNumber();
    const SizeType NumPNodes = mpPressureGeometry->PointsNumber();

    SizeType Row = NumVNodes*Dim;
    SizeType FirstCol = 0;

    double DivTerm = 0.0;

    for (SizeType i = 0; i < NumPNodes; ++i)
    {
        for (SizeType j = 0; j < NumVNodes; ++j)
        {
            for (SizeType d = 0; d < Dim; ++d)
            {
                DivTerm = Weight * Np[i] * DNv_DX(j,d);
                rLHS(Row,FirstCol + d) += DivTerm; // Divergence term
                rLHS(FirstCol + d,Row) -= DivTerm; // Gradient term
            }

            // Update column index
            FirstCol += Dim;
        }
        // Update matrix indices
        FirstCol = 0;
        Row += 1;
    }
}

void PrimitiveVarTaylorHoodElement::EvaluateConvection(Vector &rResult,
                                         const array_1d<double, 3> &rConvVel,
                                         const ShapeDerivativesType &DNv_DX)
{
    const SizeType Dim = this->GetGeometry().WorkingSpaceDimension();
    const SizeType NumVNodes = this->GetGeometry().PointsNumber();

    if(rResult.size() != NumVNodes) rResult.resize(NumVNodes);

    for (SizeType i = 0; i < NumVNodes; i++)
    {
        rResult[i] = rConvVel[0]*DNv_DX(i,0);
        for(SizeType k = 1; k < Dim; k++)
            rResult[i] += rConvVel[k]*DNv_DX(i,k);
    }
}

}
