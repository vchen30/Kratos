// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix
//                   Klaus B. Sautter
//


// Project includes
#include "includes/define.h"
#include "custom_conditions/line_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************
    
    LineLoadCondition::LineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************
    
    LineLoadCondition::LineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************
    
    Condition::Pointer LineLoadCondition::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const
    {
        return boost::make_shared<LineLoadCondition>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************
    
    Condition::Pointer LineLoadCondition::Create( 
        IndexType NewId, 
        NodesArrayType const& ThisNodes,  
        PropertiesType::Pointer pProperties 
        ) const
    {
        return boost::make_shared<LineLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************
    
    LineLoadCondition::~LineLoadCondition()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void LineLoadCondition::CalculateAll( 
        MatrixType& rLeftHandSideMatrix, 
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag 
        )
    {
        KRATOS_TRY;
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int block_size = this->GetBlockSize();

        // Resizing as needed the LHS
        unsigned int mat_size = number_of_nodes * block_size;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != mat_size )
            {
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );
            }

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
        }


        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size( ) != mat_size )
            {
                rRightHandSideVector.resize( mat_size, false );
            }

            noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
        }


        IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(GetGeometry());
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(integration_method);
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(integration_method);

        


         // Vector with a loading applied to the elemnt
        array_1d<double, 3 > line_load = ZeroVector(3);
        if( this->Has( LINE_LOAD ) )
        {
            noalias(line_load) = this->GetValue( LINE_LOAD );
        }
        
        for ( unsigned int point_number = 0; point_number < integration_points.size(); point_number++ )
        {                   
            const double det_j = GetGeometry().DeterminantOfJacobian( integration_points[point_number] );

            const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j); 

            //set move_mesh_flag = false ! otherwise RHS for reaction forces uses new length in integration_weight
            
            for (unsigned int ii = 0; ii < number_of_nodes; ++ii)
            {
                unsigned int base = ii * block_size;
                
                for(unsigned int k = 0; k < dimension; ++k)
                {
                    rRightHandSideVector[base + k] += integration_weight * Ncontainer( point_number, ii ) * line_load[k];
                }
            }


        } 

       if (this->HasRotDof()) this->CalculateAndAddWorkEquivalentNodalForcesLineLoad(line_load,rRightHandSideVector);

        KRATOS_WATCH(rRightHandSideVector);
        KRATOS_CATCH( "" )
    }


    void LineLoadCondition::CalculateAndAddWorkEquivalentNodalForcesLineLoad(
        const Vector& ForceInput, VectorType& rRightHandSideVector) const
    {
        KRATOS_TRY;
        const int dimension = this->GetGeometry().WorkingSpaceDimension();
        //calculate orthogonal load vector
        Vector GeometricOrientation = ZeroVector(dimension);
        GeometricOrientation[0] = this->GetGeometry()[1].X() 
            - this->GetGeometry()[0].X();
        GeometricOrientation[1] = this->GetGeometry()[1].Y() 
            - this->GetGeometry()[0].Y();
        if (dimension == 3)
        {
            GeometricOrientation[2] = this->GetGeometry()[1].Z() 
                - this->GetGeometry()[0].Z();
        }

        const double VectorNormA = MathUtils<double>::Norm(GeometricOrientation);
        if (VectorNormA != 0.00) GeometricOrientation /= VectorNormA;

        Vector LineLoadDir = ZeroVector(dimension);
        for (int i = 0; i < dimension; ++i)
        {
            LineLoadDir[i] = ForceInput[i];
        }

        const double VectorNormB = MathUtils<double>::Norm(LineLoadDir);
        if (VectorNormB != 0.00) LineLoadDir /= VectorNormB;

        double cosAngle = 0.00;
        for (int i = 0; i < dimension; ++i)
        {
            cosAngle += LineLoadDir[i] * GeometricOrientation[i];
        }

        const double sinAngle = std::sqrt(1.00 - (cosAngle*cosAngle));
        const double NormForceVectorOrth = sinAngle * VectorNormB;


        Vector NodeA = ZeroVector(dimension);
        NodeA[0] = this->GetGeometry()[0].X();
        NodeA[1] = this->GetGeometry()[0].Y();
        if (dimension == 3)	NodeA[2] = this->GetGeometry()[0].Z();

        Vector NodeB = ZeroVector(dimension);
        NodeB = NodeA + LineLoadDir;

        Vector NodeC = ZeroVector(dimension);
        NodeC = NodeA + (GeometricOrientation*cosAngle);

        Vector LoadOrthogonalDir = ZeroVector(dimension);
        LoadOrthogonalDir = NodeB - NodeC;
        const double VectorNormC = MathUtils<double>::Norm(LoadOrthogonalDir);
        if(VectorNormC != 0.00) LoadOrthogonalDir /= VectorNormC;



        // now caluclate respective work equivilent nodal moments

        const double CustomMoment = (NormForceVectorOrth *
            VectorNormA*VectorNormA) / 12.00;

        Vector MomentNodeA = ZeroVector(dimension);
        MomentNodeA = MathUtils<double>::CrossProduct(GeometricOrientation,
            LoadOrthogonalDir);
        MomentNodeA *= CustomMoment;

        if(dimension == 3)
        {        
            for (int i = 0; i < dimension; ++i)
            {
                rRightHandSideVector[(1 * dimension) + i] += MomentNodeA[i];
                rRightHandSideVector[(3 * dimension) + i] -= MomentNodeA[i];
            }
        }
        else if(dimension == 2)
        {
            rRightHandSideVector[2] += MomentNodeA[2];
            rRightHandSideVector[5] -= MomentNodeA[2];
        }
        else KRATOS_ERROR << "the conditions only works for 2D and 3D elements";


        KRATOS_CATCH("")            
    }



} // Namespace Kratos