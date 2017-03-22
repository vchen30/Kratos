/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        KratosParticleMechanicsApplication $
//   Last modified by:    $Author:                 IIaconeta $
//   Date:                $Date:               November 2016 $
//   Revision:            $Revision:                     0.0 $
//
//


// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_UP2.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "particle_mechanics_application.h"


#include <omp.h>
#include <sstream>

namespace Kratos
{
    
///**
 //* Flags related to the element computation
 //*/
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_RHS_VECTOR,                 0 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_LHS_MATRIX,                 1 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_RHS_VECTOR_WITH_COMPONENTS, 2 );
//KRATOS_CREATE_LOCAL_FLAG( UpdatedLagrangian, COMPUTE_LHS_MATRIX_WITH_COMPONENTS, 3 );

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

	UpdatedLagrangianUP2::UpdatedLagrangianUP2( )
		: UpdatedLagrangian( )
	{
	  //DO NOT CALL IT: only needed for Register and Serialization!!!
	}
//******************************CONSTRUCTOR*******************************************
//************************************************************************************
    UpdatedLagrangianUP2::UpdatedLagrangianUP2( IndexType NewId, GeometryType::Pointer pGeometry )
            : UpdatedLagrangian( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

    UpdatedLagrangianUP2::UpdatedLagrangianUP2( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : UpdatedLagrangian( NewId, pGeometry, pProperties )
    {
    mFinalizedStep = true; 


    }
//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

    UpdatedLagrangianUP2::UpdatedLagrangianUP2( UpdatedLagrangianUP2 const& rOther)
        :UpdatedLagrangian(rOther)
        //,mDeformationGradientF0(rOther.mDeformationGradientF0)
        //,mDeterminantF0(rOther.mDeterminantF0)
    {
    }

//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

    UpdatedLagrangianUP2&  UpdatedLagrangianUP2::operator=(UpdatedLagrangianUP2 const& rOther)
    {
        //UpdatedLagrangian::operator=(rOther);

        //mDeformationGradientF0.clear();
        //mDeformationGradientF0 = rOther.mDeformationGradientF0;
    

        mDeterminantF0 = rOther.mDeterminantF0;
        

        return *this;
    }
//*********************************OPERATIONS*****************************************
//************************************************************************************

    Element::Pointer UpdatedLagrangianUP2::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new UpdatedLagrangianUP2( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }
//************************************CLONE*******************************************
//************************************************************************************

    Element::Pointer UpdatedLagrangianUP2::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
    {

        UpdatedLagrangianUP2 NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

        //-----------//

        
        


        NewElement.mConstitutiveLawVector = mConstitutiveLawVector->Clone();



        //-----------//


        
        NewElement.mDeformationGradientF0 = mDeformationGradientF0;
        
        
        NewElement.mDeterminantF0 = mDeterminantF0;
        
        return Element::Pointer( new UpdatedLagrangianUP2(NewElement) );
    }
//*******************************DESTRUCTOR*******************************************
//************************************************************************************
    UpdatedLagrangianUP2::~UpdatedLagrangianUP2()
    {
    }


//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::Initialize()
    {
        KRATOS_TRY
        
        
        array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
       
        

        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        
        //Constitutive Law initialisation

        //if ( mConstitutiveLawVector.size() != 1 )  
        //{
            //mConstitutiveLawVector.resize( 1 );
           
        //}
        
        mDeterminantF0 = 1;
        
        mDeformationGradientF0 = identity_matrix<double> (dim);
              
        
        //Compute jacobian inverses 
        
        Matrix J0 = ZeroMatrix(dim, dim);
        
        J0 = this->MPMJacobian(J0, xg);    
        
        //calculating and storing inverse and the determinant of the jacobian 
        MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );
        
        Matrix j = ZeroMatrix(dim,dim);
        j = this->MPMJacobian(j,xg);
        double detj;
        MathUtils<double>::InvertMatrix( j, mInverseJ, detj );
        
        InitializeMaterial();  
        
        
        //double MP_KineticEnergy = 0.0;
        //double MP_StrainEnergy = 0.0;
                           
        //for(unsigned int k = 0;k<3;k++)
            //{
            //MP_KineticEnergy += 0.5 * this->GetValue(MP_MASS) * this->GetValue(MP_VELOCITY)[k] * this->GetValue(MP_VELOCITY)[k] ;
            //}
        //for(unsigned int j = 0; j < this->GetValue(MP_CAUCHY_STRESS_VECTOR).size(); j++)
            //{
            //MP_StrainEnergy +=  0.5 * this->GetValue(MP_VOLUME) * this->GetValue(MP_CAUCHY_STRESS_VECTOR)[j] * this->GetValue(MP_ALMANSI_STRAIN_VECTOR)[j];
            //}
        
        //this->GetValue(MP_KINETIC_ENERGY) = MP_KineticEnergy;
        //this->GetValue(MP_STRAIN_ENERGY) = MP_StrainEnergy;
        this->GetValue(MP_DENSITY) = GetProperties()[DENSITY];
        

        //std::cout<<"The element is initialized"<<std::endl;


        KRATOS_CATCH( "" )
    }
//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

    UpdatedLagrangian::InitializeGeneralVariables(rVariables,rCurrentProcessInfo);

    //stabilization factor
    
    const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int voigtsize  = 3;
            
        if( dimension == 3 )
        {
            voigtsize  = 6;
        }
        rVariables.detF  = 1;

        rVariables.detF0 = 1;
        
        rVariables.detFT = 1;

        rVariables.detJ = 1;

        rVariables.B.resize( voigtsize , number_of_nodes * dimension );

        rVariables.F.resize( dimension, dimension );

        rVariables.F0.resize( dimension, dimension );
        
        rVariables.FT.resize( dimension, dimension );
        
        rVariables.StrainVector.resize( voigtsize );

        //I need to information in 3D
        
        rVariables.ConstitutiveMatrix.resize( 6, 6 );

        rVariables.StressVector.resize( 6 );

        rVariables.DN_DX.resize( number_of_nodes, dimension );
        rVariables.DN_De.resize( number_of_nodes, dimension );
        
        array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
        
        rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);
        
		//if (this->Id() == 17693)
		//{
			//double node1ID = this->GetGeometry()[0].Id();
			//double node2ID = this->GetGeometry()[1].Id();
			//double node3ID = this->GetGeometry()[2].Id();
			//double node4ID = this->GetGeometry()[3].Id();
			//std::cout<<" rVariables.N "<< rVariables.N<<std::endl;
			//std::cout<<" xg "<< xg<<std::endl;
			//std::cout<<" node1ID "<< node1ID<<std::endl;
			//std::cout<<" node2ID "<< node2ID<<std::endl;
			//std::cout<<" node3ID "<< node3ID<<std::endl;
			//std::cout<<" node4ID "<< node4ID<<std::endl;
			
		//}
		
		
        //reading shape functions local gradients
        rVariables.DN_De = this->MPMShapeFunctionsLocalGradients( rVariables.DN_De);
        
        
        //**********************************************************************************************************************
        
        //CurrentDisp is the variable unknown. It represents the nodal delta displacement. When it is predicted is equal to zero.
        
        rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
        
        
        //calculating the current jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n+1/d£]
        rVariables.j = this->MPMJacobianDelta( rVariables.j, xg, rVariables.CurrentDisp);
        
                
        //calculating the reference jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n/d£]
        rVariables.J = this->MPMJacobian( rVariables.J, xg);
        
        double StabilizationFactor = 1.0;
		if( GetProperties().Has(STABILIZATION_FACTOR) ){
		  StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
		}
		else if( rCurrentProcessInfo.Has(STABILIZATION_FACTOR) ){
		  StabilizationFactor = rCurrentProcessInfo[STABILIZATION_FACTOR];
		}
		GetProperties().SetValue(STABILIZATION_FACTOR, StabilizationFactor);
        
        //std::cout<<"The general variables are initialized"<<std::endl;
        //*************************************************************************************************************************
    
    


	
    KRATOS_CATCH( "" )

    }
//************************************************************************************
//************************************************************************************
	
    
 /**
     * The position of the Gauss points/Material points is updated
     */

    void UpdatedLagrangianUP2::UpdateGaussPoint( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        
        array_1d<double,3> xg = this->GetValue(GAUSS_COORD);
        array_1d<double,3> MP_PreviousAcceleration = this->GetValue(MP_ACCELERATION);
        array_1d<double,3> MP_PreviousVelocity = this->GetValue(MP_VELOCITY);
        //double MP_Mass = this->GetValue(MP_MASS);
        array_1d<double,3> delta_xg = ZeroVector(3);
        array_1d<double,3> MP_Acceleration = ZeroVector(3);
        array_1d<double,3> MP_Velocity = ZeroVector(3);
        double MP_Pressure = 0.0;
        //double DeltaTime = rCurrentProcessInfo[DELTA_TIME];
        
        
        rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);
        //int MP_number = this->GetValue(MP_NUMBER);
        
        //double total_nodal_mass = 0.0;
        //for ( unsigned int i = 0; i < number_of_nodes; i++ )
        //{
            //total_nodal_mass += GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);
        //}
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {   
            if (rVariables.N[i] > 1e-16)
            {
            array_1d<double, 3 > & NodalAcceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
            array_1d<double, 3 > & NodalVelocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            double NodalMass = GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0);
            array_1d<double,3> NodalMomentum = NodalMass * NodalVelocity;
            array_1d<double,3> NodalInertia = NodalMass * NodalAcceleration;
            
            //if (this->Id() == 8752)// || this->Id() == 1513)
            //{
                //std::cout<< "Nodal ID "<< GetGeometry()[i].Id()<<std::endl;
                //std::cout<< "NodalAcceleration "<<NodalAcceleration<<std::endl;
                //std::cout<< "NodalVelocity "<<NodalVelocity<<std::endl;
                //std::cout<< "NodalMass "<<NodalMass<<std::endl;
                
                //std::cout<< "rVariables.N "<<rVariables.N<<std::endl;
            //}
            double NodalPressure = GetGeometry()[i].GetSolutionStepValue(PRESSURE, 0);
            MP_Pressure += rVariables.N[i] * NodalPressure;
				
            
            
            for ( unsigned int j = 0; j < dimension; j++ )
                {   
                
                delta_xg[j] += rVariables.N[i] * rVariables.CurrentDisp(i,j);
                MP_Acceleration[j] += rVariables.N[i] * NodalAcceleration[j];
                MP_Velocity[j] += rVariables.N[i] * NodalVelocity[j];
                
                //MP_Acceleration[j] +=NodalInertia[j]/(rVariables.N[i] * MP_Mass * MP_number);//
                //MP_Velocity[j] += NodalMomentum[j]/(rVariables.N[i] * MP_Mass * MP_number);
                //MP_Velocity[j] += DeltaTime * rVariables.N[i] * NodalAcceleration[j];////
                
                
            
                
                }
            }
            
        }
     
        
        //**************************************************************************************************************************
        //Another way to update the MP velocity (see paper Guilkey and Weiss, 2003) 
        //MP_Velocity = MP_PreviousVelocity + 0.5 * DeltaTime * (MP_Acceleration + MP_PreviousAcceleration);
        //MP_Acceleration = 4/(DeltaTime * DeltaTime) * delta_xg - 4/DeltaTime * MP_PreviousVelocity;
        //MP_Velocity = 2/DeltaTime * delta_xg - MP_PreviousVelocity;
        
        
        this -> SetValue(MP_PRESSURE,MP_Pressure ); 
        
        
        //if(this->Id() == 13290)
		//{
                  
			//std::cout<<" MP_Pressure_FINALIZE " <<  MP_Pressure<<std::endl;
        //}

        this -> SetValue(MP_VELOCITY,MP_Velocity ); 
        
        
        const array_1d<double,3>& new_xg = xg + delta_xg ;   
        
        //Update the MP Position           
        this -> SetValue(GAUSS_COORD,new_xg);
        
        //Update the MP Acceleration
        this -> SetValue(MP_ACCELERATION,MP_Acceleration);
        
        array_1d<double,3>& MP_Displacement = this->GetValue(MP_DISPLACEMENT);
       
        MP_Displacement += delta_xg;
        
        //Update the MP Displacement
        this -> SetValue(MP_DISPLACEMENT,MP_Displacement );  
        
        
        
        //if (this->Id() == 1518 || this->Id() == 1513)
        
        //{
            //std::cout<<" MP position "<<this->Id()<<this -> GetValue(GAUSS_COORD)<<std::endl;
            //std::cout<<" delta_xg "<<this->Id()<<delta_xg<<std::endl;
            
            //std::cout<<" MP_Velocity "<<this->Id()<<this -> GetValue(MP_VELOCITY)<<std::endl;
            
            //std::cout<<" MP_Acceleration "<<this->Id()<<this -> GetValue(MP_ACCELERATION)<<std::endl;
            
        //}
        
        KRATOS_CATCH( "" )
    }
//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangian::SetGeneralVariables(GeneralVariables& rVariables,
            //ConstitutiveLaw::Parameters& rValues)
    //{   
        ////Variables.detF is the determinant of the incremental total deformation gradient
        //rVariables.detF  = MathUtils<double>::Det(rVariables.F);

        //if(rVariables.detF<0){
            
            //std::cout<<" Element: "<<this->Id()<<std::endl;
            //std::cout<<" Element position "<<this->GetValue(GAUSS_COORD)<<std::endl;
            //unsigned int number_of_nodes = GetGeometry().PointsNumber();

            //for ( unsigned int i = 0; i < number_of_nodes; i++ )
            //{
            //array_1d<double, 3> &CurrentPosition  = GetGeometry()[i].Coordinates();
            //array_1d<double, 3> & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            //array_1d<double, 3> & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            ////array_1d<double, 3> PreviousPosition  = CurrentPosition - (CurrentDisplacement-PreviousDisplacement);         
            //std::cout<<" NODE ["<<GetGeometry()[i].Id()<<"]: (Current position: "<<CurrentPosition<<") "<<std::endl;
            //std::cout<<" ---Current Disp: "<<CurrentDisplacement<<" (Previour Disp: "<<PreviousDisplacement<<")"<<std::endl;
            //}

            //for ( unsigned int i = 0; i < number_of_nodes; i++ )
            //{
                //if( GetGeometry()[i].SolutionStepsDataHas(CONTACT_FORCE) ){
                //array_1d<double, 3 > & PreContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE,1);
                //array_1d<double, 3 > & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
                //std::cout<<" ---Contact_Force: (Pre:"<<PreContactForce<<", Cur:"<<ContactForce<<") "<<std::endl;
                //}
                //else{
                //std::cout<<" ---Contact_Force: NULL "<<std::endl;
                //}   
            //}
        
            //KRATOS_THROW_ERROR( std::invalid_argument," MPM UPDATED LAGRANGIAN DISPLACEMENT ELEMENT INVERTED: |F|<0  detF = ", rVariables.detF )
        //}
        
        //rVariables.detFT = rVariables.detF * rVariables.detF0;
        //rVariables.FT    = prod( rVariables.F, rVariables.F0 );

        
        //rValues.SetDeterminantF(rVariables.detFT);
        //rValues.SetDeformationGradientF(rVariables.FT);
        //rValues.SetStrainVector(rVariables.StrainVector);
        //rValues.SetStressVector(rVariables.StressVector);
        //rValues.SetConstitutiveMatrix(rVariables.ConstitutiveMatrix);
        //rValues.SetShapeFunctionsDerivatives(rVariables.DN_DX);
        //rValues.SetShapeFunctionsValues(rVariables.N);
        
        
        ////std::cout<<"The general variables are set"<<std::endl;

    //}
//************************************************************************************
//*****************check size of LHS and RHS matrices*********************************

    void UpdatedLagrangianUP2::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector,
            Flags& rCalculationFlags)

    {

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

    if ( rCalculationFlags.Is(UpdatedLagrangian::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(UpdatedLagrangian::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        
        //std::cout<<" AAAAAAAAAAAAAAAAAAAAAAAAAa "<<std::endl;
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
        //std::cout<<" rRightHandSideVector "<<rRightHandSideVector<<std::endl;
    }
        //std::cout<<"The system matrices are initialized"<<std::endl;
    }

//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangian::CalculateElementalSystem( LocalSystemComponents& rLocalSystem,
                                 //ProcessInfo& rCurrentProcessInfo)
    //{
        //KRATOS_TRY
        
        ////create and initialize element variables:
        //GeneralVariables Variables;
        
        //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);
        

        ////create constitutive law parameters:
        //ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        
       
        ////set constitutive law flags:
        //Flags &ConstitutiveLawOptions=Values.GetOptions();
        
        ////std::cout<<"in CalculateElementalSystem 5"<<std::endl;
        //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRAIN);
        
        //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
        
        //ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        
        
        ////auxiliary terms
        //Vector VolumeForce;
        
        
        ////compute element kinematics B, F, DN_DX ...
        //this->CalculateKinematics(Variables,rCurrentProcessInfo);
        
        ////set general variables to constitutivelaw parameters
        //this->SetGeneralVariables(Variables,Values);
        
        //mConstitutiveLawVector->CalculateMaterialResponse(Values, Variables.StressMeasure);
        
        ////this->SetValue(MP_CAUCHY_STRESS_VECTOR, Variables.StressVector);
        ////std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;
        ////std::cout<<"Variables.StressVector in the element "<<Variables.StressVector<<std::endl;
        
        ////this->SetValue(MP_ALMANSI_STRAIN_VECTOR, Variables.StrainVector);
        ////double EquivalentPlasticStrain = mConstitutiveLawVector->GetValue( PLASTIC_STRAIN, EquivalentPlasticStrain );
        ////this->SetValue(MP_EQUIVALENT_PLASTIC_STRAIN, EquivalentPlasticStrain);
        ////at the first iteration I recover the previous state of stress and strain
        //if(rCurrentProcessInfo[NL_ITERATION_NUMBER] == 1)
        //{
            //this->SetValue(PREVIOUS_MP_CAUCHY_STRESS_VECTOR, Variables.StressVector);
            //this->SetValue(PREVIOUS_MP_ALMANSI_STRAIN_VECTOR, Variables.StrainVector);
        //}
        ////the MP density is updated
        //double MP_Density = (GetProperties()[DENSITY]) / Variables.detFT;
        ////if(this->Id() == 1786 || this->Id() == 1836)
            ////{
                ////std::cout<<"density "<<this->Id() << GetProperties()[DENSITY]<<std::endl;
            ////}
        
        ////the integration weight is evaluated
        //double MP_Volume = this->GetValue(MP_MASS)/this->GetValue(MP_DENSITY);
        
        //this->SetValue(MP_DENSITY, MP_Density);
        //this->SetValue(MP_VOLUME, MP_Volume);
        
        
        
        
        //if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
        //{
        
        ////contributions to stiffness matrix calculated on the reference config
        //this->CalculateAndAddLHS ( rLocalSystem, Variables, MP_Volume );
        
        //}

        //if ( rLocalSystem.CalculationFlags.Is(UpdatedLagrangian::COMPUTE_RHS_VECTOR) ) //calculation of the vector is required
        //{
        ////contribution to external forces
        //VolumeForce  = this->CalculateVolumeForce( VolumeForce, Variables );
        
        //this->CalculateAndAddRHS ( rLocalSystem, Variables, VolumeForce, MP_Volume );
        
        //}

        
               
        //KRATOS_CATCH( "" )
    //}
//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


    void UpdatedLagrangianUP2::CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)

    {
        KRATOS_TRY
		
        //const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        
        //Define the stress measure
        rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

        //Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
        Matrix InvJ;
        
        MathUtils<double>::InvertMatrix( rVariables.J, InvJ, rVariables.detJ);
        
        
        
        //Calculating the inverse of the jacobian and the parameters needed [d£/(dx_n+1)]
        Matrix Invj;
        MathUtils<double>::InvertMatrix( rVariables.j, Invj, rVariables.detJ ); //overwrites detJ

        //Compute cartesian derivatives [dN/dx_n+1]        
        rVariables.DN_DX = prod( rVariables.DN_De, Invj); //overwrites DX now is the current position dx
        
        
        
        //Deformation Gradient F [(dx_n+1 - dx_n)/dx_n] to be updated in constitutive law parameter as total deformation gradient
        //the increment of total deformation gradient can be evaluated in 2 ways.
        //1 way.
        noalias( rVariables.F ) = prod( rVariables.j, InvJ);
        
        //2 way by means of the gradient of nodal displacement: using this second expression quadratic convergence is not guarantee

        //Matrix I=identity_matrix<double>( dimension );
        
        //Matrix GradientDisp = ZeroMatrix(dimension, dimension);
        //rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
        //GradientDisp = prod(trans(rVariables.CurrentDisp),rVariables.DN_DX);
        
        ////REMEMBER THAT USING JUST ONLY THE FIRST ORDER TERM SOME ISSUES CAN COME UP WHEN FOR PROBLEMS WITH LOTS OF ROTATIONAL MOTION(slender cantilever beam??)
        //noalias( rVariables.F ) = (I + GradientDisp);
        //if (this->Id() == 365)
        //{
            //std::cout<<" AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "<<std::endl;
            //std::cout<<"rVariables.CurrentDisp in calculate kinematic "<<rVariables.CurrentDisp<<std::endl;
            //std::cout<<" AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "<<std::endl;
        //}
        
        
        
        
        
        //Determinant of the Deformation Gradient F_n
                        
        rVariables.detF0 = mDeterminantF0;
        rVariables.F0    = mDeformationGradientF0;
        
        //if(this->Id() == 365)
        //{   
        
        //std::cout<<"rVariables.DN_DX "<<this->Id()<<rVariables.DN_DX<<std::endl;
        //std::cout<<"rVariables.DN_De "<<this->Id()<<rVariables.DN_De<<std::endl;
        //std::cout<<"rVariables.J "<<this->Id()<<rVariables.J<<std::endl;
        //std::cout<<"rVariables.j "<<this->Id()<<rVariables.j<<std::endl;
        //std::cout<<"Invj "<<this->Id()<<Invj<<std::endl;
        //}
        
        //Compute the deformation matrix B
        this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX);

        
        KRATOS_CATCH( "" )
    }
//************************************************************************************

    void UpdatedLagrangianUP2::CalculateDeformationMatrix(Matrix& rB,
            Matrix& rF,
            Matrix& rDN_DX)
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
        
        rB.clear(); //set all components to zero

        if( dimension == 2 )
        {

            for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {
                unsigned int index = 2 * i;

                rB( 0, index + 0 ) = rDN_DX( i, 0 );
                rB( 1, index + 1 ) = rDN_DX( i, 1 );
                rB( 2, index + 0 ) = rDN_DX( i, 1 );
                rB( 2, index + 1 ) = rDN_DX( i, 0 );

            }
            //if(this->Id() == 365)
        //{
			//std::cout<<"rB "<< this->Id()<< rB<<std::endl;
		//}

        }
        
        else if( dimension == 3 )
        {

            for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {
                unsigned int index = 3 * i;

                rB( 0, index + 0 ) = rDN_DX( i, 0 );
                rB( 1, index + 1 ) = rDN_DX( i, 1 );
                rB( 2, index + 2 ) = rDN_DX( i, 2 );

                rB( 3, index + 0 ) = rDN_DX( i, 1 );
                rB( 3, index + 1 ) = rDN_DX( i, 0 );

                rB( 4, index + 1 ) = rDN_DX( i, 2 );
                rB( 4, index + 2 ) = rDN_DX( i, 1 );

                rB( 5, index + 0 ) = rDN_DX( i, 2 );
                rB( 5, index + 2 ) = rDN_DX( i, 0 );

            }
        }
        else
        {

            KRATOS_THROW_ERROR( std::invalid_argument, "something is wrong with the dimension", "" )

        }

        KRATOS_CATCH( "" )
    }
    
////************************************************************************************
////************************************************************************************

    void UpdatedLagrangianUP2::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
    {   
        // In the Initialize of each time step the nodal initial conditions are evaluated
        //1. first of all I need to evaluate the MP momentum and MP_inertia
                
        
        
        //int MP_bool = this->GetValue(MP_BOOL);
        
        //std::cout<<" in InitializeSolutionStep2"<<std::endl;
            unsigned int dimension = GetGeometry().WorkingSpaceDimension();
            const unsigned int number_of_nodes = GetGeometry().PointsNumber();
            array_1d<double,3>& xg = this->GetValue(GAUSS_COORD);
            GeneralVariables Variables;
            //this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);
            
            
            
            Matrix J0 = ZeroMatrix(dimension, dimension);
            
            J0 = this->MPMJacobian(J0, xg);    
            
            //calculating and storing inverse and the determinant of the jacobian 
            MathUtils<double>::InvertMatrix( J0, mInverseJ0, mDeterminantJ0 );
            
            Variables.N = this->MPMShapeFunctionPointValues(Variables.N, xg);
            
            mConstitutiveLawVector->InitializeSolutionStep( GetProperties(),
                    GetGeometry(), Variables.N, rCurrentProcessInfo );

            mFinalizedStep = false;
            
            
            
            array_1d<double,3>& MP_Velocity = this->GetValue(MP_VELOCITY);
            array_1d<double,3>& MP_Acceleration = this->GetValue(MP_ACCELERATION);
            double MP_Pressure = this->GetValue(MP_PRESSURE);
            array_1d<double,3>& AUX_MP_Velocity = this->GetValue(AUX_MP_VELOCITY);
            array_1d<double,3>& AUX_MP_Acceleration = this->GetValue(AUX_MP_ACCELERATION);
            double AUX_MP_Pressure = this->GetValue(AUX_MP_PRESSURE);
            double MP_Mass = this->GetValue(MP_MASS);
            array_1d<double,3> MP_Momentum;
            array_1d<double,3> MP_Inertia;
            //double MP_MPressure;
            array_1d<double,3> NodalMomentum;
            array_1d<double,3> NodalInertia;
            double NodalMPressure;
            
           for (unsigned int j=0;j<number_of_nodes;j++)
            {   
                //these are the values of nodal velocity and nodal acceleration evaluated in the initialize solution step
                array_1d<double, 3 > & NodalAcceleration = GetGeometry()[j].FastGetSolutionStepValue(ACCELERATION,1);
                array_1d<double, 3 > & NodalVelocity = GetGeometry()[j].FastGetSolutionStepValue(VELOCITY,1);
                double & NodalPressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE,1);
                //if(this->Id() == 13290)
                //{
                  
                    //std::cout<<" AUX_MP_Pressure_BEFORE " <<  AUX_MP_Pressure<<std::endl;
                //}
                AUX_MP_Pressure += Variables.N[j] * NodalPressure;
                                
                //std::cout<<"NodalVelocity "<< GetGeometry()[j].Id()<<std::endl;
                for (unsigned int k = 0; k < dimension; k++)
                {
                AUX_MP_Velocity[k] += Variables.N[j] * NodalVelocity[k];
                AUX_MP_Acceleration[k] += Variables.N[j] * NodalAcceleration[k];
                }
            }
              //if(this->Id() == 13290)
                //{
                    
                    ////std::cout<<" MP_Velocity " <<  MP_Velocity<<std::endl;
                    ////std::cout<<" MP_Acceleration " <<  MP_Acceleration<<std::endl;
                    ////std::cout<<" AUX_MP_Velocity " <<  AUX_MP_Velocity<<std::endl;
                    ////std::cout<<" AUX_MP_Acceleration " <<  AUX_MP_Acceleration<<std::endl;
                    //std::cout<<" AUX_MP_Pressure_AFTER " <<  AUX_MP_Pressure<<std::endl;
                //}
           
            
            // Here MP contribution in terms of momentum, inertia and mass are added        
            //if(this->Id() == 13290)
                //{
                  
                    //std::cout<<" MP_Pressure " <<  MP_Pressure<<std::endl;
                //}
            for ( unsigned int i = 0; i < number_of_nodes; i++ )
            {  
				
				NodalMPressure =  Variables.N[i] * (MP_Pressure - AUX_MP_Pressure) * MP_Mass;
				
                for (unsigned int j = 0; j < dimension; j++)
                {
                    NodalMomentum[j] = Variables.N[i] * (MP_Velocity[j] - AUX_MP_Velocity[j]) * MP_Mass;
                    NodalInertia[j] = Variables.N[i] * (MP_Acceleration[j] - AUX_MP_Acceleration[j]) * MP_Mass;
                    
                } 
                GetGeometry()[i].GetSolutionStepValue(NODAL_MOMENTUM, 0) += NodalMomentum;
                GetGeometry()[i].GetSolutionStepValue(NODAL_INERTIA, 0) += NodalInertia;
                GetGeometry()[i].GetSolutionStepValue(NODAL_MPRESSURE, 0) += NodalMPressure;
                GetGeometry()[i].GetSolutionStepValue(NODAL_MASS, 0) += Variables.N[i] * MP_Mass;
                
            }
            
            AUX_MP_Velocity.clear();
            AUX_MP_Acceleration.clear();
            AUX_MP_Pressure = 0.0;
        
        
         
        
        
    }
//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::CalculateAndAddRHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, Vector& rVolumeForce, double& rIntegrationWeight)
    {

	const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;

      mElementScalingNumber = 1.0;

      if (YoungModulus > 1e-5 ) {
         mElementScalingNumber = 100.0/BulkModulus; 
      }
      else {
         mElementScalingNumber = 0.01;
      }

        //contribution of the internal and external forces
        //contribution of the internal and external forces
    VectorType& rRightHandSideVector = rLocalSystem.GetRightHandSideVector(); 
      
    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian  
    
    ThisElementGeneralVariables ElementVariables;
    CalculateThisElementGeneralVariables( ElementVariables, rVariables);
      
    // operation performed: rRightHandSideVector += ExtForce*IntegrationWeight
    CalculateAndAddExternalForces( rRightHandSideVector, rVariables, rVolumeForce, rIntegrationWeight );

    // operation performed: rRightHandSideVector -= IntForce*IntegrationWeight
    CalculateAndAddInternalForces( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

    // operation performed: rRightHandSideVector -= PressureForceBalance*IntegrationWeight
    CalculateAndAddPressureForces( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);

    // operation performed: rRightHandSideVector -= Stabilized Pressure Forces
    CalculateAndAddStabilizedPressure( rRightHandSideVector, rVariables, ElementVariables, rIntegrationWeight);
    
    //std::cout<<" rRightHandSideVector " << rRightHandSideVector<<std::endl;
    
    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;
    
        
    }
//************************************************************************************
//*********************Calculate the contribution of external force*******************

    void UpdatedLagrangianUP2::CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
            GeneralVariables& rVariables,
            Vector& rVolumeForce,
            double& rIntegrationWeight)

    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        
        

        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            int indexup = dimension * i + i;
            
            for ( unsigned int j = 0; j < dimension; j++ )
            {
            rRightHandSideVector[indexup + j] += rVariables.N[i] * rVolumeForce[j]; 
            
            }
            
        }
        
        
        
        KRATOS_CATCH( "" )
    }
//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
            GeneralVariables & rVariables,
            ThisElementGeneralVariables& rElementVariables, 
            double& rIntegrationWeight)
    {
        KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		Vector StressVector = rElementVariables.StressVector;
		VectorType Fh=rRightHandSideVector;
      //Vector InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), StressVector );
		
		//in the largedisplacement UP element is defined Vector
        VectorType InternalForces = rIntegrationWeight * prod( trans( rVariables.B ), StressVector );
        
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
       {
        unsigned int indexup = dimension * i + i;
        unsigned int indexu  = dimension * i;

			for ( unsigned int j = 0; j < dimension; j++ )
			{
				rRightHandSideVector[indexup + j] -= InternalForces[indexu + j];
			}
        }
        //std::cout<< "StressVector "<<this->Id()<<" "<<StressVector<<std::endl;
        //std::cout<<"rRightHandSideVector "<<this->Id()<<" "<<rRightHandSideVector<<std::endl;
      //if(this->Id() == 8034 || this->Id() == 8035|| this->Id() == 9045|| this->Id() == 9061|| this->Id() == 9062)
        //{
			//std::cout<< "StressVector "<<this->Id()<<" "<<StressVector<<std::endl;
			//std::cout<<"rRightHandSideVector "<<this->Id()<<" "<<rRightHandSideVector<<std::endl;
		//}         
		//if(this->Id() == 8183)
		//{
		//std::cout<<" FInt "<<rRightHandSideVector-Fh<<std::endl;
		//}

        KRATOS_CATCH( "" )
    }
    
//******************************************************************************************************************
//******************************************************************************************************************
	double& UpdatedLagrangianUP2::CalculatePUCoefficient(double& rCoefficient, GeneralVariables & rVariables)
	{
	  KRATOS_TRY
	   
		//Mechanical volumetric:
		
		//Constitutive A:
		//rCoefficient = 0.5*(rVariables.detF0*rVariables.detF0-1)/rVariables.detF0); //(J²-1)/2
		
		//Constitutive B:
		rCoefficient = (std::log(rVariables.detF0)/rVariables.detF0);  //(ln(J)) for a free energy function of 0.5K ln^2(J)
	  
		//Thermal volumetric:

		
		return rCoefficient;

	  KRATOS_CATCH( "" )
	}


//************************************************************************************
//************************************************************************************

	double& UpdatedLagrangianUP2::CalculatePUDeltaCoefficient(double &rDeltaCoefficient, GeneralVariables & rVariables)
	{

	  KRATOS_TRY

		//Mechanical volumetric:

		//Constitutive A:
		//rDeltaCoefficient = (rVariables.detF0*rVariables.detF0 + 1)/(rVariables.detF0*rVariables.detF0); //(J²-1)/2

		//Constitutive B:
		rDeltaCoefficient = (1.0-std::log(rVariables.detF0))/(rVariables.detF0*rVariables.detF0);   //(ln(J)) for a free energy function of 0.5K ln^2(J)

		
		return rDeltaCoefficient;


	    KRATOS_CATCH( "" )

	}   
 
//************************************************************************************
//************************************************************************************
// I changed the term of the pressure: as after a first prediction of the pressure
// variable, this term is evaluated in the return mapping of the DP plastic model.
// After plasticity the pressure is evaluated directly on the integration point
// here called "NewPressure"
//************************************************************************************
//************************************************************************************

	void UpdatedLagrangianUP2::CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
									   GeneralVariables & rVariables,
									   ThisElementGeneralVariables& rElementVariables, 
									   double& rIntegrationWeight)
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		unsigned int indexp = dimension;

		VectorType Fh=rRightHandSideVector;

		

		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			for ( unsigned int j = 0; j < number_of_nodes; j++ )
			{

				double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
               
			//*********************************************************
			//mixed formulation considering a variation of p: 1st term
			rRightHandSideVector[indexp] += rVariables.N[i] * rVariables.N[j] * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF) * mElementScalingNumber;
			//*********************************************************
			}
			
			//rRightHandSideVector[indexp] -=  Coefficient * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);
			
			//********************************************************
			//mixed formulation considering a variation of p: 2nd term
			rRightHandSideVector[indexp] -=  rVariables.N[i] * rElementVariables.ElementalMeanStress * rIntegrationWeight / (rVariables.detF0/rVariables.detF)* mElementScalingNumber;
			//********************************************************
			indexp += (dimension + 1);
			
		}
		
		//if(this->Id() == 8183)
		//{
		//std::cout<<" Fpres "<<rRightHandSideVector-Fh<<std::endl;
		//}

		// std::cout<<std::endl;
		// std::cout<<" Coefficient " <<Coefficient<<" F0 "<<rVariables.detF0<<std::endl;
		// std::cout<<" Fpres "<<rRightHandSideVector-Fh<<std::endl;

		KRATOS_CATCH( "" )
	}
//************************************************************************************
//************************************************************************************

 	void UpdatedLagrangianUP2::CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
			GeneralVariables & rVariables,
			ThisElementGeneralVariables& rElementVariables, 
			double& rIntegrationWeight)
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		unsigned int indexp = dimension;

		VectorType Fh=rRightHandSideVector;
		// std::cout<<" Element "<<this->Id()<<" "<<std::endl;

		//use of this variable for the complete parameter:
		double AlphaStabilization  = 8.0; 
		double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
		AlphaStabilization *= StabilizationFactor;
		
		const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
		const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

		double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
		double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;
		
		AlphaStabilization=(AlphaStabilization/(36.0*LameMu));
		AlphaStabilization *= BulkModulus;
		
		//Experimental
		// if(LameMu < rVariables.ConstitutiveMatrix(2,2))
		//   LameMu = rVariables.ConstitutiveMatrix(2,2);

		//double NewPressure = mConstitutiveLawVector->GetValue(MP_PRESSURE, NewPressure );
		double consistent = 1;
		//double FactorValue = 8.0; //JMR deffault value	 
		//if( dimension == 3 )
		  //FactorValue = 10.0; //JMC deffault value
		  
		  
		//NON FUNZIONA  
		//consistent = AlphaStabilization*FactorValue / LameMu;
		////double stabilization_term = 1;

		//for ( unsigned int i = 0; i < number_of_nodes; i++ )
		//{
			////for ( unsigned int j = 0; j < number_of_nodes; j++ )
			////{

				////double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
				  

				////rRightHandSideVector[indexp] += consistent * ( -1/9)* Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D
			////}
				//rRightHandSideVector[indexp] += consistent * (rVariables.N[i]*rVariables.N[j] -1/9)* NewPressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D
		


			//indexp += (dimension + 1);
		//}
		for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {

            double& Pressure = GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
			
	    if( dimension == 2 ){ //consistent 2D

	      consistent=(-1)*AlphaStabilization;//*FactorValue/(36.0*LameMu);
	      if(i==j)
                consistent=2*AlphaStabilization;//*FactorValue/(36.0*LameMu);

	      //rRightHandSideVector[indexp] += consistent * (rVariables.N[i]*rVariables.N[j] -1/9)* Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D
	      rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF)* mElementScalingNumber; //2D
	      
	    }
	    else{

	      consistent=(-1)*AlphaStabilization;//*FactorValue/(80.0*LameMu);
	      if(i==j)
                consistent=3*AlphaStabilization;//*FactorValue/(80.0*LameMu);

	      rRightHandSideVector[indexp] += consistent * Pressure * rIntegrationWeight / (rVariables.detF0/rVariables.detF);//* mElementScalingNumber; //3D

	    }

	    // std::cout<<" Pressure "<<Pressure<<std::endl;
        }


        indexp += (dimension + 1);
    }

//********************************************************************************************************************************************
		
	

		//if(this->Id() == 8183)
		//{
		//std::cout<<" FpStab "<<rRightHandSideVector-Fh<<std::endl;
		//}




		// std::cout<<std::endl;
		// std::cout<<" IntegrationWeight "<<rIntegrationWeight<<" detF "<<rVariables.detF0<<std::endl;
		// std::cout<<" FpStab "<<rRightHandSideVector-Fh<<std::endl;

		KRATOS_CATCH( "" )
	}
//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::CalculateAndAddLHS(LocalSystemComponents& rLocalSystem, GeneralVariables& rVariables, double& rIntegrationWeight)
    {

	// Scaling constant for the second equation.
      const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
      const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

      double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;

      mElementScalingNumber = 1.0;

      if (YoungModulus > 1e-5 ) {
         mElementScalingNumber = 100.0/BulkModulus; 
      }
      else {
         mElementScalingNumber = 0.01;
      }


      //contributions of the stiffness matrix calculated on the reference configuration
    //contributions of the stiffness matrix calculated on the reference configuration
    MatrixType& rLeftHandSideMatrix = rLocalSystem.GetLeftHandSideMatrix();
    rVariables.detF0   *= rVariables.detF;
    double DeterminantF = rVariables.detF;
    rVariables.detF = 1; //in order to simplify updated and spatial lagrangian
    // operation performed: add Km to the rLefsHandSideMatrix


	ThisElementGeneralVariables ElementVariables;
    this->CalculateThisElementGeneralVariables( ElementVariables, rVariables);
    MatrixType Kh;//=rLeftHandSideMatrix;
    //respect to the current configuration n+1
    CalculateAndAddKuum( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );
    
    //std::cout<<" Kuum " << rLeftHandSideMatrix<<std::endl;
    
    Kh=rLeftHandSideMatrix;
    //std::cout<<"111111111111111111111"<<std::endl;
    // operation performed: add Kg to the rLefsHandSideMatrix
    CalculateAndAddKuug( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );
    
    //std::cout<<" Kuug " << rLeftHandSideMatrix - Kh<<std::endl;
    Kh=rLeftHandSideMatrix;
	//std::cout<<"222222222222222222222"<<std::endl;
    // operation performed: add Kup to the rLefsHandSideMatrix
    CalculateAndAddKup( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );
    //std::cout<<" Kup " << rLeftHandSideMatrix - Kh<<std::endl;
    
    Kh=rLeftHandSideMatrix;
	//std::cout<<"33333333333333333333"<<std::endl;
    // operation performed: add Kpu to the rLefsHandSideMatrix
    CalculateAndAddKpu( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight);
    //std::cout<<" Kpu " << rLeftHandSideMatrix - Kh<<std::endl;
    Kh=rLeftHandSideMatrix;
	//std::cout<<"444444444444444444444"<<std::endl;
    // operation performed: add Kpp to the rLefsHandSideMatrix
    CalculateAndAddKpp( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );
    //std::cout<<" Kpp " << rLeftHandSideMatrix - Kh<<std::endl;
    Kh=rLeftHandSideMatrix;
	//std::cout<<"55555555555555555555"<<std::endl;
    // operation performed: add Kpp Stab to the rLefsHandSideMatrix
    CalculateAndAddKppStab( rLeftHandSideMatrix, rVariables, ElementVariables, rIntegrationWeight );
    //std::cout<<" KppStab " << rLeftHandSideMatrix - Kh<<std::endl;
    
	//std::cout<<"66666666666666666666"<<std::endl;
	
	//std::cout<<" rLeftHandSideMatrix " << rLeftHandSideMatrix<<std::endl;
    rVariables.detF     = DeterminantF;
    rVariables.detF0   /= rVariables.detF;
    
    //KRATOS_WATCH( rLeftHandSideMatrix )

        
    }
//************************************************************************************
//************************************************************************************

    //void UpdatedLagrangianUP2::CalculateAndAddKuum(MatrixType& rLeftHandSideMatrix,
            //GeneralVariables& rVariables,
            //ThisElementGeneralVariables& rElementVariables,
            //double& rIntegrationWeight
                                                      //)
    //{
        //KRATOS_TRY
        ////std::stringstream ss;
        
        ////unsigned int number_of_nodes = GetGeometry().size();
        ////unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        ////unsigned int voigtsize  = 3;
        ////unsigned int MatSize = number_of_nodes * dimension;
       
        ////Matrix temp = ZeroMatrix(voigtsize, MatSize);
        
        
        ////temp = prod( rVariables.ConstitutiveMatrix, rVariables.B );
        
        
               
        //Matrix Kuum = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( rVariables.ConstitutiveMatrix, rVariables.B ) ) ); 
        ////assemble into rk the material uu contribution:
		//const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		//unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		 //// MatrixType Kh=rLeftHandSideMatrix;

		//unsigned int indexi = 0;
		//unsigned int indexj = 0;
		//for ( unsigned int i = 0; i < number_of_nodes; i++ )
		//{
			//for ( unsigned int idim = 0; idim < dimension ; idim ++)
			//{
				//indexj=0;
				//for ( unsigned int j = 0; j < number_of_nodes; j++ )
				//{
					//for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
					//{
						//rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuum(indexi,indexj);
						//indexj++;
					//}
				//}
				//indexi++;
			//}
		//}

        ////std::cout << ss.str();

        //KRATOS_CATCH( "" )
    //}

   void UpdatedLagrangianUP2::CalculateAndAddKuum ( MatrixType& rLeftHandSideMatrix,
         GeneralVariables & rVariables,
         ThisElementGeneralVariables& rElementVariables, 
         double& rIntegrationWeight)
   {
      KRATOS_TRY

      //assemble into rk the material uu contribution:
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();
	  MatrixType Kh=rLeftHandSideMatrix;

      Matrix ECConstitutiveMatrix = rVariables.ConstitutiveMatrix;
      Matrix ConstitutiveMatrix = ZeroMatrix(rElementVariables.voigtsize,rElementVariables.voigtsize); 
      
      ECConstitutiveMatrix = prod( rElementVariables.DeviatoricTensor, ECConstitutiveMatrix);
	  //if(this->Id() == 19289)
        //{
			//std::cout<<" ECConstitutiveMatrix 0 "<<ECConstitutiveMatrix<<std::endl;
			
		//} 
      //ECConstitutiveMatrix = prod( rElementVariables.DeviatoricTensor, ECConstitutiveMatrix);
	  //double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

      //for (unsigned int i = 0; i < 3; ++i) {
        //for (unsigned int j = 0; j < 3 ; ++j) {
           //ECConstitutiveMatrix(i,j)  -= BulkModulus;
         //}
      //}
	  
	  //std::cout<<" DeviatoricTensor "<<rElementVariables.DeviatoricTensor<<std::endl;
	  //std::cout<<" ConstitutiveMatrix in Kuum1 "<<rVariables.ConstitutiveMatrix<<std::endl;
	  //std::cout<<" ECConstitutiveMatrix 1 "<<ECConstitutiveMatrix<<std::endl; 
	  //std::cout<<" deviatoric constitutive matrix "<<ECConstitutiveMatrix<<std::endl;    


      // a. 2 ( pElem - p Nodal) I4S
      //for ( unsigned int i = 0; i < rElementVariables.voigtsize; i++) {
         
         //ECConstitutiveMatrix(i,i) += 2.0 * ( rElementVariables.ElementalMeanStress - rElementVariables.NodalMeanStress);
      //}
      for ( unsigned int i = 0; i < 6; i++) {
         double voigtnumber = 1.0;
         if ( i > 2)
            voigtnumber = 0.500;
         ECConstitutiveMatrix(i,i) += 2.0*voigtnumber * ( rElementVariables.ElementalMeanStress - rElementVariables.NodalMeanStress);
      }
      
      //if(this->Id() == 19289)
        //{
			//std::cout<<" pg-p "<<rElementVariables.ElementalMeanStress - rElementVariables.NodalMeanStress<<std::endl;
			//std::cout<<" ECConstitutiveMatrix 1 "<<ECConstitutiveMatrix<<std::endl;
			
		//} 
      //if(this->Id() == 9794)
        //{
			//std::cout<<" ECConstitutiveMatrix 2"<<ECConstitutiveMatrix<<std::endl;
			
		//} 
      //std::cout<<"2( pElem - p Nodal) "<<2.0 * ( rElementVariables.ElementalMeanStress - rElementVariables.NodalMeanStress)<<std::endl;
      //std::cout<<" ECConstitutiveMatrix 2 "<<ECConstitutiveMatrix<<std::endl;    

      // b. p 1 cross 1
      //if(rElementVariables.voigtsize==3){
	  //double alfa = 1/2.0;
      //for (unsigned int i = 0; i < 2; i++) {
         //for (unsigned int j = 0; j < 2; j++) {
            //ECConstitutiveMatrix(i,j) += (rElementVariables.NodalMeanStress - 2 * alfa * rElementVariables.ElementalMeanStress) ;
         //}
      //}
	  //}
	  //else if(rElementVariables.voigtsize==6){
	  //double alfa = 1/3.0;
      //for (unsigned int i = 0; i < 3; i++) {
         //for (unsigned int j = 0; j < 3; j++) {
            //ECConstitutiveMatrix(i,j) += (rElementVariables.NodalMeanStress- 2 * alfa * rElementVariables.ElementalMeanStress) ;
         //}
      //}
	  //}
	  double alfa = 1/3.0;
	  for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
            ECConstitutiveMatrix(i,j) += (rElementVariables.NodalMeanStress- 2 * alfa * rElementVariables.ElementalMeanStress) ;
         }
      }
      //if(this->Id() == 19289)
        //{
			//std::cout<<" pg-2 alfa p  "<<rElementVariables.NodalMeanStress- 2 * alfa * rElementVariables.ElementalMeanStress<<std::endl;
			//std::cout<<" ECConstitutiveMatrix 2 "<<ECConstitutiveMatrix<<std::endl;
			
		//} 
	  //std::cout<<" NodalMeanStress "<<rElementVariables.NodalMeanStress<<std::endl; 
	  //std::cout<<" ECConstitutiveMatrix 3 "<<ECConstitutiveMatrix<<std::endl; 

      // c . sigma time m
      //if(rElementVariables.voigtsize==3){
	  //double alfa = 1/2.0;
      //for (unsigned int i = 0; i < 2; i++) {
         //for (unsigned int j = 0; j < 3; j++) {
           //double aux_num = 1.0;
           //if ( j < 2){
			   //aux_num = 0.5;
		   //}
           //ECConstitutiveMatrix(i,j) -= (2 * alfa) * aux_num * (rVariables.StressVector(j)-rElementVariables.ElementalMeanStress);
         //}
      //}
      
      //for (unsigned int i = 0; i < 3; i++) {
         //for (unsigned int j = 0; j < 2; j++) {
           //double aux_num = 1.0;
           //if ( i < 2){
			   //aux_num = 0.5;
		   //}
           //ECConstitutiveMatrix(i,j) -= (2 * alfa) * aux_num * (rVariables.StressVector(i)-rElementVariables.ElementalMeanStress);
         //}
      //}
	  //}
      //else if(rElementVariables.voigtsize==6){
	  //double alfa = 1/3.0;
      //for (unsigned int i = 0; i < 3; i++) {
         //for (unsigned int j = 0; j < 6; j++) {
           //double aux_num = 1.0;
           //if ( j < 3){
			   //double aux_num = 0.5;
		   //}
           //ECConstitutiveMatrix(i,j) -= (2 * alfa) * aux_num * (rVariables.StressVector(j)-rElementVariables.ElementalMeanStress);
         //}
      //}
      
      //for (unsigned int i = 0; i < 6; i++) {
         //for (unsigned int j = 0; j < 3; j++) {
           //double aux_num = 1.0;
           //if ( i < 3){
			   //aux_num = 0.5;
		   //}
           //ECConstitutiveMatrix(i,j) -= (2 * alfa) * aux_num * (rVariables.StressVector(i)-rElementVariables.ElementalMeanStress);
         //}
      //}
	  //}
	  
	  //I define the deviatoric stress
	  
	  Vector DevAuxStress = ZeroVector(6);
	  DevAuxStress = rVariables.StressVector; 
      for (unsigned int i = 0; i < 3; i++){
         DevAuxStress(i) += ( - rElementVariables.ElementalMeanStress);
	  }
	  
	  Vector Identity = ZeroVector(6);
      
      for (unsigned int i = 0; i < 3 ; ++i) {
		  Identity(i) = 1.0; 
	  }
	  
	  
	  for (unsigned int i = 0; i < 6; i++) {
         for (unsigned int j = 0; j < 6; j++) {
           //double voigtnumber = 1.0;
           double aux_num = 0.5;
           if ( j < 3){
              //voigtnumber = 1.00;
              
              }
           ECConstitutiveMatrix(i,j) -= 2.0 * alfa * aux_num * DevAuxStress(i) * Identity(j);//(2.0/3.0) *voigtnumber* rVariables.StressVector(j);
           ECConstitutiveMatrix(i,j) -= 2.0 * alfa * aux_num * DevAuxStress(j) * Identity(i);
         }
      }
      //if(this->Id() == 19289)
        //{
			//std::cout<<" DevAuxStress 3"<<DevAuxStress<<std::endl;
			//std::cout<<" ECConstitutiveMatrix 3"<<ECConstitutiveMatrix<<std::endl;
			
		//} 
	  //for (unsigned int i = 0; i < 6; i++) {
         //for (unsigned int j = 0; j < 3; j++) {
           //double aux_num = 1.0;
           //if ( i < 3){
			   //aux_num = 0.5;
		   //}
           //ECConstitutiveMatrix(i,j) -= 2.0 * alfa * aux_num * DevAuxStress(i);
         //}
      //}
	  
	 //if(this->Id() == 9794)
        //{
			//std::cout<<" ECConstitutiveMatrix 4"<<ECConstitutiveMatrix<<std::endl;
			
		//} 
      
     if ( rElementVariables.voigtsize == 6)
      {
         ConstitutiveMatrix = ECConstitutiveMatrix;
      }
      else
      {
         int indexi , indexj; 
         for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = 0; j < 3 ; j++)
            {
               indexi = i;
               indexj = j;
               if ( i == 2)
                  indexi += 1;
               if ( j == 2)
                  indexj += 1;

               ConstitutiveMatrix(i,j) = ECConstitutiveMatrix(indexi, indexj);
            }
         }
      }
	  //std::cout<<" ConstitutiveMatrix in Kuum2 "<<ConstitutiveMatrix<<std::endl; 
      Matrix Kuu = prod( trans( rVariables.B ),  rIntegrationWeight * Matrix( prod( ConstitutiveMatrix, rVariables.B ) ) ); 
	  //std::cout<< " rVariables.B "<< rVariables.B<<std::endl;
	  //std::cout<< " ConstitutiveMatrix "<< ConstitutiveMatrix<<std::endl;
      // MatrixType Kh=rLeftHandSideMatrix;

      unsigned int indexi = 0;
      unsigned int indexj  = 0;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int idim = 0; idim < dimension ; idim ++)
         {
            indexj=0;
            for ( unsigned int j = 0; j < number_of_nodes; j++ )
            {
               for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
               {
                  rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuu(indexi,indexj);
                  indexj++;
               }
            }
            indexi++;
         }
      }
// */
	//if(this->Id() == 8188)
		//{
		//std::cout<<" Kmat "<<rLeftHandSideMatrix-Kh<<std::endl;
		//}
      // std::cout<<std::endl;
      //std::cout<<" Kmat "<<rLeftHandSideMatrix-Kh<<std::endl;

      KRATOS_CATCH( "" )

   }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
            GeneralVariables& rVariables,
            ThisElementGeneralVariables& rElementVariables, 
            double& rIntegrationWeight)

    {
        KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        int size = number_of_nodes * dimension;
        
        MatrixType Kh=rLeftHandSideMatrix;
        
        Vector StressVector = rElementVariables.StressVector;
        //if(this->Id() == 8188)
        //{
        //std::cout<<" StressVector of element use in GeoMatrix "<< rElementVariables.StressVector<<std::endl;
        //std::cout<<" StressVector use in GeoMatrix "<< rVariables.StressVector<<std::endl;
		//}
        Matrix StressTensor = MathUtils<double>::StressVectorToTensor(StressVector );
        
        Matrix ReducedKg = prod( rVariables.DN_DX, rIntegrationWeight * Matrix( prod( StressTensor, trans( rVariables.DN_DX ) ) ) ); //to be optimized
        
        Matrix Kuug = zero_matrix<double> (size);
        MathUtils<double>::ExpandAndAddReducedMatrix( Kuug, ReducedKg, dimension );
        
        unsigned int indexi = 0;
		unsigned int indexj = 0;
		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			for ( unsigned int idim = 0; idim < dimension ; idim ++)
			{
				indexj=0;
				for ( unsigned int j = 0; j < number_of_nodes; j++ )
				{
					for ( unsigned int jdim = 0; jdim < dimension ; jdim ++)
					{
						rLeftHandSideMatrix(indexi+i,indexj+j)+=Kuug(indexi,indexj);
						indexj++;
					}
				}
				indexi++;
			}
		}
        //if(this->Id() == 19289)
        //{
			////std::cout<<" ReducedKg "<<ReducedKg<<std::endl;
			////std::cout<<" Kuug "<<Kuug<<std::endl;
			//std::cout<<" Kgeo "<<rLeftHandSideMatrix-Kh<<std::endl;
		//}
		//std::cout<<" Kgeo "<<rLeftHandSideMatrix-Kh<<std::endl;
        KRATOS_CATCH( "" )
    }
    
//************************************************************************************
//************************************************************************************

	void UpdatedLagrangianUP2::CalculateAndAddKup (MatrixType& rLeftHandSideMatrix,
			GeneralVariables& rVariables,
			ThisElementGeneralVariables& rElementVariables, 
			double& rIntegrationWeight)
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		
		//unsigned int voigtsize  = 3;
		
            
        //if( dimension == 3 )
        //{
            //voigtsize  = 6;
        //}
		
		//here I call a public method of CL class which evaluate the additional term due to the 
		//derivative of the deviatoric term respect to the pressure
		
		//Matrix MixedUPTerm = mConstitutiveLawVector->GetValue(MIXED_UP_TERM, MixedUPTerm );
		//std::cout<<"MixedUPTerm "<<MixedUPTerm<<std::endl;
		//Vector MixedUPTermVector = MathUtils<double>::StressTensorToVector( MixedUPTerm, voigtsize );
		 
		MatrixType Kh=rLeftHandSideMatrix;
		
		//Matrix ShapeFunction = ZeroMatrix(1,3);
		//Matrix IdentityMatrix = ZeroMatrix(3,1);
		//for(unsigned int i=0; i<3 ; i++)
		//{
			//ShapeFunction(0,i) = rVariables.N[i];
			//IdentityMatrix(i,0) = 1;
		//}
		
		//Matrix TestMatrix = rIntegrationWeight * prod(trans(rVariables.B ), Matrix(prod(IdentityMatrix,ShapeFunction )));
		
		
		
		//for ( unsigned int i = 0; i < number_of_nodes; i++ )
		//{
			////double A1 = prod(trans( rVariables.DN_DX( i , k ) ), MixedUPTermVector);
			//unsigned int indexp  = dimension;
			//unsigned int indexup = dimension * i + i;
			//unsigned int indexi = i * 2;
			//for ( unsigned int j = 0; j < number_of_nodes; j++ )
			//{
				//for ( unsigned int k = 0; k < dimension; k++ )
				//{
				//rLeftHandSideMatrix(indexup + k,indexp) +=  TestMatrix(indexi + k,j);
				//}
				//indexp += (dimension + 1);
			//}
			
		//}
		//contributions to stiffness matrix calculated on the reference configuration
		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			//double A1 = prod(trans( rVariables.DN_DX( i , k ) ), MixedUPTermVector);
			unsigned int indexp  = dimension;
			unsigned int indexup = dimension * i + i;
			for ( unsigned int j = 0; j < number_of_nodes; j++ )
			{
				
				for ( unsigned int k = 0; k < dimension; k++ )
				{
					rLeftHandSideMatrix(indexup+k,indexp) +=  rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight;// * rVariables.detF;
					
					//rLeftHandSideMatrix(indexup+k,indexp) +=  MixedUPTerm(i,k) * rVariables.DN_DX ( i , k ) *  rVariables.N[j] * rIntegrationWeight;
					
				}
				indexp += (dimension + 1);
			}
		}
		//if(this->Id() == 19289)
		//{
		//std::cout<<" Kup "<<rLeftHandSideMatrix-Kh<<std::endl;
		//}
		 // std::cout<<std::endl;
		 //std::cout<<" Kup "<<rLeftHandSideMatrix-Kh<<std::endl;

		KRATOS_CATCH( "" )
	}    

//************************************************************************************
//************************************************************************************

	//void UpdatedLagrangianUP2::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
			//GeneralVariables& rVariables,
			//ThisElementGeneralVariables& rElementVariables, 
			//double& rIntegrationWeight)

	//{
		//KRATOS_TRY

		////repasar

		//const unsigned int number_of_nodes = GetGeometry().size();
		//const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		//// MatrixType Kh=rLeftHandSideMatrix;

		////contributions to stiffness matrix calculated on the reference configuration
		//unsigned int indexp = dimension;

		//double DeltaCoefficient = 0;
		//DeltaCoefficient = this->CalculatePUDeltaCoefficient( DeltaCoefficient, rVariables ); 


		//for ( unsigned int i = 0; i < number_of_nodes; i++ )
		//{
			//for ( unsigned int j = 0; j < number_of_nodes; j++ )
			//{
				//int indexup= dimension*j + j;
				//for ( unsigned int k = 0; k < dimension; k++ )
				//{
			  //rLeftHandSideMatrix(indexp,indexup+k) +=  DeltaCoefficient  * rVariables.N[i] * rVariables.DN_DX ( j , k ) * rIntegrationWeight * rVariables.detF;

					////std::cout<<" value ("<<indexp<<","<<indexup+k<<") "<<(2*detF) * rN[i] * rDN_DX ( j , k ) * rIntegrationWeight<<std::endl;
				//}
			//}
			//indexp += (dimension + 1);
		//}


		//// std::cout<<std::endl;
		//// std::cout<<" Kpu "<<rLeftHandSideMatrix-Kh<<std::endl;


		//KRATOS_CATCH( "" )
	//}
	void UpdatedLagrangianUP2::CalculateAndAddKpu (MatrixType& rLeftHandSideMatrix,
         GeneralVariables& rVariables,
         ThisElementGeneralVariables& rElementVariables, 
         double& rIntegrationWeight)

   {
      KRATOS_TRY

     

      const unsigned int number_of_nodes = GetGeometry().size();
      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      Matrix Kpu;
      MatrixType Kh=rLeftHandSideMatrix;
		
     
      Matrix ECConstitutiveMatrix = rVariables.ConstitutiveMatrix;
      
	  Matrix Identity = ZeroMatrix(1,6);
      
      for (unsigned int i = 0; i < 3 ; ++i) {
		  Identity(0,i) = 1.0; 
	  }
	  Matrix AuxConstitutiveVector = prod( Identity, ECConstitutiveMatrix);
	  AuxConstitutiveVector *= ( 1.0/ 3.0 );
	  
      
      for (unsigned int i = 0; i < 6; i++) 
      {
		  double alfa = 1/3.0; 
		  AuxConstitutiveVector(0, i) += 2.0 * alfa * rVariables.StressVector(i);	
	  
	  
	  }
      
	    
      Matrix AuxPressureVector = rElementVariables.ElementalMeanStress * Identity;
      //std::cout<<"AuxPressureVector "<<AuxPressureVector<<std::endl;
      // - ElemenPressure * 1
      for (unsigned int i = 0; i < 6; i++)
      {
          AuxConstitutiveVector(0,i) -= AuxPressureVector(0,i);
	  }
	  
	 
       
     
     Matrix ConstitutiveVector = ZeroMatrix( 1, rElementVariables.voigtsize);
      if ( rElementVariables.voigtsize == 6)
      {
         ConstitutiveVector = AuxConstitutiveVector;
      }
      else
      {
         int indexi; 
         for (unsigned int i = 0; i < 3; i++) {
            indexi = i;
            if ( i == 2)
               indexi += 1;

            ConstitutiveVector(0,i) = AuxConstitutiveVector(0,indexi);
         }
      }
		
	  
      Kpu = prod( ConstitutiveVector, rVariables.B);
      
     
      

      Matrix Kpu2 = ZeroMatrix( number_of_nodes, number_of_nodes*dimension);
      for (unsigned int i = 0; i < number_of_nodes;  i++) {
         for (unsigned int j = 0; j < number_of_nodes*dimension; j++) {
            Kpu2 ( i, j ) = rVariables.N[i] * Kpu(0,j);
         }
      }
      
	  
      Kpu2 *= rIntegrationWeight;
      
      unsigned int indexp = dimension;
      for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
         for ( unsigned int j = 0; j < number_of_nodes; j++ )
         {
            int indexup= dimension*j + j;
            for ( unsigned int k = 0; k < dimension; k++ )
            {
               rLeftHandSideMatrix(indexp,indexup+k) +=  Kpu2( i, j*2 + k)*mElementScalingNumber;
            }
         }
         indexp += (dimension + 1);
      }
		
		
      //std::cout<<" Kpu "<<rLeftHandSideMatrix-Kh<<std::endl;
     //std::cout<<"38888888888888888888"<<std::endl;	
     //std::cout<<" rLeftHandSideMatrix "<< this->Id()<<" "<<rLeftHandSideMatrix<<std::endl;
// */
      KRATOS_CATCH( "" )
   }
//************************************************************************************
//************************************************************************************

	void UpdatedLagrangianUP2::CalculateAndAddKpp (MatrixType& rLeftHandSideMatrix,
			GeneralVariables& rVariables,
			ThisElementGeneralVariables& rElementVariables, 
			double& rIntegrationWeight)
	{
		KRATOS_TRY


		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		//double BulkModulus= GetProperties()[YOUNG_MODULUS]/(3*(1-2*GetProperties()[POISSON_RATIO]));

		MatrixType Kh=rLeftHandSideMatrix;

		//contributions to stiffness matrix calculated on the reference configuration
		unsigned int indexpi = dimension;
		//double consistent = 1.0;

		for ( unsigned int i = 0; i < number_of_nodes; i++ )
		{
			unsigned int indexpj = dimension;
			for ( unsigned int j = 0; j < number_of_nodes; j++ )
		  {
			 
			// consistent=1;
			// if(indexpi==indexpj)
			//   consistent=2;

			// if( dimension == 2 ){ //consistent 2D
			  
			//   rLeftHandSideMatrix(indexpi,indexpj)  -= consistent * ((1.0)/(BulkModulus)) * (1.0/12.0) * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D

			// }
			// else{

			//   rLeftHandSideMatrix(indexpi,indexpj)  -= consistent * ((1.0)/(BulkModulus)) * (1.0/20.0) * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //3D
			  
			// }

			//rLeftHandSideMatrix(indexpi,indexpj)  -= ((1.0)/(BulkModulus)) * rVariables.N[i] * rVariables.N[j] * rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D-3D
			
			//********************************************************
			//mixed formulation considering a variation of p: 2nd term
			rLeftHandSideMatrix(indexpi,indexpj)  -= rVariables.N[i] * rVariables.N[j] * rIntegrationWeight / (rVariables.detF0/rVariables.detF)*mElementScalingNumber;
			
				indexpj += (dimension + 1);
		  }

			indexpi += (dimension + 1);
		}
		
		
		//for ( unsigned int i = 0; i < number_of_nodes; i++ )
		//{
			//unsigned int indexpj = dimension;
			//for ( unsigned int j = 0; j < number_of_nodes; j++ )
		  //{
			  //rLeftHandSideMatrix(indexpj,indexpi)  -= ((1.0)/(BulkModulus)) * rVariables.N[i] * rIntegrationWeight / (rVariables.detF0/rVariables.detF);
		
			  //indexpj += (dimension + 1);
		
		  //}
		  
		    //indexpi += (dimension + 1);
		//}
		
		//if(this->Id() == 19289)
		//{
		//std::cout<<" Kpp "<<rLeftHandSideMatrix-Kh<<std::endl;
		//}
	
		

		// std::cout<<std::endl;
		//std::cout<<" Kpp "<<rLeftHandSideMatrix-Kh<<std::endl;

		KRATOS_CATCH( "" )
	}
//************************************************************************************
//************************************************************************************
// I changed the constant matrix in the stabilized term:
// as in MPM the position of the integration points does not coincide with the
// position of the Gauss points the first matrix is substitute with the product of the
// shape function values of each integration point
//************************************************************************************
//************************************************************************************

	void UpdatedLagrangianUP2::CalculateAndAddKppStab (MatrixType& rLeftHandSideMatrix,
			GeneralVariables & rVariables,
			ThisElementGeneralVariables& rElementVariables, 
			double& rIntegrationWeight)
	{
		KRATOS_TRY

		//repasar

		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		MatrixType Kh=rLeftHandSideMatrix;

		//contributions to stiffness matrix calculated on the reference configuration
		unsigned int indexpi = dimension;
		
		double AlphaStabilization  = 8.0; 
		double StabilizationFactor = GetProperties()[STABILIZATION_FACTOR];
		AlphaStabilization *= StabilizationFactor;

		const double& YoungModulus          = GetProperties()[YOUNG_MODULUS];
		const double& PoissonCoefficient    = GetProperties()[POISSON_RATIO];

		double LameMu =  YoungModulus/(2*(1+PoissonCoefficient));
		double BulkModulus = YoungModulus/ 3.0 / ( 1.0 - 2.0*PoissonCoefficient) ;
		
		
		//Experimental
		// if(LameMu < rVariables.ConstitutiveMatrix(2,2))
		//   LameMu = rVariables.ConstitutiveMatrix(2,2);

		double consistent = 1.0;
		AlphaStabilization=(AlphaStabilization/(36.0*LameMu));
        AlphaStabilization *= BulkModulus;
        
        //std::cout<<" AlphaStabilization "<<AlphaStabilization<<std::endl;
        //std::cout<<" BulkModulus "<<BulkModulus<<std::endl;
        //std::cout<<" LameMu "<<LameMu<<std::endl;
        
		//double FactorValue = 8.0; //JMR deffault value	
		//if( dimension == 3 )
		  //FactorValue = 10.0; //JMC deffault value
		
		//consistent = AlphaStabilization*FactorValue/LameMu;
		
		for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int indexpj = dimension;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
	  {

	    if( dimension == 2 ){ //consistent 2D
	      
	      consistent=(-1)*AlphaStabilization;//*FactorValue/(36.0*LameMu);
	      if(indexpi==indexpj)
                consistent=2*AlphaStabilization;//*FactorValue/(36.0*LameMu);

	      //rLeftHandSideMatrix(indexpi,indexpj) -= consistent *(rVariables.N[i]*rVariables.N[j] -1/9)*rIntegrationWeight / (rVariables.detF0/rVariables.detF); //2D
	      rLeftHandSideMatrix(indexpi,indexpj) -= consistent *rIntegrationWeight / (rVariables.detF0/rVariables.detF)*mElementScalingNumber; //2D
	    
	    }
	    else{

	      consistent=(-1)*AlphaStabilization;//*FactorValue/(80.0*LameMu);
	      if(indexpi==indexpj)
                consistent=3*AlphaStabilization;//*FactorValue/(80.0*LameMu);
	      
	      rLeftHandSideMatrix(indexpi,indexpj) -= consistent * rIntegrationWeight / (rVariables.detF0/rVariables.detF);//*mElementScalingNumber; //3D

	    }


            indexpj += (dimension + 1);
	  }

        indexpi += (dimension + 1);
      }
      
//******************************************************************************************************************************************************      
      
     
      //if(this->Id() == 19289)
		//{
		//std::cout<<" KppStab "<<rLeftHandSideMatrix-Kh<<std::endl;
		//}
      
      

		// std::cout<<std::endl;
		//std::cout<<" KppStab "<<rLeftHandSideMatrix-Kh<<std::endl;

		KRATOS_CATCH( "" )
	}


//************************************CALCULATE VOLUME CHANGE*************************
//************************************************************************************

    double& UpdatedLagrangianUP2::CalculateVolumeChange( double& rVolumeChange, GeneralVariables& rVariables )
    {
        KRATOS_TRY
          
        rVolumeChange = 1.0 / (rVariables.detF * rVariables.detF0); 

        return rVolumeChange;

        KRATOS_CATCH( "" )
    }


//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int element_size          = number_of_nodes * dimension + number_of_nodes;

    if ( rResult.size() != element_size )
        rResult.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension + i;
        rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();

        if( dimension == 3)
        {
            rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
            rResult[index + 3] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
        }
        else
        {
            rResult[index + 2] = GetGeometry()[i].GetDof( PRESSURE ).EquationId();
        }

    }
        
    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo )
    {
        rElementalDofList.resize( 0 );

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );

        rElementalDofList.push_back( GetGeometry()[i].pGetDof( PRESSURE ));
    }
        //std::cout<< "ElementalDofList.size() "<<rElementalDofList.size()<<std::endl;
        //std::cout<<" in GetDofList of derived class"<<std::endl;
    }



//************************************************************************************
//*******************DAMPING MATRIX***************************************************

    //void UpdatedLagrangian::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
    //{
        //KRATOS_TRY

        ////0.-Initialize the DampingMatrix:
        //unsigned int number_of_nodes = GetGeometry().size();
        //unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        ////resizing as needed the LHS
        //unsigned int MatSize = number_of_nodes * dimension;

        //if ( rDampingMatrix.size1() != MatSize )
            //rDampingMatrix.resize( MatSize, MatSize, false );

        //noalias( rDampingMatrix ) = ZeroMatrix( MatSize, MatSize );


        ////1.-Calculate StiffnessMatrix:

        //MatrixType StiffnessMatrix  = Matrix();

        //this->CalculateLeftHandSide( StiffnessMatrix, rCurrentProcessInfo );

        ////2.-Calculate MassMatrix:

        //MatrixType MassMatrix  = Matrix();

        //this->CalculateMassMatrix ( MassMatrix, rCurrentProcessInfo );
        
        
        ////3.-Get Damping Coeffitients (RAYLEIGH_ALPHA, RAYLEIGH_BETA)
        //double alpha = 0;
        //if( GetProperties().Has(RAYLEIGH_ALPHA) ){
          //alpha = GetProperties()[RAYLEIGH_ALPHA];
        //}
        //else if( rCurrentProcessInfo.Has(RAYLEIGH_ALPHA) ){
          //alpha = rCurrentProcessInfo[RAYLEIGH_ALPHA];
        //}

        //double beta  = 0;
        //if( GetProperties().Has(RAYLEIGH_BETA) ){
          //beta = GetProperties()[RAYLEIGH_BETA];
        //}
        //else if( rCurrentProcessInfo.Has(RAYLEIGH_BETA) ){
          //beta = rCurrentProcessInfo[RAYLEIGH_BETA];
        //}

        ////4.-Compose the Damping Matrix:
       
        ////Rayleigh Damping Matrix: alpha*M + beta*K
        //rDampingMatrix  = alpha * MassMatrix;
        //rDampingMatrix += beta  * StiffnessMatrix;
        ////std::cout<<" rDampingMatrix "<<rDampingMatrix<<std::endl;

        //KRATOS_CATCH( "" )
    //}
//************************************************************************************
//****************MASS MATRIX*********************************************************

     void UpdatedLagrangianUP2::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //I need to call the values of the shape function for the single element
        GeneralVariables Variables;
        this->InitializeGeneralVariables(Variables,rCurrentProcessInfo);

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        unsigned int MatSize = number_of_nodes * dimension + number_of_nodes;

        if ( rMassMatrix.size1() != MatSize )
            rMassMatrix.resize( MatSize, MatSize, false );

        rMassMatrix = ZeroMatrix( MatSize, MatSize );

        double TotalMass = 0;

        //TOTAL MASS OF ONE MP ELEMENT
        
        TotalMass = this->GetValue(MP_MASS);
        
        //LUMPED MATRIX
        
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
        double temp = Variables.N[i] * TotalMass;
        unsigned int indexup = i * dimension + i;
        for ( unsigned int j = 0; j < dimension; j++ )
        {
            //unsigned int index = i * dimension + j;
            rMassMatrix( indexup+j , indexup+j ) = temp;
        }
        }
        //for ( unsigned int i = 0; i < number_of_nodes; i++ )
      	//{
      	  //unsigned int indexupi = dimension * i + i;

      	  //for ( unsigned int j = 0; j < number_of_nodes; j++ )
      	    //{
      	      //unsigned int indexupj = dimension * j + j;

      	      //for ( unsigned int k = 0; k < dimension; k++ )
      		//{
      		  //rMassMatrix( indexupi+k , indexupj+k ) += Variables.N[i] * Variables.N[j] * TotalMass;
      		//}
      	    //}
      	//}
        
        

        //std::cout<<"rMassMatrix "<<rMassMatrix<<std::endl;

        KRATOS_CATCH( "" )
    }

//************************************************************************************
//************************************************************************************

 

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::GetValuesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
		unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
        {
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
            values[index + 3] = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
        }
        else
        {
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( PRESSURE, Step );
        }

    }
    }


//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::GetFirstDerivativesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );
        if ( dimension == 3 )
        {
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
            values[index + 3] = 0;
        }
        else
        {
            values[index + 2] = 0;
        }
    }
    }

//************************************************************************************
//************************************************************************************

    void UpdatedLagrangianUP2::GetSecondDerivativesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       element_size    = number_of_nodes * dimension + number_of_nodes;

    if ( values.size() != element_size ) values.resize( element_size, false );


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension + i;
        values[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
        {
            values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
            values[index + 3] = 0;
        }
        else
        {
            values[index + 2] = 0;
        }
    }
    }
//************************************************************************************
//************************************************************************************
    void UpdatedLagrangianUP2::GetHistoricalVariables( GeneralVariables& rVariables )
    {
        //Deformation Gradient F ( set to identity )
        unsigned int size =  rVariables.F.size1();
        rVariables.detF  = 1;
        rVariables.F     = IdentityMatrix(size);
        
        rVariables.detF0 = mDeterminantF0;
        rVariables.F0    = mDeformationGradientF0;

    }
    
////************************************************************************************

    void UpdatedLagrangianUP2::FinalizeStepVariables( GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo)
    {
    UpdatedLagrangian::FinalizeStepVariables( rVariables, rCurrentProcessInfo);
    ThisElementGeneralVariables ElementVariables;
    CalculateThisElementGeneralVariables( ElementVariables, rVariables);
    
    this->SetValue(MP_CAUCHY_STRESS_VECTOR, ElementVariables.StressVector);
    
    
    }    
//*************************************************************************************
   void UpdatedLagrangianUP2::CalculateThisElementGeneralVariables( ThisElementGeneralVariables& rElementGeneralVariables, const GeneralVariables & rVariables)
   {
		
      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
      unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      rElementGeneralVariables.voigtsize = 3;
      if ( dimension == 3)
         rElementGeneralVariables.voigtsize = 6;

      rElementGeneralVariables.NodalMeanStress = 0;
      for (unsigned int i = 0; i < number_of_nodes; i++)
         rElementGeneralVariables.NodalMeanStress += GetGeometry()[i].GetSolutionStepValue( PRESSURE ) * rVariables.N[i];

      rElementGeneralVariables.ElementalMeanStress = 0;
      for (unsigned int i = 0; i < 3 ; i++)
         rElementGeneralVariables.ElementalMeanStress += rVariables.StressVector(i);
      rElementGeneralVariables.ElementalMeanStress /= 3.0;

      Vector AuxStress = ZeroVector(6);
      AuxStress = rVariables.StressVector; 
      //if(this->Id() == 443)
	  //{
		  
		  //std::cout<<"AuxStress 1"<<AuxStress<<std::endl;
		  
	  //}
      for (unsigned int i = 0; i < 3; i++)
         AuxStress(i) += ( rElementGeneralVariables.NodalMeanStress - rElementGeneralVariables.ElementalMeanStress);

      
      //if(this->Id() == 443)
	  //{
		  
		  //std::cout<<"AuxStress 2"<<AuxStress<<std::endl;
		  //std::cout<<"p"<<rElementGeneralVariables.NodalMeanStress <<std::endl;
		  //std::cout<<"p(u)"<<rElementGeneralVariables.ElementalMeanStress <<std::endl;
		  //std::cout<<"p-p(u)"<<rElementGeneralVariables.NodalMeanStress - rElementGeneralVariables.ElementalMeanStress<<std::endl;
		  
	  //}
      rElementGeneralVariables.StressVector = ZeroVector(rElementGeneralVariables.voigtsize);

      if ( rElementGeneralVariables.voigtsize == 6) {
         rElementGeneralVariables.StressVector = AuxStress; 
      }
      else {
         rElementGeneralVariables.StressVector(0) = AuxStress(0);
         rElementGeneralVariables.StressVector(1) = AuxStress(1);
         rElementGeneralVariables.StressVector(2) = AuxStress(3);
      }
	  //if(this->Id() == 443)
	  //{
		  
		  //std::cout<<"AuxStress 3"<<AuxStress<<std::endl;
		  //std::cout<<"rElementGeneralVariables.StressVector "<<rElementGeneralVariables.StressVector<<std::endl;
	  //}
      rElementGeneralVariables.DeviatoricTensor = ZeroMatrix(6,6);
      for (unsigned int i = 0; i < 6; i++)
         rElementGeneralVariables.DeviatoricTensor(i,i) = 1;

      for (unsigned int i = 0; i < 3; i++) {
         for (unsigned int j = 0; j < 3; j++) {
         rElementGeneralVariables.DeviatoricTensor(i,j) -= 1.0/3.0;
         }
      }



   }
//************************************************************************************
//************************************************************************************


//************************************************************************************
//************************************************************************************
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
    int  UpdatedLagrangianUP2::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        int correct = 0;

		correct = UpdatedLagrangian::Check(rCurrentProcessInfo);


		//verify compatibility with the constitutive law
		//ConstitutiveLaw::Features LawFeatures;
		//this->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetLawFeatures(LawFeatures);

		//if(LawFeatures.mOptions.IsNot(ConstitutiveLaw::U_P_LAW))
			//KRATOS_THROW_ERROR( std::logic_error, "constitutive law is not compatible with the U-P element type ", " Large Displacements U_P" )

		//verify that the variables are correctly initialized

		if ( PRESSURE.Key() == 0 )
			KRATOS_THROW_ERROR( std::invalid_argument, "PRESSURE has Key zero! (check if the application is correctly registered", "" )

		return correct;
       
       

        KRATOS_CATCH( "" );
    }

    


    void UpdatedLagrangianUP2::save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
        //int IntMethod = int(mThisIntegrationMethod);
        //rSerializer.save("IntegrationMethod",IntMethod);
        rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
        rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
        rSerializer.save("DeterminantF0",mDeterminantF0);
        

    }

    void UpdatedLagrangianUP2::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
        rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
        rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
        rSerializer.load("DeterminantF0",mDeterminantF0);
        //int IntMethod;
        //rSerializer.load("IntegrationMethod",IntMethod);
        //mThisIntegrationMethod = IntegrationMethod(IntMethod);
        

    }





} // Namespace Kratos

