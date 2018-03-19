#include "includes/define.h"
#include <string>
#include "includes/constitutive_law.h"
#include "custom_constitutive/zarate_law.hpp"
#include "femdem3d_element.hpp"
#include "romfemdem3d_element.hpp"
#include "includes/element.h"
#include "includes/node.h"
#include "fem_to_dem_application_variables.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"
//#include "solid_mechanics_application_variables.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{
	//***********************DEFAULT CONSTRUCTOR******************************************
	//************************************************************************************

	RomFemDem3DElement::RomFemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: FemDem3DElement(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}
	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	RomFemDem3DElement::RomFemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		:FemDem3DElement(NewId, pGeometry, pProperties)
	{
		//BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
		mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
	}

	//******************************COPY CONSTRUCTOR**************************************
	//************************************************************************************

	RomFemDem3DElement::RomFemDem3DElement(RomFemDem3DElement const& rOther)
		:FemDem3DElement(rOther)
	{
		//ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
		//IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
	}

	//*******************************ASSIGMENT OPERATOR***********************************
	//************************************************************************************

	RomFemDem3DElement&  RomFemDem3DElement::operator=(RomFemDem3DElement const& rOther)
	{
		//ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

		FemDem3DElement::operator=(rOther);
		return *this;
	}

	//*********************************OPERATIONS*****************************************
	//************************************************************************************

	Element::Pointer RomFemDem3DElement::Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
	{
		//NEEDED TO CREATE AN ELEMENT   
		return Element::Pointer(new RomFemDem3DElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
	}


	//************************************CLONE*******************************************
	//************************************************************************************

	Element::Pointer RomFemDem3DElement::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
	{

		//YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
		//ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

		RomFemDem3DElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

		return Element::Pointer(new RomFemDem3DElement(NewElement));
	}

	//*******************************DESTRUCTOR*******************************************
	//************************************************************************************

	RomFemDem3DElement::~RomFemDem3DElement()
	{
	}


	void RomFemDem3DElement::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
	{
		//*****************************
		KRATOS_TRY

		bool is_active = true;
		if (this->IsDefined(ACTIVE))
		{
			is_active = this->Is(ACTIVE);
		}

		// Inactive elements can have negative determinant of the Jacobian
		if (is_active == true)
		{
			//1.-Initialize sizes for the system components:
			const unsigned int number_of_nodes = GetGeometry().size();
			const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
			unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

			Vector StrainVector(voigt_size);
			noalias(StrainVector) = ZeroVector(voigt_size);
			Vector StressVector(voigt_size);
			noalias(StressVector) = ZeroVector(voigt_size);
			Matrix ConstitutiveMatrix(voigt_size, voigt_size);
			noalias(ConstitutiveMatrix) = ZeroMatrix(voigt_size, voigt_size);
			Matrix B(voigt_size, dimension*number_of_nodes);
			noalias(B) = ZeroMatrix(voigt_size, dimension*number_of_nodes);
			Matrix DN_DX(number_of_nodes, dimension);
			noalias(DN_DX) = ZeroMatrix(number_of_nodes, dimension);


			//deffault values for the infinitessimal theory
			double detF = 1;
			Matrix F(dimension, dimension);
			noalias(F) = identity_matrix<double>(dimension);

			//3.-Calculate elemental system:

			//reading integration points
			const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

			//get the shape functions [N] (for the order of the default integration method)
			const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

			//get the shape functions parent coodinates derivative [dN/d�] (for the order of the default integration method)
			const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);

			//calculate delta position (here coincides with the current displacement)
			Matrix DeltaPosition(number_of_nodes, dimension);
			noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
			DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

			//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
			GeometryType::JacobiansType J;
			J.resize(1, false);
			J[0].resize(dimension, dimension, false);
			noalias(J[0]) = ZeroMatrix(dimension, dimension);
			J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);

			// Loop Over Integration Points
			for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
			{
				Matrix InvJ(dimension, dimension);
				noalias(InvJ) = ZeroMatrix(dimension, dimension);
				double detJ = 0;
				MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

				if (detJ < 0)
				{
					this->Set(ACTIVE, false); // element alone inside a crack
					detJ = fabs(detJ);
				}

				if (detJ < 0) KRATOS_THROW_ERROR(std::invalid_argument, " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ)

				//compute cartesian derivatives for this integration point  [dN/dx_n]
				noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

				//set shape functions for this integration point
				Vector N = row(Ncontainer, PointNumber);

				//b.-compute infinitessimal strainof the composite
				this->CalculateInfinitesimalStrain(StrainVector, DN_DX);
				this->SetValue(STRAIN_VECTOR, StrainVector);

                // Compute predictive stresses for the concrete and steel
                this->CalculatePredictiveStresses(StrainVector);

				this->CalculateDeformationMatrix(B, DN_DX);
				this->SetBMatrix(B);

			}
		}
		
		KRATOS_CATCH("")
	}

    void RomFemDem3DElement::CalculatePredictiveStresses(const Vector& StrainVector)
    {
        double Ec,Es,nuc,nus;
        Ec  = this->GetProperties()[YOUNG_MODULUS];
        Es  = this->GetProperties()[YOUNG_MODULUS_STEEL];
        nuc = this->GetProperties()[POISSON_RATIO];
        nus = this->GetProperties()[POISSON_RATIO_STEEL];

        const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int voigt_size = dimension * (dimension + 1) * 0.5;
        Matrix ConstitutiveMatrixConcrete = ZeroMatrix(voigt_size, voigt_size);
        Matrix ConstitutiveMatrixSteel    = ZeroMatrix(voigt_size, voigt_size);

        // Elastic C
        this->CalculateConstitutiveMatrix(ConstitutiveMatrixConcrete, Ec, nuc);
        this->CalculateConstitutiveMatrix(ConstitutiveMatrixSteel, Es, nus);

        Vector StressVectorConcrete = prod(ConstitutiveMatrixConcrete, StrainVector);

        Vector StressVectorSteel;
        if (this->GetProperties()[STEEL_VOLUMETRIC_PART] > 0.0)
        {
            StressVectorSteel = prod(ConstitutiveMatrixSteel, StrainVector);
        }
        else StressVectorSteel = ZeroVector(voigt_size);
        

        // Predictive Stresses
        this->SetValue(CONCRETE_STRESS_VECTOR, StressVectorConcrete);
        this->SetValue(STEEL_STRESS_VECTOR, StressVectorSteel);

    }

	void RomFemDem3DElement::CalculateAverageStressOnEdge(Vector& rAverageVector, const std::vector<Element*> VectorOfElems)
	{
        // Only averages the stress over the concrete part!!!!
		Vector CurrentElementStress = this->GetValue(CONCRETE_STRESS_VECTOR);
		rAverageVector = CurrentElementStress;
		int counter = 0;

		for (int elem = 0; elem < VectorOfElems.size(); elem++)
		{
			// Only take into account the active elements
			bool is_active = true;
			if (VectorOfElems[elem]->IsDefined(ACTIVE))
			{
				is_active = VectorOfElems[elem]->Is(ACTIVE);
			}

			if (is_active == true)
			{
				rAverageVector += VectorOfElems[elem]->GetValue(CONCRETE_STRESS_VECTOR);
				counter++;
			}
		}
		rAverageVector /= (counter + 1);
	}


	void RomFemDem3DElement::CalculateLocalSystem (MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
		unsigned int system_size = number_of_nodes * dimension;
		if (rLeftHandSideMatrix.size1() != system_size) rLeftHandSideMatrix.resize(system_size, system_size, false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(system_size, system_size); 
																 
		if (rRightHandSideVector.size() != system_size) rRightHandSideVector.resize(system_size, false);
		noalias(rRightHandSideVector) = ZeroVector(system_size); 

		Matrix DeltaPosition(number_of_nodes, dimension);
		noalias(DeltaPosition) = ZeroMatrix(number_of_nodes, dimension);
		DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);

		GeometryType::JacobiansType J;
		J.resize(1, false);
		J[0].resize(dimension, dimension, false);
		noalias(J[0]) = ZeroMatrix(dimension, dimension);
		J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);
		
		for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
		{
			const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);
			Vector N = row(Ncontainer, PointNumber);

			double detJ = 0;
			Matrix InvJ(dimension, dimension);
			noalias(InvJ) = ZeroMatrix(dimension, dimension);
			MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

			double IntegrationWeight = integration_points[PointNumber].Weight() * detJ;
			const Matrix& B = this->GetBMatrix();
			Vector IntegratedStressVectorConcrete = ZeroVector(voigt_size);
			Vector DamagesOnEdges = ZeroVector(6);
			
			// Loop over edges of the element
			for (int edge = 0; edge < 6; edge++)
			{
				std::vector<Element*> EdgeNeighbours = this->GetEdgeNeighbourElements(edge);

				Vector AverageStressVectorConcrete, AverageStrainVectorConcrete, IntegratedStressVectorOnEdge;
				this->CalculateAverageStressOnEdge(AverageStressVectorConcrete, EdgeNeighbours);
				this->CalculateAverageStrainOnEdge(AverageStrainVectorConcrete, EdgeNeighbours);

				double DamageEdge = 0.0;
				double Lchar = this->Get_l_char(edge);

                // Integrate the stress on edge
				this->IntegrateStressDamageMechanics(IntegratedStressVectorOnEdge, DamageEdge,
					AverageStrainVectorConcrete, AverageStressVectorConcrete, edge, Lchar );
				
				this->Set_NonConvergeddamages(DamageEdge, edge);
				DamagesOnEdges[edge] = DamageEdge;

			} // End loop over edges

            // Compute elemental damage
			double damage_element = this->CalculateElementalDamage(DamagesOnEdges);
			if (damage_element >= 0.999) { damage_element = 0.999; }
			this->Set_NonConvergeddamage(damage_element);
			
			const Vector& StressVectorConcrete = this->GetValue(CONCRETE_STRESS_VECTOR);
			IntegratedStressVectorConcrete = (1 - damage_element)*StressVectorConcrete;
			this->SetIntegratedStressVector(IntegratedStressVectorConcrete);

            // Linear elastic const matrix concrete
			Matrix ConstitutiveMatrixConcrete = ZeroMatrix(voigt_size, voigt_size);
			double Ec  = this->GetProperties()[YOUNG_MODULUS];
			double nuc = this->GetProperties()[POISSON_RATIO];
			this->CalculateConstitutiveMatrix(ConstitutiveMatrixConcrete, Ec, nuc);

            // Linear elastic const matrix steel
			Matrix ConstitutiveMatrixSteel = ZeroMatrix(voigt_size, voigt_size);
			double Es  = this->GetProperties()[YOUNG_MODULUS_STEEL];
			double nus = this->GetProperties()[POISSON_RATIO_STEEL];
			this->CalculateConstitutiveMatrix(ConstitutiveMatrixSteel, Es, nus);

            double k = this->GetProperties()[STEEL_VOLUMETRIC_PART];

            Matrix CompositeTangentMatrix = k*ConstitutiveMatrixSteel + (1.0-k)*(1.0-damage_element)*ConstitutiveMatrixConcrete;

			noalias(rLeftHandSideMatrix) += prod(trans(B), IntegrationWeight * Matrix(prod(CompositeTangentMatrix, B))); // LHS

			Vector VolumeForce = ZeroVector(dimension);
			VolumeForce = this->CalculateVolumeForce(VolumeForce, N);

			// RHS Volumetric load
			for (unsigned int i = 0; i < number_of_nodes; i++)
			{
				int index = dimension * i;
				for (unsigned int j = 0; j < dimension; j++)
				{
					rRightHandSideVector[index + j] += IntegrationWeight * N[i] * VolumeForce[j];
				}
			}

			//compute and add internal forces (RHS = rRightHandSideVector = Fext - Fint)
            Vector SteelStressVector = this->GetValue(STEEL_STRESS_VECTOR);
            Vector CompositeStressVector = k*SteelStressVector + (1.0-k)*IntegratedStressVectorConcrete;
			noalias(rRightHandSideVector) -= IntegrationWeight * prod(trans(B), CompositeStressVector);

			// Add nodal DEM forces
			Vector NodalRHS = ZeroVector(system_size);
			this->AddDEMContactForces(NodalRHS);
		
			// Add nodal contact forces from the DEM
			noalias(rRightHandSideVector) += NodalRHS;

		}
		KRATOS_CATCH("")
		//*****************************
	}

    Vector& RomFemDem3DElement::CalculateVolumeForce(Vector& rVolumeForce, const Vector& rN)
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		if (rVolumeForce.size() != dimension)
			rVolumeForce.resize(dimension, false);

		noalias(rVolumeForce) = ZeroVector(dimension);

		for (unsigned int j = 0; j < number_of_nodes; j++)
		{
			if (GetGeometry()[j].SolutionStepsDataHas(VOLUME_ACCELERATION)) { // it must be checked once at the begining only
				array_1d<double, 3 >& VolumeAcceleration = GetGeometry()[j].FastGetSolutionStepValue(VOLUME_ACCELERATION);
				for (unsigned int i = 0; i < dimension; i++)
					rVolumeForce[i] += rN[j] * VolumeAcceleration[i];
			}
		}
        double k = this->GetProperties()[STEEL_VOLUMETRIC_PART];
		rVolumeForce *= (GetProperties()[DENSITY]*(1.0-k) + GetProperties()[DENSITY_STEEL]*k);

		return rVolumeForce;

		KRATOS_CATCH("")
	}

    // 	TENSOR VARIABLES
	void RomFemDem3DElement::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo)
	{
    	const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    	const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

        if ( rOutput[0].size2() != dimension )
            rOutput[0].resize( dimension, dimension, false );

		if (rVariable == CONCRETE_STRESS_TENSOR)
		{
			rOutput[0] = MathUtils<double>::StressVectorToTensor(this->GetValue(CONCRETE_STRESS_VECTOR));
		}

		if (rVariable == STEEL_STRESS_TENSOR)
		{
			if (this->GetProperties()[STEEL_VOLUMETRIC_PART] > 0.0)
			{
				rOutput[0] = MathUtils<double>::StressVectorToTensor(this->GetValue(STEEL_STRESS_VECTOR));			
			}
			else
			{
				Matrix dummy;
				rOutput[0] = dummy;
			}
		}

		if (rVariable == STRAIN_TENSOR)
		{
			rOutput[0] =  MathUtils<double>::StrainVectorToTensor(this->GetValue(STRAIN_VECTOR));
		}

		if (rVariable == CONCRETE_STRESS_TENSOR_INTEGRATED)
		{
			rOutput[0] =  MathUtils<double>::StressVectorToTensor(this->GetIntegratedStressVector());
		}

    }

	// Tensor variables
	void RomFemDem3DElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rValues,
		const ProcessInfo& rCurrentProcessInfo )
	{
		if (rVariable == STRAIN_TENSOR)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
		if (rVariable == STEEL_STRESS_TENSOR)
		{

			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
		if (rVariable == CONCRETE_STRESS_TENSOR)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
		if (rVariable == CONCRETE_STRESS_TENSOR_INTEGRATED)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
	}






} // Element