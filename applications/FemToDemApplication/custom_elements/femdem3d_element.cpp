#include "includes/define.h"
#include <string>
#include "includes/constitutive_law.h"
#include "custom_constitutive/zarate_law.hpp"
#include "femdem3d_element.hpp"
#include "includes/element.h"
#include "includes/node.h"
#include "fem_to_dem_application_variables.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"
#include "solid_mechanics_application_variables.h"
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos
{
	//***********************DEFAULT CONSTRUCTOR******************************************
	//************************************************************************************

	FemDem3DElement::FemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: SmallDisplacementElement(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}
	//******************************CONSTRUCTOR*******************************************
	//************************************************************************************

	FemDem3DElement::FemDem3DElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		:SmallDisplacementElement(NewId, pGeometry, pProperties)
	{
		//BY DEFAULT, THE GEOMETRY WILL DEFINE THE INTEGRATION METHOD
		mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
	}

	//******************************COPY CONSTRUCTOR**************************************
	//************************************************************************************

	FemDem3DElement::FemDem3DElement(FemDem3DElement const& rOther)
		:SmallDisplacementElement(rOther)
	{
		//ALL MEMBER VARIABLES THAT MUST BE KEPT AFTER COPYING AN ELEMENT HAVE TO BE DEFINED HERE
		//IF NO ASSIGMENT OPERATOR IS DEFINED THE COPY CONSTRUCTOR WILL DEFINE IT BY DEFFAULT
	}

	//*******************************ASSIGMENT OPERATOR***********************************
	//************************************************************************************

	FemDem3DElement&  FemDem3DElement::operator=(FemDem3DElement const& rOther)
	{
		//ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

		SmallDisplacementElement::operator=(rOther);
		return *this;
	}

	//*********************************OPERATIONS*****************************************
	//************************************************************************************

	Element::Pointer FemDem3DElement::Create(IndexType NewId, NodesArrayType const& rThisNodes, PropertiesType::Pointer pProperties) const
	{
		//NEEDED TO CREATE AN ELEMENT   
		return Element::Pointer(new FemDem3DElement(NewId, GetGeometry().Create(rThisNodes), pProperties));
	}


	//************************************CLONE*******************************************
	//************************************************************************************

	Element::Pointer FemDem3DElement::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
	{

		//YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
		//ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

		FemDem3DElement NewElement(NewId, GetGeometry().Create(rThisNodes), pGetProperties());

		return Element::Pointer(new FemDem3DElement(NewElement));
	}

	//*******************************DESTRUCTOR*******************************************
	//************************************************************************************

	FemDem3DElement::~FemDem3DElement()
	{
	}

	void FemDem3DElement::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
	{
		//Geometry< Node < 3 > >& Nodes = this->GetGeometry();

		//Node<3>& pNode = Nodes[0];
		//WeakPointerVector< Element >& rneigh_el = pNode.GetValue(NEIGHBOUR_ELEMENTS); //elementos vecinos del nodo
		//KRATOS_WATCH(Nodes[0].Id())
		////KRATOS_WATCH(rneigh_el.size())
		//KRATOS_WATCH(rneigh_el[0].Id())
		//KRATOS_WATCH(rneigh_el[1].Id())
		//KRATOS_WATCH(rneigh_el[2].Id())
		//KRATOS_WATCH(rneigh_el[3].Id())
		//KRATOS_WATCH(rneigh_el[3].Id())

		//std::vector<Element> aux;
		//std::vector <std::vector<Element>> AUX;
		//aux.push_back(rneigh_el[2]);
		//AUX.push_back(aux);

		//mEdgeNeighboursContainer.push_back(aux);

		this->ComputeEdgeNeighbours(rCurrentProcessInfo);


	}

	void FemDem3DElement::ComputeEdgeNeighbours(ProcessInfo& rCurrentProcessInfo)
	{
		std::vector <std::vector<Element*>> EdgeNeighboursContainer;

		Geometry< Node < 3 > >& NodesCurrentElement = this->GetGeometry();

		Node<3>& pNode0 = NodesCurrentElement[0];
		Node<3>& pNode1 = NodesCurrentElement[1];
		Node<3>& pNode2 = NodesCurrentElement[2];
		Node<3>& pNode3 = NodesCurrentElement[3];
		
		// Neighbour elements of each node of the current element 
		WeakPointerVector< Element >& NeighNode0 = pNode0.GetValue(NEIGHBOUR_ELEMENTS);
		WeakPointerVector< Element >& NeighNode1 = pNode1.GetValue(NEIGHBOUR_ELEMENTS);
		WeakPointerVector< Element >& NeighNode2 = pNode2.GetValue(NEIGHBOUR_ELEMENTS);
		WeakPointerVector< Element >& NeighNode3 = pNode3.GetValue(NEIGHBOUR_ELEMENTS);

		// Nodal neighbours container
		std::vector<WeakPointerVector< Element >> NodalNeighbours;
		NodalNeighbours.push_back(NeighNode0);
		NodalNeighbours.push_back(NeighNode1);
		NodalNeighbours.push_back(NeighNode2);
		NodalNeighbours.push_back(NeighNode3);

		// Aux indexes
		Matrix NodeIndexes = ZeroMatrix(6, 2);
		this->SetNodeIndexes(NodeIndexes);

		// Loop over EDGES to assign the elements that share that edge -> Fill mEdgeNeighboursContainer
		for (int edge = 0 ; edge < 4; edge++)
		{
			int NodeIndex1 = NodeIndexes(edge, 0);
			int NodeIndex2 = NodeIndexes(edge, 1);

			// Neigh elements of local node 1 and 2  // 
			WeakPointerVector< Element >& NeighOfNode1 = NodalNeighbours[NodeIndex1];
			WeakPointerVector< Element >& NeighOfNode2 = NodalNeighbours[NodeIndex2];

			int NodeId1 = NodesCurrentElement[NodeIndex1].Id();
			int NodeId2 = NodesCurrentElement[NodeIndex2].Id();

			std::vector<Element*> EdgeSharedElementsNode1;
			// Loop over neigh elements of the node 1
			for (int neigh_elem = 0; neigh_elem < NeighOfNode1.size();neigh_elem++)
			{
				//std::vector<Element> EdgeSharedElementsNode1;

				// Nodes of the neigh element
				Geometry< Node < 3 > >& NodesNeighElem = NeighOfNode1[neigh_elem].GetGeometry();

				// Loop over the nodes of the neigh element
				for (int neigh_elem_node = 0; neigh_elem_node < 4 ; neigh_elem_node++)
				{
					int NeighElementNodeId = NodesNeighElem[neigh_elem_node].Id();

					if (NeighElementNodeId == NodeId2 & this->Id() != NeighOfNode1[neigh_elem].Id())
					{
						EdgeSharedElementsNode1.push_back(&NeighOfNode1[neigh_elem]); // ( [] returns a Element object!!)
					} 
				}
			}

			std::vector<Element*> EdgeSharedElementsNode2;
			// Loop over neigh elements of the node 2
			for (int neigh_elem = 0; neigh_elem < NeighOfNode2.size();neigh_elem++)
			{
				//std::vector<Element> EdgeSharedElementsNode2;

				// Nodes of the neigh element
				Geometry< Node < 3 > >& NodesNeighElem = NeighOfNode2[neigh_elem].GetGeometry();

				// Loop over the nodes of the neigh element
				for (int neigh_elem_node = 0; neigh_elem_node < 4 ;neigh_elem_node++)
				{
					int NeighElementNodeId = NodesNeighElem[neigh_elem_node].Id();

					if (NeighElementNodeId == NodeId1 & this->Id() != NeighOfNode2[neigh_elem].Id())
					{
						EdgeSharedElementsNode2.push_back(&NeighOfNode2[neigh_elem]);
					}
				}
			}

			// Let's create the vector of neighbour elements for this edge
			std::vector<Element*> EdgeSharedElements = EdgeSharedElementsNode1;

			// Add the neigh elements from the node 2
			for (int i = 0; i < EdgeSharedElementsNode2.size(); i++)
			{
				int aux = 0;

				for (int j = 0; j < EdgeSharedElements.size(); j++)
				{
					if (EdgeSharedElementsNode2[i]->Id() == EdgeSharedElements[j]->Id()) aux++;
				}

				if (aux == 0) EdgeSharedElements.push_back(EdgeSharedElementsNode2[i]);
			}

			EdgeNeighboursContainer.push_back(EdgeSharedElements);



			KRATOS_WATCH(this->Id())
			KRATOS_WATCH(NodeId1)
			KRATOS_WATCH(NodeId2)
			//KRATOS_WATCH(EdgeSharedElementsNode1[0].Id())
			//KRATOS_WATCH(EdgeSharedElementsNode2[0].Id())
			KRATOS_WATCH(EdgeSharedElementsNode1.size())
			KRATOS_WATCH(EdgeSharedElementsNode2.size())
			//KRATOS_WATCH(EdgeSharedElements.size())

			for (int i=0;i<EdgeSharedElements.size();i++)
			{
				KRATOS_WATCH(EdgeSharedElements[i]->Id())
			}

		} // End loop edges

		// Storages the information inside the element
		this->SaveEdgeNeighboursContainer(EdgeNeighboursContainer);
	}

	void FemDem3DElement::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
	{
		//double CurrentfSigma = 0.0, damage_element = 0.0;
		

		// //Loop over edges
		// for (int cont = 0; cont < 3; cont++)
		// {
		// 	this->Set_Convergeddamages(this->Get_NonConvergeddamages(cont), cont);
		// 	this->SetConverged_f_sigmas(this->Get_NonConvergedf_sigma(cont), cont);
		// 	CurrentfSigma = this->GetConverged_f_sigmas(cont);
		// 	if (CurrentfSigma > this->Get_threshold(cont)) { this->Set_threshold(CurrentfSigma, cont); }
		// } // End Loop over edges

		// damage_element = this->Get_NonConvergeddamage();
		// this->Set_Convergeddamage(damage_element);

		// if (damage_element > 0.0) 
		// {
		// 	this->SetValue(IS_DAMAGED, 1);
		// }
		
		// if (damage_element >= 0.98)
		// {
		// 	this->Set(ACTIVE, false);
		// 	double old_threshold = this->GetValue(STRESS_THRESHOLD);
		// 	this->SetValue(INITIAL_THRESHOLD, old_threshold);
		// }

		// this->ResetNonConvergedVars();
		// this->SetToZeroIteration();

		// // computation of the equivalent damage threshold and damage of the element for AMR mapping
		// Vector thresholds = this->GetThresholds();
		
		// Vector TwoMinValues;
		// this->Get2MaxValues(TwoMinValues, thresholds[0], thresholds[1], thresholds[2]);  // todo ojo con la funcion modificada
		// double EqThreshold = 0.5*(TwoMinValues[0] + TwoMinValues[1]);  // El menor o mayor?? TODO


		// this->SetValue(STRESS_THRESHOLD, EqThreshold); // AMR
		// this->Set_threshold(EqThreshold);
		// this->SetValue(DAMAGE_ELEMENT, damage_element);


		// // Reset the nodal force flag for the next time step
		// Geometry< Node < 3 > >& NodesElement = this->GetGeometry();
		// for (int i = 0; i < 3; i++)
		// {
		// 	#pragma omp critical
		// 	{
		// 		NodesElement[i].SetValue(NODAL_FORCE_APPLIED, false);
		// 	}
		// }
	}

	

	void FemDem3DElement::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
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


				if (detJ<0) KRATOS_THROW_ERROR(std::invalid_argument, " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ)

				//compute cartesian derivatives for this integration point  [dN/dx_n]
				noalias(DN_DX) = prod(DN_De[PointNumber], InvJ);

				//set shape functions for this integration point
				Vector N = row(Ncontainer, PointNumber);


				//b.-compute infinitessimal strain
				this->CalculateInfinitesimalStrain(StrainVector, DN_DX);
				this->SetValue(STRAIN_VECTOR, StrainVector);

				ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

				//set constitutive law variables: (it passes only references to this local variables)
				Values.SetStrainVector(StrainVector);
				Values.SetStressVector(StressVector);
				Values.SetConstitutiveMatrix(ConstitutiveMatrix);
				Values.SetShapeFunctionsDerivatives(DN_DX);
				Values.SetShapeFunctionsValues(N);
				//values to be set:
				Values.SetDeterminantF(detF);
				Values.SetDeformationGradientF(F);

				//set constitutive law flags:
				Flags &ConstitutiveLawOptions = Values.GetOptions();

				//compute stress and constitutive matrix
				ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS);
				ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

				//CALL THE CONSTITUTIVE LAW (for this integration point)
				//(after calling the constitutive law StressVector and ConstitutiveMatrix are set and can be used)
				mConstitutiveLawVector[PointNumber]->CalculateMaterialResponseCauchy(Values);  
				this->SetValue(STRESS_VECTOR, Values.GetStressVector());

				this->CalculateDeformationMatrix(B, DN_DX);
				this->SetBMatrix(B);

			}
		}
		
		KRATOS_CATCH("")

	}

	
	void FemDem3DElement::CalculateLocalSystem (MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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

			Vector IntegratedStressVector = ZeroVector(voigt_size);

			// Find Neighbour Elements
			WeakPointerVector< Element >& elem_neigb = this->GetValue(NEIGHBOUR_ELEMENTS);
			if (elem_neigb.size() == 0) { KRATOS_THROW_ERROR(std::invalid_argument, " Neighbour Elements not calculated --> size = ", elem_neigb.size()) }

			const Vector& StressVector = this->GetValue(STRESS_VECTOR);
			
			Matrix ConstitutiveMatrix = ZeroMatrix(voigt_size, voigt_size);
			double E  = this->GetProperties()[YOUNG_MODULUS];
			double nu = this->GetProperties()[POISSON_RATIO];
			this->CalculateConstitutiveMatrix(ConstitutiveMatrix, E, nu);

			noalias(rLeftHandSideMatrix) += prod(trans(B), IntegrationWeight * Matrix(prod(ConstitutiveMatrix, B))); //LHS

			Vector VolumeForce = ZeroVector(dimension);
			VolumeForce = this->CalculateVolumeForce(VolumeForce, N);

			// RHS
			for (unsigned int i = 0; i < number_of_nodes; i++)
			{
				int index = dimension * i;
				for (unsigned int j = 0; j < dimension; j++)
				{
					rRightHandSideVector[index + j] += IntegrationWeight * N[i] * VolumeForce[j];
				}
			}

			//compute and add internal forces (RHS = rRightHandSideVector = Fext - Fint)
			noalias(rRightHandSideVector) -= IntegrationWeight * prod(trans(B), StressVector);

		}
		KRATOS_CATCH("")
		//*****************************
	}

	void FemDem3DElement::CalculateDeformationMatrix(Matrix& rB, const Matrix& rDN_DX)
	{
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int voigt_size = dimension * (dimension + 1) * 0.5;

		/*if (rB.size1() != voigt_size || rB.size2() != dimension*number_of_nodes)
			rB.resize(voigt_size, dimension*number_of_nodes, false);*/

		for (unsigned int i = 0; i < number_of_nodes; i++)
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

	void FemDem3DElement::CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix, const double &rYoungModulus,
		const double &rPoissonCoefficient)
	{
		rConstitutiveMatrix.clear();

		// 3D linear elastic constitutive matrix
		rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus*(1.0-rPoissonCoefficient)/((1.0+rPoissonCoefficient)*(1.0-2.0*rPoissonCoefficient)));
		rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );
		rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 );

		rConstitutiveMatrix ( 3 , 3 ) = rConstitutiveMatrix ( 0 , 0 )*(1.0-2.0*rPoissonCoefficient)/(2.0*(1.0-rPoissonCoefficient));
		rConstitutiveMatrix ( 4 , 4 ) = rConstitutiveMatrix ( 3 , 3 );
		rConstitutiveMatrix ( 5 , 5 ) = rConstitutiveMatrix ( 3 , 3 );

		rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient/(1.0-rPoissonCoefficient);
		rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

		rConstitutiveMatrix ( 0 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
		rConstitutiveMatrix ( 2 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

		rConstitutiveMatrix ( 1 , 2 ) = rConstitutiveMatrix ( 0 , 1 );
		rConstitutiveMatrix ( 2 , 1 ) = rConstitutiveMatrix ( 0 , 1 );
	}

	void FemDem3DElement::CalculateDN_DX(Matrix& rDN_DX, int PointNumber)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		//reading integration points
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		//get the shape functions [N] (for the order of the default integration method)
		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues(mThisIntegrationMethod);

		//get the shape functions parent coordinates derivative [dN/d�] (for the order of the default integration method)
		const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients(mThisIntegrationMethod);
		//calculate delta position (here coincides with the current displacement)
		Matrix DeltaPosition = ZeroMatrix(number_of_nodes, dimension);
		DeltaPosition = this->CalculateDeltaPosition(DeltaPosition);
		//KRATOS_WATCH(DeltaPosition)
		//calculating the reference jacobian from cartesian coordinates to parent coordinates for all integration points [dx_n/d�]
		GeometryType::JacobiansType J;
		J.resize(1, false);
		J[0] = ZeroMatrix(1, 1);
		J = GetGeometry().Jacobian(J, mThisIntegrationMethod, DeltaPosition);
		//a.-compute element kinematics

		//calculating the inverse of the jacobian for this integration point[d�/dx_n]
		Matrix InvJ = ZeroMatrix(dimension, dimension);
		double detJ = 0;
		MathUtils<double>::InvertMatrix(J[PointNumber], InvJ, detJ);

		if (detJ < 0)
			KRATOS_THROW_ERROR(std::invalid_argument, " SMALL DISPLACEMENT ELEMENT INVERTED: |J|<0 ) detJ = ", detJ)

			//compute cartesian derivatives for this integration point  [dN/dx_n]
			rDN_DX = prod(DN_De[PointNumber], InvJ);
	}

	void FemDem3DElement::CalculateInfinitesimalStrain(Vector& rStrainVector, const Matrix& rDN_DX)
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		Matrix H = zero_matrix<double>(dimension); //[dU/dx_n]
		

		for (unsigned int i = 0; i < number_of_nodes; i++)
		{

			array_1d<double, 3 > & Displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);

            H ( 0 , 0 ) += Displacement[0]*rDN_DX ( i , 0 );
            H ( 0 , 1 ) += Displacement[0]*rDN_DX ( i , 1 );
            H ( 0 , 2 ) += Displacement[0]*rDN_DX ( i , 2 );
            H ( 1 , 0 ) += Displacement[1]*rDN_DX ( i , 0 );
            H ( 1 , 1 ) += Displacement[1]*rDN_DX ( i , 1 );
            H ( 1 , 2 ) += Displacement[1]*rDN_DX ( i , 2 );
            H ( 2 , 0 ) += Displacement[2]*rDN_DX ( i , 0 );
            H ( 2 , 1 ) += Displacement[2]*rDN_DX ( i , 1 );
            H ( 2 , 2 ) += Displacement[2]*rDN_DX ( i , 2 );
		}

		//Infinitesimal Strain Calculation
		if (rStrainVector.size() != 6) rStrainVector.resize(6, false);

        rStrainVector[0] = H( 0, 0 );

        rStrainVector[1] = H( 1, 1 );

        rStrainVector[2] = H( 2, 2 );

        rStrainVector[3] = ( H( 0, 1 ) + H( 1, 0 ) ); // xy

        rStrainVector[4] = ( H( 1, 2 ) + H( 2, 1 ) ); // yz

        rStrainVector[5] = ( H( 0, 2 ) + H( 2, 0 ) ); // xz
		
		KRATOS_CATCH("")
	}

	void FemDem3DElement::CalculateStressVector(Vector& rStressVector, const Matrix& rConstitutiveMAtrix, const Vector& rInfinitesimalStrainVector)
	{
		noalias(rStressVector) = prod(rConstitutiveMAtrix, rInfinitesimalStrainVector);
	}

	void FemDem3DElement::CalculatePrincipalStress(Vector& PrincipalStressVector, const Vector StressVector)
	{
		PrincipalStressVector.resize(2);
		PrincipalStressVector[0] = 0.5*(StressVector[0] + StressVector[1]) + sqrt(pow(0.5*(StressVector[0] - StressVector[1]), 2) + pow(StressVector[2], 2));
		PrincipalStressVector[1] = 0.5*(StressVector[0] + StressVector[1]) - sqrt(pow(0.5*(StressVector[0] - StressVector[1]), 2) + pow(StressVector[2], 2));
	}

	void FemDem3DElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{
	}

	void FemDem3DElement::AverageVector(Vector& rAverageVector, const Vector& v, const Vector& w)
	{
		int n = v.size();
		int m = w.size();
		if (n != m) KRATOS_ERROR << "The dimension of the vectors are different or null";
		rAverageVector.resize(n);
		for (int cont = 0;cont < n;cont++)
		{
			rAverageVector[cont] = (v[cont] + w[cont])*0.5;
		}
	}

	// Double values
	void FemDem3DElement::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == DAMAGE_ELEMENT || rVariable == IS_DAMAGED || rVariable == STRESS_THRESHOLD)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}

		
	}

	// Vector Values
	void FemDem3DElement::GetValueOnIntegrationPoints(const Variable<Vector>& rVariable,
		std::vector<Vector>& rValues,
		const ProcessInfo& rCurrentProcessInfo)
	{

		if (rVariable == STRAIN_VECTOR)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}

		if (rVariable == STRESS_VECTOR)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}

		if (rVariable == STRESS_VECTOR_INTEGRATED)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}

		
	}

	// Tensor variables
	void FemDem3DElement::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable,
		std::vector<Matrix>& rValues,
		const ProcessInfo& rCurrentProcessInfo )

	{
		if (rVariable == STRAIN_TENSOR)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}

		if (rVariable == STRESS_TENSOR)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}

		if (rVariable == STRESS_TENSOR_INTEGRATED)
		{
			CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
		}
	}

	// DOUBLE VARIABLES
	void FemDem3DElement::CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == DAMAGE_ELEMENT)
		{
			rOutput.resize(1);
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
				//rOutput[PointNumber] = double(this->Get_Convergeddamage());
				rOutput[PointNumber] = double(this->GetValue(DAMAGE_ELEMENT));
			}
		}

		if (rVariable == IS_DAMAGED)
		{
			rOutput.resize(1);
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
				rOutput[PointNumber] = double(this->GetValue(IS_DAMAGED));
			}
		}

		if (rVariable == STRESS_THRESHOLD)
		{
			rOutput.resize(1);
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) {
				rOutput[PointNumber] = double(this->GetValue(STRESS_THRESHOLD));
			}
		}
	}

	// 	VECTOR VARIABLES
	void FemDem3DElement::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo)
	{

		KRATOS_TRY

		if (rVariable == STRESS_VECTOR)
		{
			rOutput[0] = this->GetValue(STRESS_VECTOR);
		}

		if (rVariable == STRAIN_VECTOR)
		{
			rOutput[0] = this->GetValue(STRAIN_VECTOR);
		}

		if (rVariable == STRESS_VECTOR_INTEGRATED)
		{
			rOutput[0] = this->GetIntegratedStressVector();
		}

		KRATOS_CATCH("")
		
	}


	// 	TENSOR VARIABLES
	void FemDem3DElement::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo)
	{
    	const unsigned int& integration_points_number = GetGeometry().IntegrationPointsNumber( mThisIntegrationMethod );
    	const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

        if ( rOutput[0].size2() != dimension )
            rOutput[0].resize( dimension, dimension, false );

		if (rVariable == STRESS_TENSOR)
		{
			rOutput[0] = MathUtils<double>::StressVectorToTensor(this->GetValue(STRESS_VECTOR));
		}

		if (rVariable == STRAIN_TENSOR)
		{
			rOutput[0] =  MathUtils<double>::StrainVectorToTensor(this->GetValue(STRAIN_VECTOR));
		}

		if (rVariable == STRESS_TENSOR_INTEGRATED)
		{
			rOutput[0] =  MathUtils<double>::StressVectorToTensor(this->GetIntegratedStressVector());
		}

	}


	double FemDem3DElement::CalculateLchar(FemDem3DElement* CurrentElement, const Element& NeibElement, int cont)
	{
		Geometry< Node < 3 > >& NodesElem1 = CurrentElement->GetGeometry();  // 3 nodes of the Element 1
		Geometry< Node < 3 > > NodesElem2  = NeibElement.GetGeometry();      // "         " 2
		Vector Xcoord, Ycoord;
		Xcoord.resize(3);
		Ycoord.resize(3);

		// Let's find the two shared nodes between the 2 elements 
		int aux = 0;
		double l_char = 0;
		for (int cont = 0; cont < 3; cont++)
		{
			for (int cont2 = 0; cont2 < 3; cont2++)
			{
				if (NodesElem1[cont].Id() == NodesElem2[cont2].Id())
				{
					Xcoord[aux] = NodesElem1[cont].X0();
					Ycoord[aux] = NodesElem1[cont].Y0();
					aux++;                              // aux > 3 if the two elements are the same one (in fact aux == 9)
				}
			}
		} // End finding nodes

		// Computation of the l_char
		if (aux < 3) 
		{                                                                                        // It is not an edge element --> The 2 elements are not equal
			l_char = pow((pow(Xcoord[0] - Xcoord[1], 2) + pow(Ycoord[0] - Ycoord[1], 2)), 0.5);  // Length of the edge between 2 elements
			//l_char = length;                                                                   // Currently the characteristic length is the edge length (can be modified)
		}
		else  // Edge Element
		{ 
			double ElementArea = abs(this->GetGeometry().Area());
			l_char = sqrt(4 * ElementArea / sqrt(3));   // Cervera's Formula
			
		} // l_char computed

		CurrentElement->Set_l_char(l_char, cont);  // Storages the l_char of this side
		CurrentElement->IterationPlus();

		return 0.0;
	}


	void FemDem3DElement::Get2MaxValues(Vector& MaxValues, double a, double b, double c)
	{
		MaxValues.resize(2);
		Vector V;
		V.resize(3);
		V[0] = a;
		V[1] = b;
		V[2] = c;
		int n = 3, imin = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - 1; j++) {
				if (V[j] > V[j + 1]) {
					double aux = V[j];
					V[j] = V[j + 1];
					V[j + 1] = aux;
				}
			}
		}
		MaxValues[0] = V[2];
		MaxValues[1] = V[1];
	}

	void FemDem3DElement::Get2MinValues(Vector& MaxValues, double a, double b, double c)
	{
		MaxValues.resize(2);
		Vector V;
		V.resize(3);
		V[0] = a;
		V[1] = b;
		V[2] = c;
		int n = 3, imin = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - 1; j++) {
				if (V[j] > V[j + 1]) {
					double aux = V[j];
					V[j] = V[j + 1];
					V[j + 1] = aux;
				}
			}
		}
		MaxValues[0] = V[1];
		MaxValues[1] = V[0];
	}
	double FemDem3DElement::Calculate_I1_Invariant(double sigma1, double sigma2) { return sigma1 + sigma2; }
	double FemDem3DElement::Calculate_J2_Invariant(double sigma1, double sigma2)
	{
		return  (pow((sigma1 - sigma2), 2) + pow(sigma1, 2) + pow(sigma2, 2)) / 6;
	}
	double FemDem3DElement::Calculate_J3_Invariant(double sigma1, double sigma2, double I1)
	{
		return	(sigma1 - I1 / 3)*((sigma2 - I1 / 3))*(-I1 / 3);
	}
	double FemDem3DElement::Calculate_Theta_Angle(double J2, double J3)
	{
		double sint3;
		sint3 = (-3.0*sqrt(3)*J3) / (2 * J2*sqrt(J2));
		if (sint3 < -0.95) { sint3 = -1; }
		if (sint3 > 0.95) { sint3 = 1; }
		return asin(sint3) / 3;
	}

	void FemDem3DElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			bool ComputeLumpedMassMatrix = false;
		if (rCurrentProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX))
			if (rCurrentProcessInfo[COMPUTE_LUMPED_MASS_MATRIX] == true)
				ComputeLumpedMassMatrix = true;

		if (ComputeLumpedMassMatrix == false) {

			//create local system components
			LocalSystemComponents LocalSystem;

			//calculation flags   
			LocalSystem.CalculationFlags.Set(SmallDisplacementElement::COMPUTE_LHS_MATRIX);

			VectorType RightHandSideVector = Vector();

			//Initialize sizes for the system components:
			this->InitializeSystemMatrices(rMassMatrix, RightHandSideVector, LocalSystem.CalculationFlags);

			//Set Variables to Local system components
			LocalSystem.SetLeftHandSideMatrix(rMassMatrix);
			LocalSystem.SetRightHandSideVector(RightHandSideVector);

			//Calculate elemental system
			CalculateDynamicSystem(LocalSystem, rCurrentProcessInfo);

		}
		else {

			//lumped
			unsigned int dimension = GetGeometry().WorkingSpaceDimension();
			const unsigned int number_of_nodes = GetGeometry().PointsNumber();
			unsigned int MatSize = dimension * number_of_nodes;

			if (rMassMatrix.size1() != MatSize)
				rMassMatrix.resize(MatSize, MatSize, false);

			noalias(rMassMatrix) = ZeroMatrix(MatSize, MatSize);

			double TotalMass = 0;
			TotalMass = this->CalculateTotalMass(TotalMass, rCurrentProcessInfo);

			Vector LumpFact(number_of_nodes);
			noalias(LumpFact) = ZeroVector(number_of_nodes);

			LumpFact = GetGeometry().LumpingFactors(LumpFact);

			for (unsigned int i = 0; i < number_of_nodes; i++)
			{
				double temp = LumpFact[i] * TotalMass;

				for (unsigned int j = 0; j < dimension; j++)
				{
					unsigned int index = i * dimension + j;
					rMassMatrix(index, index) = temp;
				}
			}

		}

		KRATOS_CATCH("")
	}
	Vector& FemDem3DElement::CalculateVolumeForce(Vector& rVolumeForce, const Vector& rN)
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

		rVolumeForce *= GetProperties()[DENSITY];

		return rVolumeForce;

		KRATOS_CATCH("")
	}

	double FemDem3DElement::GetMaxValue(Vector Strain)
	{
		Vector V;
		int n = Strain.size();
		V.resize(n);

		for (int cont = 0;cont < n;cont++)
		{
			V[cont] = Strain[cont];
		}

		int imin = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - 1; j++) {
				if (V[j] > V[j + 1]) {
					double aux = V[j];
					V[j] = V[j + 1];
					V[j + 1] = aux;
				}
			}
		}

		return V[n - 1];
	}

	double FemDem3DElement::GetMaxAbsValue(Vector Strain)
	{
		Vector V;
		int n = Strain.size();
		V.resize(n);

		for (int cont = 0;cont < n;cont++)
		{
			V[cont] = abs(Strain[cont]);
		}

		int imin = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - 1; j++) {
				if (V[j] > V[j + 1]) {
					double aux = V[j];
					V[j] = V[j + 1];
					V[j + 1] = aux;
				}
			}
		}

		return V[n - 1];
	}

	double FemDem3DElement::GetMinAbsValue(Vector Strain)
	{
		Vector V;
		V.resize(3);
		V[0] = abs(Strain[0]);
		V[1] = abs(Strain[1]);
		V[2] = abs(Strain[2]);
		int n = 3, imin = 0;

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n - 1; j++) {
				if (V[j] > V[j + 1]) {
					double aux = V[j];
					V[j] = V[j + 1];
					V[j + 1] = aux;
				}
			}
		}

		return V[0];
	}


	// ****** Tangent Constitutive Tensor by Numerical Derivation ******
	void FemDem3DElement::PerturbateStrainComponent(const Vector& rStrainVector, Vector& PertubatedStrain, const double& perturbation, int component)
	{
		PertubatedStrain = rStrainVector;
		PertubatedStrain[component] += perturbation;
	}

	double FemDem3DElement::CalculatePerturbation(const Vector& StrainVector, int component)
	{
		double Pert = 0.0;
		if (StrainVector[component] != 0.0) Pert = (1e-5)*StrainVector[component];
		else Pert = (1e-5)*GetMinAbsValue(StrainVector);
		if (Pert < this->GetMaxAbsValue(StrainVector)*(1e-10)) { Pert = GetMaxAbsValue(StrainVector)*(1e-10); }

		return Pert;
	}

	void  FemDem3DElement::CalculateTangentTensor(Matrix& rTangentTensor, const Vector& StrainVector, 
		const Vector& IntegratedStressVector, int cont, double l_char)
	{
		rTangentTensor.resize(3, 3);

		double E =  this->GetProperties()[YOUNG_MODULUS];
		double nu = this->GetProperties()[POISSON_RATIO];
		Matrix ConstitutiveMatrix = ZeroMatrix(3, 3);
		this->CalculateConstitutiveMatrix(ConstitutiveMatrix, E, nu);

		// Perturbed Strain Vectors
		Vector PerturbedStrain1 = ZeroVector(3);
		Vector PerturbedStrain2 = ZeroVector(3);
		Vector PerturbedStrain3 = ZeroVector(3);

		Vector PerturbedIntegratedStress1 = ZeroVector(3);
		Vector PerturbedIntegratedStress2 = ZeroVector(3);
		Vector PerturbedIntegratedStress3 = ZeroVector(3);

		Vector PerturbedStress1 = ZeroVector(3);
		Vector PerturbedStress2 = ZeroVector(3);
		Vector PerturbedStress3 = ZeroVector(3);

		Vector DeltaStress1 = ZeroVector(3);
		Vector DeltaStress2 = ZeroVector(3);
		Vector DeltaStress3 = ZeroVector(3);

		// Calculation of the perturbations
		double Perturbation1 = this->CalculatePerturbation(StrainVector, 0);
		double Perturbation2 = this->CalculatePerturbation(StrainVector, 1);
		double Perturbation3 = this->CalculatePerturbation(StrainVector, 2);

		// Calculation of the Perturbed Strain Vectors
		this->PerturbateStrainComponent(StrainVector, PerturbedStrain1, Perturbation1, 0);
		this->PerturbateStrainComponent(StrainVector, PerturbedStrain2, Perturbation2, 1);
		this->PerturbateStrainComponent(StrainVector, PerturbedStrain3, Perturbation3, 2);

		// Calculation of the Perturbed Predictive Stress Vectors
		this->CalculateStressVector(PerturbedStress1, ConstitutiveMatrix, PerturbedStrain1);
		this->CalculateStressVector(PerturbedStress2, ConstitutiveMatrix, PerturbedStrain2);
		this->CalculateStressVector(PerturbedStress3, ConstitutiveMatrix, PerturbedStrain3);

		// Integration of the Perturbed Predictive Stress Vectors
		double damage1 = 0.0, damage2 = 0.0, damage3 = 0.0;

		this->TangentModifiedMohrCoulombCriterion(PerturbedIntegratedStress1, damage1,  PerturbedStress1, cont, l_char);
		this->TangentModifiedMohrCoulombCriterion(PerturbedIntegratedStress2, damage2,  PerturbedStress2, cont, l_char);
		this->TangentModifiedMohrCoulombCriterion(PerturbedIntegratedStress3, damage3,  PerturbedStress3, cont, l_char);

		DeltaStress1 = PerturbedIntegratedStress1 - IntegratedStressVector;
		DeltaStress2 = PerturbedIntegratedStress2 - IntegratedStressVector;
		DeltaStress3 = PerturbedIntegratedStress3 - IntegratedStressVector;

		for (int row = 0;row < 3;row++) // DeltaStress is the i column of the Tangent Tensor
		{
			rTangentTensor(row, 0) = DeltaStress1[row] / (Perturbation1);
			rTangentTensor(row, 1) = DeltaStress2[row] / (Perturbation2);
			rTangentTensor(row, 2) = DeltaStress3[row] / (Perturbation3);
		}
	}

	void FemDem3DElement::TangentModifiedMohrCoulombCriterion(Vector& rIntegratedStress, double& damage, const Vector& StressVector, int cont, double l_char)
	{
		rIntegratedStress.resize(3);
		Vector PrincipalStressVector = ZeroVector(2);
		this->CalculatePrincipalStress(PrincipalStressVector, StressVector);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0;
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		friction_angle = this->GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159 / 180; // In radians!
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];

		// Check input variables 
		if (friction_angle < 1e-24) { friction_angle = 32 * 3.14159 / 180; std::cout << "Friction Angle not defined, assumed equal to 32� " << std::endl; }
		if (sigma_c < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa "; }
		if (sigma_t < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa "; }
		if (Gt < 1e-24) { KRATOS_ERROR << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa "; }

		double K1, K2, K3, Rmorh, R, alpha_r, c_max, theta, c_threshold;
		R = abs(sigma_c / sigma_t);
		Rmorh = pow(tan((3.14159 / 4) + friction_angle / 2), 2);
		alpha_r = R / Rmorh;
		c_max = abs(sigma_c);

		double I1, J2, J3;
		I1 = Calculate_I1_Invariant(PrincipalStressVector[0], PrincipalStressVector[1]);
		J2 = Calculate_J2_Invariant(PrincipalStressVector[0], PrincipalStressVector[1]);
		J3 = Calculate_J3_Invariant(PrincipalStressVector[0], PrincipalStressVector[1], I1);
		K1 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r)*sin(friction_angle);
		K2 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r) / sin(friction_angle);
		K3 = 0.5*(1 + alpha_r)*sin(friction_angle) - 0.5*(1 - alpha_r);

		double n = sigma_c / sigma_t;
		//double ElementArea = this->GetGeometry().Area();
		//double l_char = sqrt(4 * ElementArea / sqrt(3));

		double A = 1.00 / (n*n*Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		double f = 0.0, F = 0.0; /// F = f-c = 0 classical definition of yield surface

		// Check Modified Mohr-Coulomb criterion
		if (PrincipalStressVector[0] == 0 && PrincipalStressVector[1] == 0) { f = 0; }
		else
		{
			theta = Calculate_Theta_Angle(J2, J3);
			f = (2.00*tan(3.14159*0.25 + friction_angle*0.5) / cos(friction_angle))*((I1*K3 / 3) + sqrt(J2)*(K1*cos(theta) - K2*sin(theta)*sin(friction_angle) / sqrt(3)));
		}

		if (this->Get_threshold(cont) == 0) 
		{ 
			this->Set_threshold(c_max, cont);
			this->SetValue(INITIAL_THRESHOLD, c_max);
		}   // 1st iteration sets threshold as c_max

		c_threshold = this->Get_threshold(cont);
		//this->Set_NonConvergedf_sigma(f, cont);

		F = f - c_threshold;

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamages(cont);
		}
		else
		{
			damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
			if (damage > 0.99) { damage = 0.99; }
		}

		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - damage);
	}


	// ******* DAMAGE MECHANICS YIELD SURFACES AND EXPONENTIAL SOFTENING ********
	void FemDem3DElement::IntegrateStressDamageMechanics(Vector& rIntegratedStress, double& damage,
		const Vector StrainVector, const Vector StressVector, int cont, double l_char)
	{
		std::string yield_surface = this->GetProperties()[YIELD_SURFACE];

		if (yield_surface      == "ModifiedMohrCoulomb") { this->ModifiedMohrCoulombCriterion(rIntegratedStress, damage, StressVector, cont, l_char); }
		else if (yield_surface == "SimoJu") { this->SimoJuCriterion(rIntegratedStress, damage, StrainVector, StressVector, cont, l_char); }
		else if (yield_surface == "Rankine") { this->RankineCriterion(rIntegratedStress, damage, StressVector, cont, l_char); }
		else if (yield_surface == "DruckerPrager") { this->DruckerPragerCriterion(rIntegratedStress, damage, StressVector, cont, l_char); }
		else if (yield_surface == "RankineFragile") { this->RankineFragileLaw(rIntegratedStress, damage, StressVector, cont, l_char); }
		else { KRATOS_ERROR << " Yield Surface not defined "; }
	}

	void FemDem3DElement::ModifiedMohrCoulombCriterion(Vector& rIntegratedStress, double& damage, const Vector& StressVector,int cont, double l_char)
	{
		rIntegratedStress.resize(3);
		Vector PrincipalStressVector = ZeroVector(2);
		this->CalculatePrincipalStress(PrincipalStressVector, StressVector);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0;
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		friction_angle = this->GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159 / 180; // In radians!
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];

		// Check input variables 
		if (friction_angle < 1e-24) { friction_angle = 32 * 3.14159 / 180; std::cout << "Friction Angle not defined, assumed equal to 32� " << std::endl; }
		if (sigma_c < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa "; }
		if (sigma_t < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa "; }
		if (Gt < 1e-24) { KRATOS_ERROR << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa "; }

		double K1, K2, K3, Rmorh, R, alpha_r, c_max, theta, c_threshold;
		R = abs(sigma_c / sigma_t);
		Rmorh = pow(tan((3.14159 / 4) + friction_angle / 2), 2);
		alpha_r = R / Rmorh;
		c_max = abs(sigma_c);

		double I1, J2, J3;
		I1 = Calculate_I1_Invariant(PrincipalStressVector[0], PrincipalStressVector[1]);
		J2 = Calculate_J2_Invariant(PrincipalStressVector[0], PrincipalStressVector[1]);
		J3 = Calculate_J3_Invariant(PrincipalStressVector[0], PrincipalStressVector[1], I1);
		K1 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r)*sin(friction_angle);
		K2 = 0.5*(1 + alpha_r) - 0.5*(1 - alpha_r) / sin(friction_angle);
		K3 = 0.5*(1 + alpha_r)*sin(friction_angle) - 0.5*(1 - alpha_r);

		double n = sigma_c / sigma_t;
		double ElementArea = this->GetGeometry().Area();
		//double l_char = sqrt(4 * ElementArea / sqrt(3));

		double A = 1.00 / (n*n*Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		double f = 0.0, F = 0.0; /// F = f-c = 0 classical definition of yield surface

		// Check Modified Mohr-Coulomb criterion
		if (PrincipalStressVector[0] == 0 && PrincipalStressVector[1] == 0) { f = 0; }
		else
		{
			theta = Calculate_Theta_Angle(J2, J3);
			f = (2.00*tan(3.14159*0.25 + friction_angle*0.5) / cos(friction_angle))*((I1*K3 / 3) + sqrt(J2)*(K1*cos(theta) - K2*sin(theta)*sin(friction_angle) / sqrt(3)));
		}

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		//KRATOS_WATCH(c_threshold)
		this->Set_NonConvergedf_sigma(f, cont);
	
		F = f - c_threshold;

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamages(cont);
		}
		else
		{
			damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
			if (damage > 0.99) { damage = 0.99; }
		}

		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - damage);
	}

	void FemDem3DElement::RankineCriterion(Vector& rIntegratedStress, double& damage, const Vector& StressVector, int cont, double l_char) 
	{
		Vector PrincipalStressVector = ZeroVector(3);
		this->CalculatePrincipalStress(PrincipalStressVector, StressVector);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0, c_max = 0.0, c_threshold = 0.0;
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		friction_angle = this->GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159 / 180; // In radians!
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];
		c_max = abs(sigma_t);

		double ElementArea = this->GetGeometry().Area();
		//double l_char = sqrt(4 * ElementArea / sqrt(3));
		double A = 1.00 / (Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		double f, F; /// F = f-c = 0 classical definition of yield surface
		f = GetMaxValue(PrincipalStressVector);

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);

		F = f - c_threshold;

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamage();
			//this->Set_NonConvergeddamage(damage);
		}
		else
		{
			damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));  // Exponential softening law
			if (damage > 0.99) { damage = 0.99; }
			//this->Set_NonConvergeddamage(damage);
		}
		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - damage);
	}

	void FemDem3DElement::DruckerPragerCriterion(Vector& rIntegratedStress, double& damage, const Vector& StressVector, int cont, double l_char)
	{
		Vector PrincipalStressVector = ZeroVector(3);
		this->CalculatePrincipalStress(PrincipalStressVector, StressVector);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0;
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		friction_angle = this->GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159 / 180; // In radians!
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];

		// Check input variables 
		if (friction_angle < 1e-24) { friction_angle = 32 * 3.14159 / 180; std::cout << "Friction Angle not defined, assumed equal to 32� " << std::endl; }
		if (sigma_c < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in compression not defined, include YIELD_STRESS_C in .mdpa "; }
		if (sigma_t < 1e-24) { KRATOS_ERROR << " ERROR: Yield stress in tension not defined, include YIELD_STRESS_T in .mdpa "; }
		if (Gt < 1e-24) { KRATOS_ERROR << " ERROR: Fracture Energy not defined in the model part, include FRAC_ENERGY_T in .mdpa "; }

		double  c_max, c_threshold;
		c_max = abs(sigma_t*(3 + sin(friction_angle)) / (3 * sin(friction_angle) - 3));

		double I1, J2;
		I1 = Calculate_I1_Invariant(PrincipalStressVector[0], PrincipalStressVector[1]);
		J2 = Calculate_J2_Invariant(PrincipalStressVector[0], PrincipalStressVector[1]);

		double ElementArea = this->GetGeometry().Area();
		//double l_char = sqrt(4 * ElementArea / sqrt(3));
		double A = 1.00 / (Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		double f, F; /// F = f-c = 0 classical definition of yield surface

		double CFL = 0.0, TEN0 = 0.0;

		// Check DruckerPrager criterion
		if (PrincipalStressVector[0] == 0 && PrincipalStressVector[1] == 0) { f = 0; }
		else
		{
			CFL = -sqrt(3)*(3 - sin(friction_angle)) / (3 * sin(friction_angle) - 3);
			TEN0 = 6 * I1*sin(friction_angle) / (sqrt(3)*(3 - sin(friction_angle))) + sqrt(J2);
			f = abs(CFL*TEN0);
		}

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);
		F = f - c_threshold;

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamage();
			//this->Set_NonConvergeddamage(damage);
		}
		else
		{
			damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
			if (damage > 0.99) { damage = 0.99; }
			//this->Set_NonConvergeddamage(damage);
		}
		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - damage);
	}

	void FemDem3DElement::SimoJuCriterion(Vector& rIntegratedStress, double& damage,  const Vector& StrainVector,  const Vector& StressVector, int cont, double l_char)
	{
		Vector PrincipalStressVector = ZeroVector(3);
		this->CalculatePrincipalStress(PrincipalStressVector, StressVector);

		double sigma_t = 0.0, E = 0.0, Gt = 0.0, sigma_c = 0.0, n = 0;
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];

		double c_max, c_threshold;
		n = abs(sigma_c / sigma_t);
		c_max = abs(sigma_c) / sqrt(E);

		double SumA = 0.0, SumB = 0.0, SumC = 0.0, ere0 = 0.0, ere1 = 0.0;
		for (int cont = 0;cont < 1;cont++)
		{
			SumA += abs(PrincipalStressVector[cont]);
			SumB += 0.5*(PrincipalStressVector[cont] + abs(PrincipalStressVector[cont]));
			SumC += 0.5*(-PrincipalStressVector[cont] + abs(PrincipalStressVector[cont]));
		}
		ere0 = SumB / SumA;
		ere1 = SumC / SumA;

		double f = 0, F = 0; /// F = f-c = 0 classical definition of yield surface

							 // Check SimoJu criterion
		if (StrainVector[0] == 0 && StrainVector[1] == 0) { f = 0; }
		else
		{
			double auxf = 0.0;
			for (int cont = 0;cont < 2;cont++)
			{
				auxf += StrainVector[cont] * StressVector[cont];  // E*S
			}
			f = sqrt(auxf);
			f *= (ere0*n + ere1);
		}

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);
		F = f - c_threshold;

		double ElementArea = this->GetGeometry().Area();
		//double l_char = sqrt(4 * ElementArea / sqrt(3));
		double A = 1.00 / (Gt*n*n*E / (l_char *pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamage();
			//this->Set_NonConvergeddamage(damage);
		}
		else
		{
			damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
			if (damage > 0.99) { damage = 0.99; }
		//	this->Set_NonConvergeddamage(damage);
		}
		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - damage);
	}

	void FemDem3DElement::RankineFragileLaw(Vector& rIntegratedStress, double& damage, const Vector& StressVector, int cont, double l_char) 
	{
		Vector PrincipalStressVector = ZeroVector(3);
		this->CalculatePrincipalStress(PrincipalStressVector, StressVector);

		double sigma_c = 0.0, sigma_t = 0.0, friction_angle = 0.0, E = 0.0, Gt = 0.0, c_max = 0.0, c_threshold = 0.0;
		sigma_c = this->GetProperties()[YIELD_STRESS_C];
		sigma_t = this->GetProperties()[YIELD_STRESS_T];
		friction_angle = this->GetProperties()[INTERNAL_FRICTION_ANGLE] * 3.14159 / 180; // In radians!
		E = this->GetProperties()[YOUNG_MODULUS];
		Gt = this->GetProperties()[FRAC_ENERGY_T];
		c_max = abs(sigma_t);

		double ElementArea = this->GetGeometry().Area();
		//double l_char = sqrt(4 * ElementArea / sqrt(3));
		double A = 1.00 / (Gt*E / (l_char *pow(sigma_c, 2)) - 0.5);
		if (A < 0) { KRATOS_THROW_ERROR(std::invalid_argument, " 'A' damage parameter lower than zero --> Increase FRAC_ENERGY_T", A) }

		double f, F; /// F = f-c = 0 classical definition of yield surface
		f = GetMaxValue(PrincipalStressVector);

		if (this->Get_threshold(cont) == 0) { this->Set_threshold(c_max, cont); }   // 1st iteration sets threshold as c_max
		c_threshold = this->Get_threshold(cont);
		this->Set_NonConvergedf_sigma(f, cont);

		F = f - c_threshold;

		if (F <= 0)  // Elastic region --> Damage is constant 
		{
			damage = this->Get_Convergeddamage();
			//this->Set_NonConvergeddamage(damage);
		}
		else
		{
			//damage = 1 - (c_max / f)*exp(A*(1 - f / c_max));            // Exponential softening law
			//if (damage > 0.99) { damage = 0.99; }
			damage = 0.98;  // Fragile  law
		}
		rIntegratedStress = StressVector;
		rIntegratedStress *= (1 - damage);
	}

} // Element