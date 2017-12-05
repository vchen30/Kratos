//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher 
//
	        

// System includes


// External includes 


// Project includes
#include "includes/checks.h"
#include "shell_utilities.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{
    namespace ShellUtilities
	{
		typedef std::size_t 			SizeType;
        typedef Properties 				PropertiesType; // TODO remove this?

		double dN_seren_dxi(const int ThisNodeNumber, const double Xi, 
							const double Eta)
		{
			// Natural derivatives of 8-node serendipity shape functions

			double return_value;
			switch (ThisNodeNumber)
			{
			case 1:
				return_value = -(-Eta + 1.0)*(-0.25*Xi + 0.25) -
					0.25*(-Eta + 1.0)*(-Eta - Xi - 1.0);
				break;
			case 2:
				return_value = (-Eta + 1.0)*(0.25*Xi + 0.25) +
					0.25*(-Eta + 1.0)*(-Eta + Xi - 1.0);
				break;
			case 3:
				return_value = (Eta + 1.0)*(0.25*Xi + 0.25) +
					0.25*(Eta + 1.0)*(Eta + Xi - 1.0);
				break;
			case 4:
				return_value = -(Eta + 1.0)*(-0.25*Xi + 0.25) -
					0.25*(Eta + 1.0)*(Eta - Xi - 1.0);
				break;
			case 5:
				return_value = -1.0*Xi*(-Eta + 1.0);
				break;
			case 6:
				return_value = -0.5*Eta*Eta + 0.5;
				break;
			case 7:
				return_value = -1.0*Xi*(Eta + 1.0);
				break;
			case 8:
				return_value = 0.5*Eta*Eta - 0.5;
				break;
			default:
				KRATOS_ERROR <<
					"Error: ELEMENT ShellThinElement3D4N, METHOD dN_seren_dxi"
					<< std::endl;
			}

			return return_value;
		}

		double dN_seren_deta(const int ThisNodeNumber, const double Xi,
							 const double Eta)
		{
			// Natural derivatives of 8-node serendipity shape functions

			double return_value;
			switch (ThisNodeNumber)
			{
			case 1:
				return_value = -(-Eta + 1.0)*(-0.25*Xi + 0.25) -
					(-0.25*Xi + 0.25)*(-Eta - Xi - 1.0);
				break;
			case 2:
				return_value = -(-Eta + 1.0)*(0.25*Xi + 0.25) -
					(0.25*Xi + 0.25)*(-Eta + Xi - 1.0);
				break;
			case 3:
				return_value = (Eta + 1.0)*(0.25*Xi + 0.25) +
					(0.25*Xi + 0.25)*(Eta + Xi - 1.0);
				break;
			case 4:
				return_value = (Eta + 1.0)*(-0.25*Xi + 0.25) +
					(-0.25*Xi + 0.25)*(Eta - Xi - 1.0);
				break;
			case 5:
				return_value = 0.5*Xi*Xi - 0.5;
				break;
			case 6:
				return_value = -1.0*Eta*(Xi + 1.0);
				break;
			case 7:
				return_value = -0.5*Xi*Xi + 0.5;
				break;
			case 8:
				return_value = -1.0*Eta*(-Xi + 1.0);
				break;
			default:
				KRATOS_ERROR <<
					"Error: ELEMENT ShellThinElement3D4N, METHOD dN_seren_dxi"
					<< std::endl;
			}

			return return_value;
		}
		
		void InterpToStandardGaussPoints(double& rV1, double& rV2,
			double& rV3)
		{
			double vg1 = rV1;
			double vg2 = rV2;
			double vg3 = rV3;
#ifdef OPT_AVERAGE_RESULTS
			rV1 = (vg1 + vg2 + vg3) / 3.0;
			rV2 = (vg1 + vg2 + vg3) / 3.0;
			rV3 = (vg1 + vg2 + vg3) / 3.0;
#else
			rV1 = (2.0*vg1) / 3.0 - vg2 / 3.0 + (2.0*vg3) / 3.0;
			rV2 = (2.0*vg1) / 3.0 + (2.0*vg2) / 3.0 - vg3 / 3.0;
			rV3 = (2.0*vg2) / 3.0 - vg1 / 3.0 + (2.0*vg3) / 3.0;
#endif // OPT_AVERAGE_RESULTS
		}

		void InterpToStandardGaussPoints(std::vector< double >& rVec)
		{
			if (rVec.size() != 3) return;
			InterpToStandardGaussPoints(rVec[0], rVec[1], rVec[2]);
		}

		void InterpToStandardGaussPoints(std::vector< array_1d<double,3> >& rVec)
		{
			if (rVec.size() != 3) return;
			for (SizeType i = 0; i < 3; i++)
				InterpToStandardGaussPoints(rVec[0][i], rVec[1][i], rVec[2][i]);
		}

		void InterpToStandardGaussPoints(std::vector< array_1d<double,6> >& rVec)
		{
			if (rVec.size() != 3) return;
			for (SizeType i = 0; i < 6; i++)
				InterpToStandardGaussPoints(rVec[0][i], rVec[1][i], rVec[2][i]);
		}

		void InterpToStandardGaussPoints(std::vector< Vector >& rVec)
		{
			if (rVec.size() != 3) return;
			SizeType ncomp = rVec[0].size();
			for (int i = 1; i < 3; i++)
				if (rVec[i].size() != ncomp)
					return;
			for (SizeType i = 0; i < ncomp; i++)
				InterpToStandardGaussPoints(rVec[0][i], rVec[1][i], rVec[2][i]);
		}

		void InterpToStandardGaussPoints(std::vector< Matrix >& rVec)
		{
			if (rVec.size() != 3) return;
			SizeType nrows = rVec[0].size1();
			SizeType ncols = rVec[0].size2();
			for (int i = 1; i < 3; i++)
				if (rVec[i].size1() != nrows || rVec[i].size2() != ncols)
					return;
			for (SizeType i = 0; i < nrows; i++)
				for (SizeType j = 0; j < ncols; j++)
					InterpToStandardGaussPoints
					(rVec[0](i, j), rVec[1](i, j), rVec[2](i, j));
		}



        void CheckVariables()
        {
            KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);
            KRATOS_CHECK_VARIABLE_KEY(ROTATION);
            KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
            KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
            KRATOS_CHECK_VARIABLE_KEY(DENSITY);
            KRATOS_CHECK_VARIABLE_KEY(SHELL_CROSS_SECTION);
            KRATOS_CHECK_VARIABLE_KEY(THICKNESS);
            KRATOS_CHECK_VARIABLE_KEY(CONSTITUTIVE_LAW);        
        }

        void CheckDofs(GeometryType& rGeom)
        {
            // verify that the dofs exist
            for (unsigned int i = 0; i < rGeom.size(); i++)
            {
                auto& r_node = rGeom[i];
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT, r_node);
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ROTATION, r_node);

                KRATOS_CHECK_DOF_IN_NODE(ROTATION_X, r_node);
                KRATOS_CHECK_DOF_IN_NODE(ROTATION_Y, r_node);
                KRATOS_CHECK_DOF_IN_NODE(ROTATION_Z, r_node);

                KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, r_node);
                KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, r_node);
                KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, r_node);

                if (r_node.GetBufferSize() < 2)
                    KRATOS_ERROR << "This Element needs at least a buffer size = 2" << std::endl;
            }
        }

        void CheckProperties(const Element* pTheElement, const ProcessInfo& rCurrentProcessInfo, const bool IsThickShell)
        {
            // check properties
			KRATOS_ERROR_IF_NOT(pTheElement->pGetProperties()) << "Properties not provided for element " 
															   << pTheElement->Id() << std::endl;

            const PropertiesType & props = pTheElement->GetProperties();

            const GeometryType& geom = pTheElement->GetGeometry(); // TODO check if this can be const

            if(props.Has(SHELL_CROSS_SECTION)) // if the user specified a cross section ...
            {
                const ShellCrossSection::Pointer & section = props[SHELL_CROSS_SECTION];
				KRATOS_ERROR_IF_NOT(section) << "SHELL_CROSS_SECTION not provided for element " << pTheElement->Id() << std::endl;
        
                section->Check(props, geom, rCurrentProcessInfo);
            }
            else if (props.Has(SHELL_ORTHOTROPIC_LAYERS))
            {
                CheckSpecificProperties(pTheElement, props, IsThickShell);
        
                // perform detailed orthotropic check later in shell_cross_section
            }
            else // ... allow the automatic creation of a homogeneous section from a material and a thickness
            {
                CheckSpecificProperties(pTheElement, props, IsThickShell);
        
                ShellCrossSection::Pointer dummySection = ShellCrossSection::Pointer(new ShellCrossSection());
                dummySection->BeginStack();
                dummySection->AddPly(props[THICKNESS], 0.0, 5, pTheElement->pGetProperties());
                dummySection->EndStack();
                dummySection->SetSectionBehavior(ShellCrossSection::Thick);
                dummySection->Check(props, geom, rCurrentProcessInfo);
            }

        }

        void CheckSpecificProperties(const Element* pTheElement, const PropertiesType & rProps, const bool IsThickShell)
        {
			KRATOS_ERROR_IF_NOT(rProps.Has(CONSTITUTIVE_LAW)) << "CONSTITUTIVE_LAW not provided for element " 
															  << pTheElement->Id() << std::endl;

			KRATOS_ERROR_IF_NOT(rProps[CONSTITUTIVE_LAW]) << "CONSTITUTIVE_LAW not provided for element " 
									  					  << pTheElement->Id() << std::endl;

            KRATOS_ERROR_IF_NOT(rProps.Has(THICKNESS)) << "THICKNESS not provided for element " 
													   << pTheElement->Id() << std::endl;
            KRATOS_ERROR_IF(rProps[THICKNESS] <= 0.0) << "wrong THICKNESS value provided for element " 
													  << pTheElement->Id() << std::endl;

            KRATOS_ERROR_IF_NOT(rProps.Has(DENSITY)) << "DENSITY not provided for element " 
													 << pTheElement->Id() << std::endl;
            KRATOS_ERROR_IF(rProps[DENSITY] < 0.0) << "wrong DENSITY value provided for element " 
												   << pTheElement->Id() << std::endl;

            if(IsThickShell)
            {
                // Check constitutive law has been verified with Stenberg stabilization
                // applicable for 5-parameter shells only.
				bool stenberg_stabilization_suitable = false;
				claw->GetValue(STENBERG_SHEAR_STABILIZATION_SUITABLE, stenberg_stabilization_suitable);
                if (!stenberg_stabilization_suitable)
                {
                    std::cout << "\nWARNING:\nThe current constitutive law has not been checked with Stenberg shear stabilization."
                        << "\nPlease check results carefully."
                        << std::endl;
                }
            }
        }  
    } // namespace ShellUtilities
}  // namespace Kratos.


