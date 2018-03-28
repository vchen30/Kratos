//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_INTERFACE_PREPROCESSOR_H_INCLUDED )
#define  KRATOS_INTERFACE_PREPROCESSOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_2d_4.h"

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_utilities/parallel_fill_communicator.h"
#endif


namespace Kratos
{
  ///@addtogroup ApplicationNameApplication
  ///@{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}
  ///@name Kratos Classes
  ///@{

  /// This class creates the mapper conditions
  /** Detail class definition.
  */
  class InterfacePreprocessor
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of InterfacePreprocessor
        KRATOS_CLASS_POINTER_DEFINITION(InterfacePreprocessor);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        InterfacePreprocessor(ModelPart& rModelPart, const ModelPart::Pointer& pModelPart) :
            mrModelPart(rModelPart), mpInterfaceModelPart(pModelPart)
        {

        }

        /// Destructor.
        virtual ~InterfacePreprocessor(){}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        void GenerateInterfacePart(Parameters InterfaceParameters)
        {
            CheckAndValidateParameters(InterfaceParameters);

            mpInterfaceModelPart->GetMesh().Clear(); // Clear the ModelPart


            // mDimension = mrModelPart.GetProcessInfo[DOMAIN_SIZE]

            // create dummy properties // TODO is this ok? Or is there a case where I cannot do this?
            Properties::Pointer dummy_properties = Properties::Pointer( new Properties() );
            mpInterfaceModelPart->AddProperties(dummy_properties);


            AddNodes();
            CreateConditions(InterfaceParameters);

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
            if (MapperUtilities::TotalProcesses() > 1)
            {
                // TODO check if this is actually necessary
                // Set the MPICommunicator
                std::cout << "Doing the ParallelFillCommunicator stuff" << std::endl;
                ParallelFillCommunicator parallel_fill_communicator(*mpInterfaceModelPart);
                parallel_fill_communicator.Execute();
            }
#endif

        }


        ///@}
        ///@name Access
        ///@{


        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        virtual std::string Info() const
        {
            std::stringstream buffer;
            buffer << "InterfacePreprocessor" ;
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "InterfacePreprocessor";}

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const {}


        ///@}
        ///@name Friends
        ///@{


        ///@}

    protected:
        ///@name Protected static Member Variables
        ///@{


        ///@}
        ///@name Protected member Variables
        ///@{


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


        ///@}
        ///@name Protected  Access
        ///@{


        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}

    private:
        ///@name Static Member Variables
        ///@{


        ///@}
        ///@name Member Variables
        ///@{
        ModelPart& mrModelPart;
        ModelPart::Pointer mpInterfaceModelPart;
        int mDimension;

        Parameters mDefaultParameters = Parameters( R"(
        {
            "mapper_condition_name" : "",
            "use_nodes"      : true
        }  )" );
        // Use "use_nodes" if the geometry on the destination is not needed


        ///@}
        ///@name Private Operators
        ///@{


        ///@}
        ///@name Private Operations
        ///@{

        void CheckAndValidateParameters(Parameters InterfaceParameters)
        {
            InterfaceParameters.RecursivelyValidateAndAssignDefaults(mDefaultParameters);

            KRATOS_ERROR_IF(InterfaceParameters["mapper_condition_name"].GetString() == "")
                << "Condition name for Interface-ModelPart not specified" << std::endl;
        }

        void AddNodes()
        {
            // Store pointers to all interface nodes
            unsigned int nodes_counter = 0;

            // TODO ptr iterator? => see builderandsolver
            const int num_nodes_local = mrModelPart.GetCommunicator().LocalMesh().NumberOfNodes();

            mpInterfaceModelPart->Nodes().reserve(num_nodes_local);

            for (ModelPart::NodesContainerType::const_iterator node_it = mrModelPart.NodesBegin(); node_it != mrModelPart.NodesEnd(); ++node_it)
            {
                mpInterfaceModelPart->Nodes().push_back( *(node_it.base()) );
                ++nodes_counter;
            }
        }

        void CreateConditions(Parameters InterfaceParameters)
        {
            //     // Generate Conditions from original the edges that can be considered interface
            // for (ModelPart::ElementsContainerType::const_iterator it_elem = rOriginPart.ElementsBegin(); it_elem != rOriginPart.ElementsEnd(); it_elem++)
            // {
            //     for (unsigned int i_edge = 0; i_edge < (*it_elem).GetGeometry().EdgesNumber(); i_edge++)
            //     {
            //         unsigned int count = 0;
            //         const unsigned int number_of_points = (*it_elem).GetGeometry().Edges()[i_edge].PointsNumber();
            //         for (unsigned int i_node = 0; i_node < number_of_points; i_node++)
            //         {
            //             if ((*it_elem).GetGeometry().Edges()[i_edge][i_node].IsDefined(INTERFACE) == true)
            //             {
            //                 if ((*it_elem).GetGeometry().Edges()[i_edge][i_node].Is(INTERFACE) == true)
            //                 {
            //                     count++;
            //                 }
            //             }
            //         }

            //         if (count == number_of_points)
            //         {
            //             std::string Edgecondition_name = condition_name;
            //             cond_id += 1; // NOTE: To paralellize be careful with this ID
            //             if (number_of_points == 2)
            //             {
            //                 Edgecondition_name.append("Condition2D2N");
            //                 Edgecondition_name.append(final_string);

            //                 CreateNewCondition(rInterfacePart, *(it_elem.base()), (*it_elem).GetGeometry().Edges()[i_edge], cond_id, Edgecondition_name);
            //                 cond_counter ++;
            //             }
            //             else
            //             {
            //                 if (simplest_geometry == false)
            //                 {
            //                     Edgecondition_name.append("Condition2D3N");
            //                     Edgecondition_name.append(final_string);

            //                     CreateNewCondition(rInterfacePart, *(it_elem.base()), (*it_elem).GetGeometry().Edges()[i_edge], cond_id, Edgecondition_name);
            //                     cond_counter ++;
            //                 }
            //                 else
            //                 {
            //                     Edgecondition_name.append("Condition2D2N");
            //                     Edgecondition_name.append(final_string);

            //                     Line2D2< Node<3> > lin_1((*it_elem).GetGeometry().Edges()[i_edge](0), (*it_elem).GetGeometry().Edges()[i_edge](1));
            //                     CreateNewCondition(rInterfacePart, *(it_elem.base()), lin_1, cond_id, Edgecondition_name);
            //                     cond_counter ++;
            //                     cond_id += 1;
            //                     Line2D2< Node<3> > lin_2((*it_elem).GetGeometry().Edges()[i_edge](1), (*it_elem).GetGeometry().Edges()[i_edge](2));
            //                     CreateNewCondition(rInterfacePart, *(it_elem.base()), lin_2, cond_id, Edgecondition_name);
            //                     cond_counter ++;
            //                 }
            //             }
            //         }
            //     }


            // TODO Reserve the size of the conditions, since I know it beforehand


            std::string condition_name = InterfaceParameters["mapper_condition_name"].GetString();
            condition_name.append("3D"); // TODO how to select 2D or 3D? Jordi

            unsigned int cond_counter = 0;

            if (InterfaceParameters["use_nodes"].GetBool()) // For mappers which don't depend on the destination geometry, e.g. Nearest Neighbor Mapper or Nearest Element
            {
                Condition::NodesArrayType temp_condition_nodes;
                condition_name.append("1N");
                for (ModelPart::NodesContainerType::const_iterator it_node = mrModelPart.GetCommunicator().LocalMesh().NodesBegin();
                    it_node != mrModelPart.GetCommunicator().LocalMesh().NodesEnd(); ++it_node)
                {
                    temp_condition_nodes.clear();
                    temp_condition_nodes.push_back( *(it_node).base()); // TODO make more efficient, w/o reallocation of memory all the time!
                    CreateNewNodeCondition(temp_condition_nodes, cond_counter, condition_name);
                    ++cond_counter;
                }
            }
            else // using the "Geometry" of the Conditions or Elements that are used to construct the Interface
            {
                // Use pGetGeometry to not copy the geometry!

            }
        }

        void CreateNewNodeCondition(
            Condition::NodesArrayType& rNode, // This is always only one node!
            const unsigned int CondId,
            const std::string ConditionName
            )
        {
            KRATOS_TRY;

            Condition const & rCondition = KratosComponents<Condition>::Get(ConditionName);
            Condition::Pointer p_cond = Condition::Pointer(rCondition.Create(CondId, rNode, mpInterfaceModelPart->pGetProperties(0)));
            mpInterfaceModelPart->AddCondition(p_cond);

            KRATOS_CATCH("");
        }

        void CreateNewGeometryCondition(
            Condition::NodesArrayType& rNode, // This is always only one node!
            const unsigned int CondId,
            const std::string ConditionName
            )
        {
            KRATOS_TRY;

            Condition const & rCondition = KratosComponents<Condition>::Get(ConditionName);
            Condition::Pointer p_cond = Condition::Pointer(rCondition.Create(CondId, rNode, mpInterfaceModelPart->pGetProperties(0)));
            mpInterfaceModelPart->AddCondition(p_cond);

            KRATOS_CATCH("");
        }




        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        //   InterfacePreprocessor& operator=(InterfacePreprocessor const& rOther){}

        /// Copy constructor.
        // InterfacePreprocessor(InterfacePreprocessor const& rOther){}


        ///@}

    }; // Class InterfacePreprocessor

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream,
// 				    InterfacePreprocessor& rThis){}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const InterfacePreprocessor& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_PREPROCESSOR_H_INCLUDED  defined


