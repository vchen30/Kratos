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

#if !defined(KRATOS_INTERFACE_PREPROCESS_H_INCLUDED )
#define  KRATOS_INTERFACE_PREPROCESS_H_INCLUDED



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
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

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
  
  /// Short class definition.
  /** Detail class definition.
  */
  class InterfacePreprocess
    {
    public:
        ///@name Type Definitions
        ///@{
        
        /// Pointer definition of InterfacePreprocess
        KRATOS_CLASS_POINTER_DEFINITION(InterfacePreprocess);
    
        ///@}
        ///@name Life Cycle 
        ///@{ 
        
        /// Default constructor.
        InterfacePreprocess(ModelPart& rModelPart) : 
            mrModelPart(rModelPart)
        {
            
        }

        /// Destructor.
        virtual ~InterfacePreprocess(){}
        

        ///@}
        ///@name Operators 
        ///@{
        
        
        ///@}
        ///@name Operations
        ///@{

        ModelPart::Pointer pGetInterfaceModelPart()
        {
            return mpInterfaceModelPart;
        }

        void GenerateInterfacePart(Parameters InterfaceParameters)
        {
            CheckAndValidateParameters(InterfaceParameters);      

            // Create a new / overwrite the InterfaceModelpart
            mpInterfaceModelPart = ModelPart::Pointer( new ModelPart("Mapper-Interface") );

            // mDimension = mrModelPart.GetProcessInfo[DOMAIN_SIZE]

            // create dummy properties // TODO is this ok? Or is there a case where I cannot do this?            
            Properties::Pointer dummy_properties = Properties::Pointer( new Properties() );
            mpInterfaceModelPart->AddProperties(dummy_properties);


            AddNodes();
            CreateConditions(InterfaceParameters);

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
            if (MapperUtilities::TotalProcesses() > 1)
            {
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
            buffer << "InterfacePreprocess" ;
            return buffer.str();
        }
        
        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "InterfacePreprocess";}

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
            "condition_name" : "",
            "use_nodes"      : false
        }  )" );
            
            
        ///@} 
        ///@name Private Operators
        ///@{ 
            
            
        ///@} 
        ///@name Private Operations
        ///@{ 
            
        void CheckAndValidateParameters(Parameters InterfaceParameters)
        {
            InterfaceParameters.RecursivelyValidateAndAssignDefaults(mDefaultParameters);

            KRATOS_ERROR_IF(InterfaceParameters["condition_name"].GetString() == "") 
                << "Condition name for Interface-ModelPart not specified" << std::endl;
        }

        void AddNodes()
        {
            // Store pointers to all interface nodes
            unsigned int nodes_counter = 0;
            for (ModelPart::NodesContainerType::const_iterator node_it = mrModelPart.NodesBegin(); node_it != mrModelPart.NodesEnd(); ++node_it)
            {
                mpInterfaceModelPart->Nodes().push_back( *(node_it.base()) ); // TODO resize? I know the size in advance...
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




            // std::string condition_name = InterfaceParameters["condition_name"].GetString();
            // condition_name.append("MapperCondition");

            std::string condition_name = "NearestNeighborMapperCondition3D1N";

            std::string final_string;

            unsigned int cond_counter = 0;

            if (InterfaceParameters["use_nodes"].GetBool()) // For node-based mappers, e.g. Nearest Neighbor Mapper
            {
                Condition::NodesArrayType temp_condition_nodes;
                final_string = "1N";
                for (ModelPart::NodesContainerType::const_iterator it_node = mrModelPart.GetCommunicator().LocalMesh().NodesBegin(); 
                    it_node != mrModelPart.GetCommunicator().LocalMesh().NodesEnd(); ++it_node)
                {
                    temp_condition_nodes.clear();
                    temp_condition_nodes.push_back( *(it_node).base()); // TODO make more efficient, w/o reallocation of memmory all the time!
                    CreateNewNodeCondition(temp_condition_nodes, cond_counter, condition_name);
                    ++cond_counter;
                }
            }
            else // using the "Geometry" of the Conditions or Elements that are used to construct the Interface
            {

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
        //   InterfacePreprocess& operator=(InterfacePreprocess const& rOther){}

        /// Copy constructor.
        // InterfacePreprocess(InterfacePreprocess const& rOther){}

            
        ///@}    
        
    }; // Class InterfacePreprocess 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream, 
// 				    InterfacePreprocess& rThis){}

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const InterfacePreprocess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_INTERFACE_PREPROCESS_H_INCLUDED  defined 


