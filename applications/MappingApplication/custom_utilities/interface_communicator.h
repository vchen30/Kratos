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

#if !defined(KRATOS_INTERFACE_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_INTERFACE_COMMUNICATOR_H_INCLUDED


// System includes
#include <vector>


// External includes


// Project includes
#include "includes/define.h"
#include "interface_search_structure.h"


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
  class InterfaceCommunicator
    {
    public:
      ///@name Type Definitions
      ///@{

      typedef InterfaceObjectConfigure::ContainerType TContainerType;

      /// Pointer definition of InterfaceCommunicator
      KRATOS_CLASS_POINTER_DEFINITION(InterfaceCommunicator);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      InterfaceCommunicator(ModelPart& rModelPartOrigin
                            Modelpart::Pointer pModelPartDestination);

      /// Destructor.
      virtual ~InterfaceCommunicator();


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

        void ExchangeInterfaceGeometryData()
        {
            TContainerType search_objects_origin;
            TContainerType search_objects_destination;

            CreateSearchObjectsOnOrigin(name..., search_objects_origin);
            CreateSearchObjectsOnDestination(name..., search_objects_destination);

            SelectCandidatePartitions(search_objects_destination);
            SendToOrigin();
            ReceiveFromDestination();
            SearchLocally();
            SendBacksToDestination();
            ReceiveFromOrigin();
            SelectBestSearchObject();
            AssingBestSearchObjectsToMapperConditions();
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
      virtual std::string Info() const;

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const;

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const;


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
        ModelPart& mrModelPartOrigin;
        ModelPart::Pointer mpModelPartDestination;


      ///@}
      ///@name Protected Operators
      ///@{


      ///@}
      ///@name Protected Operations
      ///@{

        void CreateSearchObjectsOnDestination(const std::string& rInterfaceSearchObjectName,
                                              TContainerType& rInterfaceSearchObjects)
        {
            // TODO select only the ones that have not found a neighbor...
            // Do this here? I gess so, since there will not be any search_structure_mpi any more
            const int num_mapper_conditions = static_cast<int>pModelPartDestination.NumberOfConditions(); // TODO local mesh

            conditions = ...

            const auto& rInterfaceSearchObject = MapperUtilities::GetInterfaceSearchObject(rInterfaceSearchObjectName);

            rInterfaceSearchObjects.resize(num_mapper_conditions);

            // TODO make parallel
            for (int i=0; i<num_mapper_conditions; ++i)
            {
                rInterfaceSearchObjects[i] = Kratos::make_shared<InterfaceSearchObject>(
                    rInterfaceSearchObject->Clone(conditions[i])); // TODO or Create ...
            }
        }

        void CreateSearchObjectsOnOrigin(const std::string& rInterfaceSearchObjectName,
                                         TContainerType& rInterfaceSearchObjects)
        {
            if (InterfaceSearchObjectType == MapperUtilities::Nodes)
            {
                const int num_nodes = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfNodes()

                nodes = ...

                rInterfaceSearchObjects.resize(num_nodes);

                // TODO make parallel
                for (int i=0; i<num_nodes; ++i)
                    rInterfaceSearchObjects[i] = Kratos::make_shared<InterfaceSearchNode>(nodes[i]);
            }
            else if (InterfaceSearchObjectType == MapperUtilities::GeometricalObjects)
            {
                // TODO use either elements or conditions
                const int num_geom_objects = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfConditions()

                conditions = ...

                rInterfaceSearchObjects.resize(num_geom_objects);

                // TODO make parallel
                for (int i=0; i<num_geom_objects; ++i)
                    rInterfaceSearchObjects[i] = Kratos::make_shared<InterfaceSearchGeometricalObject>(conditions[i]);
            }
            else
            {
                KRATOS_ERROR << "sth" << std::endl;
            }
        }

        virtual void SelectCandidatePartitions()
        {
            // exchange Bounding Boxes

            // can be parallelized if the insertion is atomic/in a critical section
            for (const auto& bounding_box : bounding_boxes)
            {
                for (auto& interface_object : rInterfaceObjects)
                {
                    if (bounding_box.IsInside(interface_object.Coordinates()))
                        interface_object.SetCandidatePartitionRank(bounding_box.GetRank())
                }
            }
        }

        virtual void AsyncSendObjectsToOrigin()
        {


        }

        virtual void SearchLocally()
        {

        }

        virtual void AsyncSendObjectsBackToDestination()
        {
            // send back only the objects that have a successful search

        }

        void SelectBestSearchObject()
        {

        }

        void AssingBestSearchObjectsToMapperConditions()
        {

        }

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


      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{


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
      InterfaceCommunicator& operator=(InterfaceCommunicator const& rOther);

      /// Copy constructor.
      InterfaceCommunicator(InterfaceCommunicator const& rOther);


      ///@}

    }; // Class InterfaceCommunicator

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    InterfaceCommunicator& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const InterfaceCommunicator& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
