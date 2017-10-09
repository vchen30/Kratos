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

#if !defined(KRATOS_MAPPER_H_INCLUDED )
#define  KRATOS_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/mapper_communicator.h"
#include "custom_utilities/mapper_utilities.h"
#include "custom_utilities/mapper_flags.h"

// For MPI-parallel Mapper
#ifdef KRATOS_USING_MPI
#include "mpi.h"
#include "custom_utilities/mapper_mpi_communicator.h"
#endif


namespace Kratos
{

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

/// Base Class for all Mappers
/** This is the base class for every mapper.
* It contains the three pure virtual functions that have to be implemented by every mapper:
* - Map: Basic function that maps a field from one ModelPart to another Modelpart
*        Mapping Direction: Origin => Destionation
* - InverseMap: This function does the opposite of the "Map" function
*               Mapping Direction: Destination => Origin
* - UpdateInterface: Called when the interface is changed. It recomputes the neighbors and 
*   other information related to the relations btw entities (node, elements,...) on the interfaces
* It is also responsible for initializing the MapperCommunicator or the MapperMPICommuniator
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/

class Mapper
{
public:

    ///@name Type Definitions
    ///@{
    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of Mapper
    KRATOS_CLASS_POINTER_DEFINITION(Mapper);

    ///@}
    ///@name Life Cycle
    ///@{


    /// Destructor.
    virtual ~Mapper()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius)
    {
        mpMapperCommunicator->UpdateInterface(MappingOptions, SearchRadius);
        
        if (mpInverseMapper) mpInverseMapper->UpdateInterface(MappingOptions, SearchRadius);

        if (MappingOptions.Is(MapperFlags::REMESHED)) ComputeNumberOfNodesAndConditions();

        UpdateInterfaceSpecific(MappingOptions);
    }

    /**
    This function maps from Origin to Destination, SCALAR version
    */
    virtual void Map(const Variable<double>& rOriginVariable,
                     const Variable<double>& rDestinationVariable,
                     Kratos::Flags MappingOptions) = 0;

    /**
    This function maps from Origin to Destination, VECTOR version
    */
    virtual void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
                     const Variable< array_1d<double, 3> >& rDestinationVariable,
                     Kratos::Flags MappingOptions) = 0;

    /**
    This function maps from Destination to Origin, SCALAR version
    Implementation can be overridden in derived classes
    */
    virtual void InverseMap(const Variable<double>& rOriginVariable,
                            const Variable<double>& rDestinationVariable,
                            Kratos::Flags MappingOptions)
    {
        // Construct the inverse mapper if it hasn't been done before
        if (!mpInverseMapper) InitializeInverseMapper();

        mpInverseMapper->Map(rDestinationVariable, rOriginVariable, MappingOptions);
    }

    /**
    This function maps from Destination to Origin, VECTOR version
    Implementation can be overridden in derived classes
    */
    virtual void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                            const Variable< array_1d<double, 3> >& rDestinationVariable,
                            Kratos::Flags MappingOptions)
    {
        // Construct the inverse mapper if it hasn't been done before
        if (!mpInverseMapper) InitializeInverseMapper();

        mpInverseMapper->Map(rDestinationVariable, rOriginVariable, MappingOptions);
    }

    MapperCommunicator::Pointer pGetMapperCommunicator()
    {
        return mpMapperCommunicator;
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
        return "Mapper";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Mapper";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }

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

    Mapper::Pointer mpInverseMapper;

    ModelPart& mrModelPartOrigin; // TODO const ref?
    ModelPart& mrModelPartDestination;

    Parameters mJsonParameters;

    MapperCommunicator::Pointer mpMapperCommunicator;

    // global, aka of the entire submodel-parts
    int mNumConditionsOrigin;
    int mNumConditionsDestination;

    int mNumNodesOrigin;
    int mNumNodesDestination;

    int mEchoLevel = 0;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{


    /**
    This function can be overridden by derived classes to do additional things
    */
    virtual void UpdateInterfaceSpecific(Kratos::Flags MappingOptions) {}

    /**
    This function must be implemented by every mapper to initialize the inverse mapper
    The InverseMapper is constructed with the order of the model_parts switched!
    */
    virtual void InitializeInverseMapper() = 0;

    // Constructor, can only be called by derived classes (actual mappers)
    Mapper(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination,
           Parameters JsonParameters) :
        mrModelPartOrigin(rModelPartOrigin),
        mrModelPartDestination(rModelPartDestination),
        mJsonParameters(JsonParameters)
    {
        mpInverseMapper.reset(); // explicitly specified to be safe
        ComputeNumberOfNodesAndConditions();

        // Create the mapper communicator
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        if (MapperUtilities::TotalProcesses() > 1)
        {
            mpMapperCommunicator = MapperCommunicator::Pointer (
                                        new MapperMPICommunicator(mrModelPartOrigin,
                                                mrModelPartDestination,
                                                mJsonParameters) );
        }
        else
        {
            InitializeSerialCommunicator();
        }

#else // serial compilation
        InitializeSerialCommunicator();
#endif
        // Access the Parameters only after the communicator is constructed,
        // bcs they are checked and validated there!
        mEchoLevel = JsonParameters["echo_level"].GetInt();
    }

    void ComputeNumberOfNodesAndConditions()
    {
        // Compute the quantities of the local model_parts
        mNumConditionsOrigin = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfConditions();
        mNumConditionsDestination = mrModelPartDestination.GetCommunicator().LocalMesh().NumberOfConditions();

        mNumNodesOrigin = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfNodes();
        mNumNodesDestination = mrModelPartDestination.GetCommunicator().LocalMesh().NumberOfNodes();

        // Compute the quantities of the global model_parts
        mrModelPartOrigin.GetCommunicator().SumAll(mNumConditionsOrigin);
        mrModelPartDestination.GetCommunicator().SumAll(mNumConditionsDestination);

        mrModelPartOrigin.GetCommunicator().SumAll(mNumNodesOrigin);
        mrModelPartDestination.GetCommunicator().SumAll(mNumNodesDestination);
    }

    void ProcessMappingOptions(const Kratos::Flags& rMappingOptions,
                               double& Factor)
    {
        if (rMappingOptions.Is(MapperFlags::SWAP_SIGN))
        {
            Factor *= (-1);
        }
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

    void InitializeSerialCommunicator()
    {
        mpMapperCommunicator = MapperCommunicator::Pointer (
                                   new MapperCommunicator(mrModelPartOrigin,
                                           mrModelPartDestination,
                                           mJsonParameters) );
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
    Mapper& operator=(Mapper const& rOther);

    /// Copy constructor.
    //Mapper(Mapper const& rOther);

    ///@}

}; // Class Mapper

}  // namespace Kratos.

#endif // KRATOS_MAPPER_H_INCLUDED  defined