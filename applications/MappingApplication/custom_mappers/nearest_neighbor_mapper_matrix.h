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

#if !defined(KRATOS_NEAREST_NEIGHBOR_MAPPER_MATRIX_H_INCLUDED )
#define  KRATOS_NEAREST_NEIGHBOR_MAPPER_MATRIX_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper_matrix_based.h"


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
/// Short class definition.
template <class TMappingMatrixBuilder,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
>
class NearestNeighborMapperMatrix : public MapperMatrixBased<TMappingMatrixBuilder, TLinearSolver>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NearestNeighborMapperMatrix
    KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborMapperMatrix);
    
    typedef TMappingMatrixBuilder TMappingMatrixBuilderType;

    typedef typename TMappingMatrixBuilderType::TDataType TDataType;
    
    typedef typename TMappingMatrixBuilderType::TSystemMatrixType TSystemMatrixType;

    typedef typename TMappingMatrixBuilderType::TSystemVectorType TSystemVectorType;

    typedef typename TMappingMatrixBuilderType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename TMappingMatrixBuilderType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename TMappingMatrixBuilderType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename TMappingMatrixBuilderType::LocalSystemVectorType LocalSystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    NearestNeighborMapperMatrix(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                                Parameters rJsonParameters) : MapperMatrixBased<TMappingMatrixBuilder, TLinearSolver>(
                                i_model_part_origin, i_model_part_destination, rJsonParameters)
    {
        KRATOS_WATCH("NearestNeighborMapperMatrix Constructor")

        this->mInterfaceParameters = Parameters( R"(
        {
            "condition_name" : "nearest_neighbor",
            "use_nodes"      : true
        }  )" );
        this->GenerateInterfaceModelPart();

        KRATOS_WATCH("Before Communicator Access")
        KRATOS_WATCH(this->mpInterfaceModelPart)
        std::cout << this->mpInterfaceModelPart << std::endl;
        // std::cout << "Rank " << this->mpInterfaceModelPart->GetCommunicator().MyPID() << ", num_local_nodes = " 
        //           << this->mpInterfaceModelPart->GetCommunicator().LocalMesh().NumberOfNodes() << std::endl;

        KRATOS_WATCH("After Communicator Access")

        this->mpMapperCommunicator->InitializeOrigin(MapperUtilities::Node_Coords);
        this->mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);
        this->mpMapperCommunicator->Initialize();

        this->ComputeMappingMatrix();
        
        
        
        // this->mpMappingMatrixBuilder->Update(this->mrModelPartOrigin, this->mpQd, PRESSURE);
    }

    /// Destructor.
    virtual ~NearestNeighborMapperMatrix() { }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable<double>& rOriginVariable,
                    const Variable<double>& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {

    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                    const Variable< array_1d<double, 3> >& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {
        
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
    virtual std::string Info() const override
    {
        return "NearestNeighborMapperMatrix";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestNeighborMapperMatrix";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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
    // NearestNeighborMapperMatrix<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(NearestNeighborMapperMatrix<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    //NearestNeighborMapperMatrix(NearestNeighborMapperMatrix const& rOther);

    ///@}

}; // Class NearestNeighborMapperMatrix

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


// inline std::istream & operator >>(std::istream& rIStream,
//                                   NearestNeighborMapperMatrix<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
// {
//     return rIStream;
// }

/// output stream function

// inline std::ostream & operator <<(std::ostream& rOStream,
//                                   const NearestNeighborMapperMatrix<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_NEAREST_NEIGHBOR_MAPPER_MATRIX_H_INCLUDED  defined