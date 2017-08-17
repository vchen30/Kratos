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

#if !defined(KRATOS_MAPPER_MATRIX_BASED_H_INCLUDED )
#define  KRATOS_MAPPER_MATRIX_BASED_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"
#include "custom_strategies/strategies/mapper_strategy.h"


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
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
>
class MapperMatrixBased : public Mapper
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of MapperMatrixBased
    KRATOS_CLASS_POINTER_DEFINITION(MapperMatrixBased);



    ///@}
    ///@name Life Cycle
    ///@{

    MapperMatrixBased(ModelPart& rModelPartOrigin, ModelPart& rModelpartDestination,
                          Parameters rJsonParameters) : Mapper(
                            rModelPartOrigin, rModelpartDestination, rJsonParameters)
    {
        // InitializeMapperStrategy();
        // MapperStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPartOrigin, rModelpartDestination, rJsonParameters);


        // @Jordi should the BuilderAndSolver be a member of this class or of the Communicator?
        // In case of mortar we need two (one for Mdo and one for Mdd). The one for Mdd is specific to Mortar and should 
        // imo not be in the communicator. In order to be consistent it might make sense to have also the Mdo-BuilderAndSolver as a
        // member of this class.
        // If we don't use a BuilderAndSolver for Mdo then the above is obsolete.
        // In any case, if the BuilderAndSolver is member of the Mapper, how can we destinguish btw serial and parallel?
        // mpCommunicator->InitializeBuilderAndSolverForMdo();

        // mpMapperCommunicator->InitilizeMappingMatrixUtility(mpMappingMatrixUtility);
    }

    /// Destructor.
    virtual ~MapperMatrixBased() { }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius) override
    {

    }

    /* This function maps from Origin to Destination */
    void Map(const Variable<double>& rOriginVariable,
             const Variable<double>& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        
    }

    /* This function maps from Origin to Destination */
    void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
             const Variable< array_1d<double, 3> >& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        
    }

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
        return "MapperMatrixBased";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperMatrixBased";
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

    // MapperStrategy::Pointer mpMapperStrategy;

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

//     void InitializeMapperStrategy()
//     {
// #ifdef KRATOS_USING_MPI // mpi-parallel compilation
//         if (mCommSize > 1)
//         {
//             mpMapperStrategy = MapperStrategy<TrilinosSparseSpaceType, 
//             LocalSpaceType, TrilinosLinearSolverType>(mrModelPartOrigin, 
//                             mrModelPartDestination);
//         }
//         else
//         {
//             mpMapperStrategy = MapperStrategy<SerialSparseSpaceType, 
//             LocalSpaceType, SerialLinearSolverType>(mrModelPartOrigin, 
//                             mrModelPartDestination);
//         }

// #else // serial compilation
//         mpMapperStrategy = MapperStrategy<SerialSparseSpaceType, 
//         LocalSpaceType, SerialLinearSolverType>(mrModelPartOrigin, 
//                         mrModelPartDestination);
// #endif


// }

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
    // MapperMatrixBased<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(MapperMatrixBased<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    //MapperMatrixBased(MapperMatrixBased const& rOther);

    ///@}

}; // Class MapperMatrixBased

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


// inline std::istream & operator >>(std::istream& rIStream,
//                                   MapperMatrixBased<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
// {
//     return rIStream;
// }

// /// output stream function

// inline std::ostream & operator <<(std::ostream& rOStream,
//                                   const MapperMatrixBased<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_MAPPER_MATRIX_BASED_H_INCLUDED  defined