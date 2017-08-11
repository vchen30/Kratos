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
//

#if !defined(KRATOS_MAPPER_STRATEGY_H_INCLUDED)
#define KRATOS_MAPPER_STRATEGY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/solving_strategy.h"

namespace Kratos
{
///@addtogroup MappingApplication
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
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class MapperStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperStrategy
    KRATOS_CLASS_POINTER_DEFINITION(MapperStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    
    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperStrategy(ModelPart& model_part_origin,
        ModelPart& model_part_destination,
        typename TBuilderAndSolverType::Pointer pMappingMatrixBuilder
    ) 
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part_destination), 
        mrModelPartOrigin(model_part_origin)
    {
        mpMappingMatrixBuilder = pMappingMatrixBuilder;

        // Initializing the Mapping Matrix and vector of quantities
        mpMdo = TSparseSpace::CreateEmptyMatrixPointer();
        mpQo = TSparseSpace::CreateEmptyVectorPointer();
        mpQd = TSparseSpace::CreateEmptyVectorPointer();
    }

    // constructor for Mortar, takes the BuilderAndSolver for Mdd
    MapperStrategy(ModelPart& model_part_origin,
        ModelPart& model_part_destination,
        typename TBuilderAndSolverType::Pointer pMappingMatrixBuilder,
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver
    ) 
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part_destination),
        mrModelPartOrigin(model_part_origin)
    {
        mIsMortar = true;

        mpMappingMatrixBuilder = pMappingMatrixBuilder;
        mpBuilderAndSolver = pBuilderAndSolver;

        // Initializing the Mapping Matrix and vector of quantities
        mpMdo = TSparseSpace::CreateEmptyMatrixPointer();
        mpMdd = TSparseSpace::CreateEmptyMatrixPointer();
        mpQo = TSparseSpace::CreateEmptyVectorPointer();
        mpQd = TSparseSpace::CreateEmptyVectorPointer();
        mpQtmp = TSparseSpace::CreateEmptyVectorPointer();
    }

    /// Destructor.
    virtual ~MapperStrategy() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /**
     * This function computes the mapping matrix
     */    
    void Initialize() override
    {
        mpMappingMatrixBuilder->BuildLHS(BaseType::GetModelPart(), mpMdo);

        // assemble the mass matrix for the mortar mapper (M_dd)
        if (mIsMortar)
        {
            mpBuilderAndSolver->BuildLHS(mrModelPartOrigin, mpMdd);
        }
    }

    template <typename T>
    void InitializeMappingStep(const Variable< T >& rOriginVariable,
                               const Variable< T >& rDestinationVariable,
                               Kratos::Flags MappingOptions,
                               const bool InverseOperation = false)
    {
        if (InverseOperation) // for conservative mapping
        {
            mpMappingMatrixBuilder->BuildRHS(BaseType::GetModelPart(), mpQd);
        }
        else
        {
            mpMappingMatrixBuilder->BuildRHS(mrModelPartOrigin, mpQo);
        }
        
    }

    template <typename T>
    void FinalizeMappingStep(const Variable< T >& rOriginVariable,
                             const Variable< T >& rDestinationVariable,
                             Kratos::Flags MappingOptions,
                             const bool InverseOperation = false)
    {
        if (InverseOperation) // for conservative mapping
        {
            Update(mrModelPartOrigin, mpMappingMatrixBuilder->GetDofSet(), 
                             mpMdo, mpQo, mpQd);
        }
        else
        {
            Update(BaseType::GetModelPart(), mpMappingMatrixBuilder->GetDofSet(), 
                             mpMdo, mpQd, mpQo);
        }
    }


    bool SolveSolutionStep(const bool InverseOperation = false) override
    {   
        if (InverseOperation) // for conservative mapping
        {
           if (mIsMortar)
            {
                mpBuilderAndSolver->SystemSolve(mpMdd, mpQtmp, mpQd); // Jordi the trilinos call also wants a modelpart!
                TSparseSpace::TransposeMult(mpMdo, mpQtmp, mpQo);
            }
            else
            {
                TSparseSpace::TransposeMult(mpMdo, mpQd, mpQo);
            } 
        }
        else
        {
            if (mIsMortar)
            {
                TSparseSpace::Mult(mpMdo, mpQo, mpQtmp);
                mpBuilderAndSolver->SystemSolve(mpMdd, mpQd, mpQtmp); // Jordi the trilinos call also wants a modelpart!
                // if this turns out to be a problem we can pass the BuilderAndSolver to the MappingMatrixBuilder to call 
                // the correct function, bcs this is different for serial and Trilinos
                // I would try to avoid duplicating the MapperStrategy at all cost!!!
            }
            else
            {
                TSparseSpace::Mult(mpMdo, mpQo, mpQd);
            }
        }

        return true;
    }

    template <typename T>
    void Map(const Variable< T >& rOriginVariable,
             const Variable< T >& rDestinationVariable,
             Kratos::Flags MappingOptions,
             const bool InverseOperation = false)
    {
        InitializeMappingStep(rOriginVariable, rDestinationVariable, 
                              MappingOptions, InverseOperation);

        SolveSolutionStep(InverseOperation);

        FinalizeMappingStep(rOriginVariable, rDestinationVariable, 
                            MappingOptions, InverseOperation);
    }



    ///@}
    ///@name Access
    ///@{

    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mMdo = *mpMdo;

        return mMdo;
    }

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
        buffer << "MapperStrategy";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const { rOStream << "MapperStrategy"; }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const {}

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

    ModelPart& mrModelPartOrigin;

    typename TBuilderAndSolverType::Pointer mpMappingMatrixBuilder; // TODO change type...?
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver; // needed for Mortar
                    
    TSystemVectorPointerType mpQo;
    TSystemVectorPointerType mpQd;
    TSystemVectorPointerType mpQtmp; // for Mortar (needed bcs Trilinos cannot multiply in place)
    TSystemMatrixPointerType mpMdo;
    TSystemMatrixPointerType mpMdd; // for Mortar

    bool mIsMortar = false; // needed to check if additional operations (for mortar-type mappers are required)

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
    MapperStrategy &operator=(MapperStrategy const &rOther) {}

    /// Copy constructor.
    MapperStrategy(MapperStrategy const &rOther) {}

    ///@}

}; // Class MapperStrategy

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_MAPPER_STRATEGY_H_INCLUDED  defined
