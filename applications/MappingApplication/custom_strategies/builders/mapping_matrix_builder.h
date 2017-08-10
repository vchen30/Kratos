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

#if !defined(KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED)
#define KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

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
class MappingMatrixBuilder : public BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MappingMatrixBuilder
    KRATOS_CLASS_POINTER_DEFINITION(MappingMatrixBuilder);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MappingMatrixBuilder() {}

    /// Destructor.
    virtual ~MappingMatrixBuilder() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /**
    This functions build the LHS (aka the Mapping Matrix Mdo) of the mapping problem
     */
    virtual void BuildLHS(typename TSchemeType::Pointer pScheme,
                          ModelPart& rModelPart,
                          TSystemMatrixType& A) override
    {   // funciton copied from "residualbased_block_builder_and_solver.h"
    //     KRATOS_TRY
    //     if (!pScheme)
    //         KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

    //     //getting the array of the conditions
    //     const int nconditions = static_cast<int>(rModelPart.Conditions().size());

    //     ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
    //     ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

    //     //contributions to the system
    //     LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

    //     //vector containing the localization in the system of the different
    //     //terms
    //     Element::EquationIdVectorType EquationId;

    //     // assemble all elements
    //     double start_build = OpenMPUtils::GetCurrentTime();

    //     #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
    //     {
    //         #pragma omp for  schedule(guided, 512)
    //         for (int k = 0; k < nconditions; k++)
    //         {
    //             ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

    //             //detect if the condition is active or not. If the user did not make any choice the element
    //             //is active by default
    //             bool condition_is_active = true;
    //             if ((it)->IsDefined(ACTIVE))
    //                 condition_is_active = (it)->Is(ACTIVE);

    //             if (condition_is_active)
    //             {
    //                 //calculate condition contribution
    //                 pScheme->Condition_Calculate_LHS_Contribution(*(it.base()), LHS_Contribution, EquationId, CurrentProcessInfo);

    //                 //assemble the condition contribution
    // #ifdef USE_LOCKS_IN_ASSEMBLY
    //                 AssembleLHS(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    // #else
    //                 AssembleLHS(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    // #endif

    //                 // clean local condition memory
    //                 pScheme->CleanMemory(*(it.base()));
    //             }
    //         }
    //     }

    //     double stop_build = OpenMPUtils::GetCurrentTime();
    //     if (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)
    //         std::cout << "build time: " << stop_build - start_build << std::endl;

    //     //for (int i = 0; i < A_size; i++)
    //     //    omp_destroy_lock(&lock_array[i]);
    //     if (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
    //     {
    //         KRATOS_WATCH("finished parallel building");
    //     }

    //     KRATOS_CATCH("")
    }

    /**
    This functions build the RHS (aka the vector of nodal quantities 
    for a given variable) of the mapping problem
     */
    virtual void BuildRHS(typename TSchemeType::Pointer pScheme,
                          ModelPart& rModelPart,
                          TSystemVectorType& b) override
    {


    }

    virtual void SetUpDofSet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& r_model_part
    ) override
    {
    }

    /**
    organises the dofset in order to speed up the building phase
     */
    virtual void SetUpSystem(
        ModelPart& r_model_part
    ) override
    {
    }

    virtual void ResizeAndInitializeVectors(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixPointerType& pA,
        TSystemVectorPointerType& pDx,
        TSystemVectorPointerType& pb,
        ElementsArrayType& rElements,
        ConditionsArrayType& rConditions,
        ProcessInfo& CurrentProcessInfo
    ) override
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MappingMatrixBuilder";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const { rOStream << "MappingMatrixBuilder"; }

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
    MappingMatrixBuilder &operator=(MappingMatrixBuilder const &rOther) {}

    /// Copy constructor.
    MappingMatrixBuilder(MappingMatrixBuilder const &rOther) {}

    ///@}

}; // Class MappingMatrixBuilder

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                MappingMatrixBuilder &rThis) {}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const MappingMatrixBuilder &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED  defined
