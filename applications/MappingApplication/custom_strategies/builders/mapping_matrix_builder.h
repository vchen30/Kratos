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
#include <unordered_map> 

// External includes

// Project includes
#include "includes/define.h"
#include "mapping_application_variables.h"

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
          class TDenseSpace  // = DenseSpace<double>,
          >
class MappingMatrixBuilder
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MappingMatrixBuilder
    KRATOS_CLASS_POINTER_DEFINITION(MappingMatrixBuilder);

    typedef std::unordered_map<int, Node<3>*> EquationIdMapType;

    typedef typename TSparseSpace::DataType TDataType;
    
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;

    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TSparseSpace::MatrixPointerType TSystemMatrixPointerType;

    typedef typename TSparseSpace::VectorPointerType TSystemVectorPointerType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;

    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MappingMatrixBuilder() {

      TSparseSpace::WhatAmI();
    }

    /// Destructor.
    virtual ~MappingMatrixBuilder() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /** 
    This function set up the structure of the system, i.e. the NodeSet,
    which relates the nodes to the equation Ids
    */
    void SetUpSystem(ModelPart& rModelPart, 
                     EquationIdMapType& EquationIdNodeMap)
    {
        // EquationIdNodeMap.clear();

        // these ids are the positions in the global vectors and matrix
        int equation_id = GetStartEquationId(rModelPart);
        std::cout << "Start Equation ID = " << equation_id << std::endl;

        for (auto& node : rModelPart.GetCommunicator().LocalMesh().Nodes())
        {
        // EquationIdNodeMap.emplace(equation_id, &node); // TODO check if this is a pointer
        node.SetValue(MAPPING_MATRIX_EQUATION_ID, equation_id); // TODO replace with sth faster?
        ++equation_id;
        }

        for (auto& node : rModelPart.GetCommunicator().LocalMesh().Nodes())
        {
            KRATOS_WATCH(node.GetValue(MAPPING_MATRIX_EQUATION_ID));
        }
    }

    virtual void UpdateSystemVector(ModelPart& rModelPart,
                        TSystemVectorPointerType b,
                        const Variable<double>& rVariable)
    {
        // Jordi how to do this conversion nicely?
        TSystemVectorType& r_b = *b;
        int index = 0;
        for (auto& node : rModelPart.Nodes())
        {
            // Jordi how to make this work?
            // r_b[index] = node.FastGetSolutionStepValue(rVariable);
            ++index;
        }
    }

    virtual void Update(ModelPart& rModelPart,
                        TSystemVectorPointerType b,
                        const Variable<double>& rVariable)
    {
        // Jordi how to do this conversion nicely?
        TSystemVectorType& r_b = *b;
        int index = 0;
        for (auto& node : rModelPart.Nodes())
        {
            // Jordi how to make this work?
            // node.FastGetSolutionStepValue(rVariable) = r_b[index];
            ++index;
        }
    }

    // /**
    // This functions build the LHS (aka the Mapping Matrix Mdo) of the mapping problem
    //  */
    // virtual void BuildLHS(ModelPart& rModelPart,
    //                       TSystemMatrixType& A) override
    // {   // funciton copied from "residualbased_block_builder_and_solver.h"
    // //     KRATOS_TRY
    // //     if (!pScheme)
    // //         KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

    // //     //getting the array of the conditions
    // //     const int nconditions = static_cast<int>(rModelPart.Conditions().size());

    // //     ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
    // //     ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

    // //     //contributions to the system
    // //     LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);

    // //     //vector containing the localization in the system of the different
    // //     //terms
    // //     Element::EquationIdVectorType EquationId;

    // //     // assemble all elements
    // //     double start_build = OpenMPUtils::GetCurrentTime();

    // //     #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
    // //     {
    // //         #pragma omp for  schedule(guided, 512)
    // //         for (int k = 0; k < nconditions; k++)
    // //         {
    // //             ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

    // //             //detect if the condition is active or not. If the user did not make any choice the element
    // //             //is active by default
    // //             bool condition_is_active = true;
    // //             if ((it)->IsDefined(ACTIVE))
    // //                 condition_is_active = (it)->Is(ACTIVE);

    // //             if (condition_is_active)
    // //             {
    // //                 //calculate condition contribution
    // //                 pScheme->Condition_Calculate_LHS_Contribution(*(it.base()), LHS_Contribution, EquationId, CurrentProcessInfo);

    // //                 //assemble the condition contribution
    // // #ifdef USE_LOCKS_IN_ASSEMBLY
    // //                 AssembleLHS(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    // // #else
    // //                 AssembleLHS(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    // // #endif

    // //                 // clean local condition memory
    // //                 pScheme->CleanMemory(*(it.base()));
    // //             }
    // //         }
    // //     }

    // //     double stop_build = OpenMPUtils::GetCurrentTime();
    // //     if (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)
    // //         std::cout << "build time: " << stop_build - start_build << std::endl;

    // //     //for (int i = 0; i < A_size; i++)
    // //     //    omp_destroy_lock(&lock_array[i]);
    // //     if (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)
    // //     {
    // //         KRATOS_WATCH("finished parallel building");
    // //     }

    // //     KRATOS_CATCH("")
    // }

    /**
    This functions build the RHS (aka the vector of nodal quantities 
    for a given variable) of the mapping problem
     */
    // template <typename T>
    // void UpdateRHS(EquationIdMapType& EquationIdNodeMap,
    //               TSystemVectorType& b
    //               const Variable< T >& rVariable)
    // {
    //     for (auto const &entry : EquationIdNodeMap)
    //     {
    //         b[entry.first] = entry.second.FastGetCurrentSolutionStepValue(rVariable);
    //     }

    //     GlobalAssembleVector(b); // TODO is this needed?
    // }

    // virtual void SetUpDofSet(
    //     ModelPart& r_model_part
    // ) override
    // {
    // }

    // /**
    // organises the dofset in order to speed up the building phase
    //  */
    // virtual void SetUpSystem(
    //     ModelPart& r_model_part
    // ) override
    // {
    // }

    // virtual void ResizeAndInitializeVectors(
    //     TSystemMatrixPointerType& pA,
    //     TSystemVectorPointerType& pDx,
    //     TSystemVectorPointerType& pb,
    //     ElementsArrayType& rElements,
    //     ConditionsArrayType& rConditions,
    //     ProcessInfo& CurrentProcessInfo
    // ) override
    // {
    // }

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

    /**
    This functions computes the start equationId for to local nodes
    aka the equation id of the global system
    */
    virtual int GetStartEquationId(ModelPart& rModelPart)
    {   
        return 0;
    }

    /**
    This function does nothing in the serial case
    It exists such that the BuildRHS only exists in this
    class and can therefore be templated
    */
    // virtual void GlobalAssembleVector(TSystemVectorType& b) {}

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
    // MappingMatrixBuilder &operator=(MappingMatrixBuilder const &rOther) {}

    /// Copy constructor.
    // MappingMatrixBuilder(MappingMatrixBuilder const &rOther) {}

    ///@}

}; // Class MappingMatrixBuilder

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream &operator>>(std::istream &rIStream,
//                                 MappingMatrixBuilder &rThis) {}

// /// output stream function
// inline std::ostream &operator<<(std::ostream &rOStream,
//                                 const MappingMatrixBuilder &rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED  defined
