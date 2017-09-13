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

#if !defined(KRATOS_UBLAS_MAPPING_MATRIX_BUILDER_H_INCLUDED)
#define KRATOS_UBLAS_MAPPING_MATRIX_BUILDER_H_INCLUDED

// System includes
#include <unordered_map> 

// External includes

// Project includes
#include "includes/define.h"
#include "mapping_matrix_builder.h"

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
class UblasMappingMatrixBuilder : public MappingMatrixBuilder<TSparseSpace, TDenseSpace>
{
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of UblasMappingMatrixBuilder
    KRATOS_CLASS_POINTER_DEFINITION(UblasMappingMatrixBuilder);

    typedef std::unordered_map<int, Node<3>*> EquationIdMapType;

    typedef typename TSparseSpace::DataType TDataType;
    
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;

    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TSparseSpace::MatrixPointerType TSystemMatrixPointerType;

    typedef typename TSparseSpace::VectorPointerType TSystemVectorPointerType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;

    typedef typename TDenseSpace::VectorType LocalSystemVectorType;
    
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > VectorComponentType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    UblasMappingMatrixBuilder() {
        KRATOS_ERROR_IF(TSparseSpace::IsDistributed()) << "WRONG SPACE!" << std::endl;        
    }

    /// Destructor.
    virtual ~UblasMappingMatrixBuilder() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
        
    void UpdateSystemVector(ModelPart& rModelPart,
                            TSystemVectorPointerType pB,
                            const Variable<double>& rVariable) override
    {
        TUpdateSystemVector(rModelPart, pB, rVariable);
    }

    void UpdateSystemVector(ModelPart& rModelPart,
                            TSystemVectorPointerType pB,
                            const VectorComponentType& rVariable) override
    {
        TUpdateSystemVector(rModelPart, pB, rVariable);
    }


    void Update(ModelPart& rModelPart,
                TSystemVectorPointerType pB,
                const Variable<double>& rVariable,
                const Kratos::Flags& MappingOptions,
                const double Factor) override
    {
        TUpdate(rModelPart, pB, rVariable, MappingOptions, Factor);
    }

    void Update(ModelPart& rModelPart,
                TSystemVectorPointerType pB,
                const VectorComponentType& rVariable,
                const Kratos::Flags& MappingOptions,
                const double Factor) override
    {
        TUpdate(rModelPart, pB, rVariable, MappingOptions, Factor);        
    }

    void ResizeAndInitializeVectors(
        TSystemMatrixPointerType& pMdo,
        TSystemVectorPointerType& pQo,
        TSystemVectorPointerType& pQd,
        const unsigned int size_origin,
        const unsigned int size_destination) override
    {
        
        if (!pMdo) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixPointerType pNewMdo = TSparseSpace::CreateEmptyMatrixPointer();
            pMdo.swap(pNewMdo);
        }
        if (!pQo) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewQo = TSparseSpace::CreateEmptyVectorPointer();
            pQo.swap(pNewQo);
        }
        if (!pQd) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemVectorPointerType pNewQd = TSparseSpace::CreateEmptyVectorPointer();
            pQd.swap(pNewQd);
        }


        TSystemMatrixType& Mdo = *pMdo;
        TSystemVectorType& Qo = *pQo;
        TSystemVectorType& Qd = *pQd;

        //resizing the system vectors and matrix
        if (Mdo.size1() == 0) //if the matrix is not initialized
        {
            Mdo.resize(size_destination, size_origin, false);
            // ConstructMatrixStructure(pScheme, A, rElements, rConditions, CurrentProcessInfo); // TODO Jordi is this needed?
        }
        else
        {
            if (Mdo.size1() != size_destination || Mdo.size2() != size_origin)
            {
                KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
                Mdo.resize(size_destination, size_origin, true);
                // ConstructMatrixStructure(pScheme, A, rElements, rConditions, CurrentProcessInfo); // TODO Jordi is this needed?
            }
        }

        if (Qo.size() != size_origin)
            Qo.resize(size_origin, false);

        if (Qd.size() != size_destination)
            Qd.resize(size_destination, false);

    }

    void BuildMappingMatrix(ModelPart::Pointer pModelPart,
                            TSystemMatrixPointerType& pA) override
    {
        TSystemMatrixType& rA = *pA;
        ProcessInfo& r_current_process_info = pModelPart->GetProcessInfo();

        // contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        
        // vector containing the localization in the system of the different terms
        Element::EquationIdVectorType equation_id;

        for (auto& r_cond : pModelPart->GetCommunicator().LocalMesh().Conditions())
        {
            r_cond.CalculateLeftHandSide(LHS_Contribution, r_current_process_info);
            r_cond.EquationIdVector(equation_id, r_current_process_info);

            AssembleLHS(rA, LHS_Contribution, equation_id);
        }
        TSparseSpace::WriteMatrixMarketMatrix("MappingMatrixSerial", rA, false);
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

    // TODO check if those functions do what they are supposed to do!
    // TODO do these two methods even make much of a difference?
    // I have to Initialize the size of the vector again anyway ...
    virtual void ClearData(TSystemMatrixPointerType& pA) override
    {
        TSystemMatrixType& rA = *pA;
        TSparseSpace::ClearData(rA);
    }
    virtual void ClearData(TSystemVectorPointerType& pB) override
    {
        TSystemVectorType& rB = *pB;
        TSparseSpace::ClearData(rB);
    }
    virtual void Clear(TSystemMatrixPointerType& pA) override
    {
        TSparseSpace::Clear(pA);
    }
    virtual void Clear(TSystemVectorPointerType& pB) override
    {
        TSparseSpace::Clear(pB);
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
        buffer << "UblasMappingMatrixBuilder";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const { rOStream << "UblasMappingMatrixBuilder"; }

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

    template< class TVarType>
    void TUpdateSystemVector(ModelPart& rModelPart,
                             TSystemVectorPointerType pB,
                             const TVarType& rVariable) 
    {
        // Jordi how to do this conversion nicely?
        TSystemVectorType& r_b = *pB;
        int index = 0;
        for (auto& node : rModelPart.Nodes())
        {
            // Jordi is this ok? => Done differently in some BuilderAndSolvers
            r_b[index] = node.FastGetSolutionStepValue(rVariable);
            ++index;
        }
        TSparseSpace::WriteMatrixMarketVector("UpdateSystemVector", r_b);
    }

    template< class TVarType>
    void TUpdate(ModelPart& rModelPart,
                 TSystemVectorPointerType pB,
                 const TVarType& rVariable,
                 const Kratos::Flags& MappingOptions,
                 const double Factor) 
    {
        // Jordi how to do this conversion nicely?
        TSystemVectorType& r_b = *pB;
        TSparseSpace::WriteMatrixMarketVector("Update", r_b);        

        int index = 0;
        for (auto& node : rModelPart.Nodes())
        {
            if (MappingOptions.Is(MapperFlags::ADD_VALUES))
            {
                // Jordi is this ok? => Done differently in some BuilderAndSolvers
                node.FastGetSolutionStepValue(rVariable) += r_b[index] * Factor;
            }
            else
            {
                // Jordi is this ok? => Done differently in some BuilderAndSolvers
                node.FastGetSolutionStepValue(rVariable) = r_b[index] * Factor;
            }
            ++index;
        }
    }

    void AssembleLHS(
        TSystemMatrixType& A,
        LocalSystemMatrixType& LHS_Contribution,
        Condition::EquationIdVectorType& EquationId
    )
    {
        const unsigned int local_size = LHS_Contribution.size1();
        int equation_id_destination;
        int equation_id_origin;
        
        for (unsigned int i = 0; i < local_size ; ++i)
        {
            for (unsigned int j = 0; j < local_size ; ++j)
            {
                equation_id_destination = EquationId[i];
                equation_id_origin = EquationId[i + local_size];

                // std::cout << "equation_id_destination: " << equation_id_destination << " ;; equation_id_origin: " 
                //           << equation_id_origin << " ;; LHS_Contribution(i,j): " << LHS_Contribution(i,j) << std::endl;

                A(equation_id_destination, equation_id_origin) += LHS_Contribution(i,j);
            }
        }


        // unsigned int local_size = LHS_Contribution.size1();

        // for (unsigned int i_local = 0; i_local < local_size; i_local++)
        // {
        //     unsigned int i_global = EquationId[i_local];

        //     for (unsigned int j_local = 0; j_local < local_size; j_local++)
        //     {
        //         unsigned int j_global = EquationId[j_local];

        //         A(i_global, j_global) += LHS_Contribution(i_local, j_local);
        //     }
        // }

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
    // UblasMappingMatrixBuilder &operator=(UblasMappingMatrixBuilder const &rOther) {}

    /// Copy constructor.
    // UblasMappingMatrixBuilder(UblasMappingMatrixBuilder const &rOther) {}

    ///@}

}; // Class UblasMappingMatrixBuilder

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream &operator>>(std::istream &rIStream,
//                                 UblasMappingMatrixBuilder &rThis) {}

// /// output stream function
// inline std::ostream &operator<<(std::ostream &rOStream,
//                                 const UblasMappingMatrixBuilder &rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_UBLAS_MAPPING_MATRIX_BUILDER_H_INCLUDED  defined
