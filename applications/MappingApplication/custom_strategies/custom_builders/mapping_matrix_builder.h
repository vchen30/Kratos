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


#if !defined(KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED)
#define KRATOS_MAPPING_MATRIX_BUILDER_H_INCLUDED

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "base_mapping_matrix_builder.h"

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
class MappingMatrixBuilder : public BaseMappingMatrixBuilder<TSparseSpace, TDenseSpace>
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

    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > VectorComponentType;

    typedef ModelPart::NodeIterator NodeIterator;
    typedef ModelPart::ConditionIterator ConditionIterator;
    typedef Condition::EquationIdVectorType EquationIdVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MappingMatrixBuilder(const int EchoLevel) : MappingMatrixBuilder<TSparseSpace, TDenseSpace>(EchoLevel)
    {
        KRATOS_ERROR_IF(TSparseSpace::IsDistributed()) << "WRONG SPACE!" << std::endl;
    }

    /// Destructor.
    virtual ~MappingMatrixBuilder() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void UpdateSystemVector(ModelPart& rModelPart,
                            TSystemVectorType& rB,
                            const Variable<double>& rVariable) override
    {
        TUpdateSystemVector(rModelPart, rB, rVariable);
    }

    void UpdateSystemVector(ModelPart& rModelPart,
                            TSystemVectorType& rB,
                            const VectorComponentType& rVariable) override
    {
        TUpdateSystemVector(rModelPart, rB, rVariable);
    }


    void Update(ModelPart& rModelPart,
                TSystemVectorType& rB,
                const Variable<double>& rVariable,
                const Kratos::Flags& MappingOptions,
                const double Factor) override
    {
        TUpdate(rModelPart, rB, rVariable, MappingOptions, Factor);
    }

    void Update(ModelPart& rModelPart,
                TSystemVectorType& rB,
                const VectorComponentType& rVariable,
                const Kratos::Flags& MappingOptions,
                const double Factor) override
    {
        TUpdate(rModelPart, rB, rVariable, MappingOptions, Factor);
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
            // ConstructMatrixStructure(pScheme, A, rElements, rConditions, CurrentProcessInfo);
            ConstructMatrixStructure();
        }
        else
        {
            if (Mdo.size1() != size_destination || Mdo.size2() != size_origin)
            {
                KRATOS_WATCH("it should not come here!!!!!!!! ... this is SLOW");
                Mdo.resize(size_destination, size_origin, true);
                // ConstructMatrixStructure(pScheme, A, rElements, rConditions, CurrentProcessInfo);
                ConstructMatrixStructure();
            }
        }

        if (Qo.size() != size_origin) Qo.resize(size_origin, false);

        if (Qd.size() != size_destination) Qd.resize(size_destination, false);

    }

    void BuildMappingMatrix(ModelPart& rModelPart,
                            TSystemMatrixType& rA) override
    {
        // contributions to the system
        LocalSystemMatrixType mapper_local_system = LocalSystemMatrixType(0,0);

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        EquationIdVectorType equation_id;

        const int num_conditions = rModelPart.NumberOfConditions();
        ConditionIterator it_begin = rModelPart.ConditionsBegin();

        // #pragma omp parallel for // TODO check if this works, i.e. if I write the same positions several times!
        for(int i = 0; i < num_conditions; i++)
        {
            ConditionIterator cond_it = it_begin + i;



            cond_it->EquationIdVector(equation_id, r_current_process_info);
            cond_it->CalculateLeftHandSide(mapper_local_system, r_current_process_info);

            Assemble(rA, mapper_local_system, equation_id);
        }

        if (this->mEchoLevel >= 1) TSparseSpace::WriteMatrixMarketMatrix("MappingMatrixSerial", rA, false); // TODO change Level to sth higher later
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
                             TSystemVectorType& rB,
                             const TVarType& rVariable)
    {
        const int num_nodes = rModelPart.NumberOfNodes();
        NodeIterator it_begin = rModelPart.NodesBegin();

        #pragma omp parallel for // Don't modify, this is best suitable for this case
        for (int i = 0; i<num_nodes; i++)
        {
            NodeIterator it = it_begin + i;
            rB[i] = it->FastGetSolutionStepValue(rVariable);
        }

        if (this->mEchoLevel >= 1) TSparseSpace::WriteMatrixMarketVector("UpdateSystemVector", rB); // TODO change Level to sth higher later
    }

    template< class TVarType>
    void TUpdate(ModelPart& rModelPart,
                 TSystemVectorType& rB,
                 const TVarType& rVariable,
                 const Kratos::Flags& MappingOptions,
                 const double Factor)
    {
        if (this->mEchoLevel >= 1) TSparseSpace::WriteMatrixMarketVector("Update", rB); // TODO change Level to sth higher later

        const int num_nodes = rModelPart.NumberOfNodes();
        NodeIterator it_begin = rModelPart.NodesBegin();

        #pragma omp parallel for // Don't modify, this is best suitable for this case
        for (int i = 0; i<num_nodes; i++)
        {
            NodeIterator it = it_begin + i;

            if (MappingOptions.Is(MapperFlags::ADD_VALUES))
                it->FastGetSolutionStepValue(rVariable) += rB[i] * Factor;
            else
                it->FastGetSolutionStepValue(rVariable) = rB[i] * Factor;
        }
    }

    void Assemble(TSystemMatrixType& rA,
                  const LocalSystemMatrixType& rLocalSystem,
                  const EquationIdVectorType& rEquationID)
    {
        /* The format of "rLocalSystem" is:
        1. Row: Weight
        2. Row: MAPPING_MATRIX_ID on Destination
        3. Row: MAPPING_MATRIX_ID on Origin */

        const unsigned int local_size = rLocalSystem.size2(); // for

        // No openmp here, is done at higher level! => but maybe protect the writing with atomic?
        for (unsigned int i = 0; i < local_size; ++i)
            rA(rLocalSystem(1,i) , rLocalSystem(2,i)) += rLocalSystem(0,i);  // Big Question: "=" or "+=" ??? TODO

        // const unsigned int local_size = LHS_Contribution.size1();
        // int equation_id_destination;
        // int equation_id_origin;

        // for (unsigned int i = 0; i < local_size ; ++i)
        // {
        //     for (unsigned int j = 0; j < local_size ; ++j)
        //     {
        //         equation_id_destination = EquationId[i];
        //         equation_id_origin = EquationId[i + local_size];

        //         // std::cout << "equation_id_destination: " << equation_id_destination << " ;; equation_id_origin: "
        //         //           << equation_id_origin << " ;; LHS_Contribution(i,j): " << LHS_Contribution(i,j) << std::endl;

        //         A(equation_id_destination, equation_id_origin) += LHS_Contribution(i,j);
        //     }
        // }


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

    void ConstructMatrixStructure()
    {
        // TODO implement this function
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
