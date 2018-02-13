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
        this->mInterfaceParameters = Parameters( R"(
        {
            "mapper_condition_name" : "NearestNeighborMapperCondition",
            "use_nodes"      : true
        }  )" );
        
        this->GenerateInterfaceModelPart();

        this->mpMapperCommunicator->InitializeOrigin(MapperUtilities::Node_Coords);
        this->mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);        
        this->mpMapperCommunicator->Initialize();

        this->ComputeMappingMatrix();
    }

    /// Destructor.
    virtual ~NearestNeighborMapperMatrix() { }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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

    void InitializeInverseMapper() override
    {
        this->mpInverseMapper = Mapper::Pointer( new NearestNeighborMapperMatrix(this->mrModelPartDestination,
                                                                                 this->mrModelPartOrigin,
                                                                                 this->mJsonParameters) );
    }

    void ExchangeInterfaceGeometryData() override
    {

        // Creating the function pointers for the InterfaceObjects
        // auto function_pointer_origin = std::bind(&GetValueOfNode,
        //                                std::placeholders::_1,
        //                                std::placeholders::_2);

        // auto function_pointer_destination = std::bind(&SetValueOfNode,
        //                                     std::placeholders::_1,
        //                                     std::placeholders::_2);

        this->mpMapperCommunicator->TransferVariableDataNEW(&GetGeometryInformation, // TODO pass by ref or not?
            &SetValueOfNode);
    }

    static Vector GetGeometryInformation(InterfaceObject* pInterfaceObject, //TODO const
                          const std::vector<double>& rShapeFunctionValues)
    {
        Node<3>* p_base_node = static_cast<InterfaceNode*>(pInterfaceObject)->pGetBase();
        KRATOS_ERROR_IF_NOT(p_base_node) << "Base Pointer is nullptr!!!" << std::endl;
        Vector equation_ids(1);
        equation_ids[0] = p_base_node->GetValue(MAPPING_MATRIX_EQUATION_ID);
        return equation_ids;
    }


    static void SetValueOfNode(InterfaceObject* pInterfaceObject,
                               const Vector& rValue)
    {
        Node<3>* p_base_node = dynamic_cast<InterfaceNode*>(pInterfaceObject)->pGetBase();
        KRATOS_ERROR_IF_NOT(p_base_node) << "Base Pointer is nullptr!!!" << std::endl;

        p_base_node->SetValue(MAPPER_NEIGHBOR_INFORMATION, rValue);
    }


    // static void SetValueOfNode(InterfaceObject* pInterfaceObject,
    //                            const Vector& rValue)
    // {
    //     // BaseMapperCondition* p_base_condition = dynamic_cast<InterfaceNode*>(pInterfaceObject)->pGetBase();
    //     // KRATOS_ERROR_IF_NOT(p_base_condition) << "Base Pointer is nullptr!!!" << std::endl;

    //     // // p_base_node->SetValue(MAPPING_MATRIX_EQUATION_ID_VECTOR, rValue);
    //     // std::vector<int> equation_ids(rValue.size());
    //     // for (unsigned int i=0; i< rValue.size(); ++i)
    //     //     equation_ids[i] = rValue[i];

    //     // p_base_condition->SetNeighborIDVector(equation_ids);
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