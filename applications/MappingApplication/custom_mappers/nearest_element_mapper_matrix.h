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

#if !defined(KRATOS_NEAREST_ELEMENT_MAPPER_MATRIX_H_INCLUDED )
#define  KRATOS_NEAREST_ELEMENT_MAPPER_MATRIX_H_INCLUDED

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
class NearestElementMapperMatrix : public MapperMatrixBased<TMappingMatrixBuilder, TLinearSolver>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NearestElementMapperMatrix
    KRATOS_CLASS_POINTER_DEFINITION(NearestElementMapperMatrix);
    
    typedef TMappingMatrixBuilder TMappingMatrixBuilderType;

    typedef typename TMappingMatrixBuilderType::TDataType TDataType;
    
    typedef typename TMappingMatrixBuilderType::TSystemMatrixType TSystemMatrixType;

    typedef typename TMappingMatrixBuilderType::TSystemVectorType TSystemVectorType;

    typedef typename TMappingMatrixBuilderType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename TMappingMatrixBuilderType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename TMappingMatrixBuilderType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename TMappingMatrixBuilderType::LocalSystemVectorType LocalSystemVectorType;

    typedef size_t SizeType;

    typedef Node<3> NodeType;

    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    NearestElementMapperMatrix(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                                Parameters rJsonParameters) : MapperMatrixBased<TMappingMatrixBuilder, TLinearSolver>(
                                i_model_part_origin, i_model_part_destination, rJsonParameters)
    {
        this->mInterfaceParameters = Parameters( R"(
        {
            "mapper_condition_name" : "NearestElementMapperCondition",
            "use_nodes"      : true
        }  )" ); // Nodes are used bcs the geometry of the destination is not necessary for this mapper
        
        this->GenerateInterfaceModelPart();

        this->mpMapperCommunicator->InitializeOrigin(MapperUtilities::Node_Coords);
        this->mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);        
        this->mpMapperCommunicator->Initialize();

        this->ComputeMappingMatrix();
    }

    /// Destructor.
    virtual ~NearestElementMapperMatrix() { }

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
        return "NearestElementMapperMatrix";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestElementMapperMatrix";
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
        this->mpInverseMapper = Mapper::Pointer( new NearestElementMapperMatrix(this->mrModelPartDestination,
                                                                                this->mrModelPartOrigin,
                                                                                this->mJsonParameters) );
    }

    void ExchangeInterfaceGeometryData() override
    {
        this->mpMapperCommunicator->TransferVariableDataNEW(&GetGeometryInformation, // TODO pass by ref or not?
            &SetValueOfNode);
    }

    static Vector GetGeometryInformation(InterfaceObject* pInterfaceObject, //TODO const
                          const std::vector<double>& rShapeFunctionValues)
    {
        GeometricalObject* p_base_geometrical_object = static_cast<InterfaceGeometryObject*>(pInterfaceObject)->pGetBase();
        KRATOS_ERROR_IF_NOT(p_base_geometrical_object) << "Base Pointer is nullptr!!!" << std::endl;
        
        GeometryType& r_geom = p_base_geometrical_object->GetGeometry();
        const SizeType num_points = r_geom.PointsNumber();

        Vector geom_information(2*num_points);

        for (SizeType i=0; i<num_points; ++i)
        {
            geom_information[2*i] = r_geom.GetPoint(i).GetValue(MAPPING_MATRIX_EQUATION_ID);
            geom_information[2*i+1] = rShapeFunctionValues[i];
        }

        return geom_information;
    }

    static void SetValueOfNode(InterfaceObject* pInterfaceObject,
                               const Vector& rValue)
    {
        Node<3>* p_base_node = dynamic_cast<InterfaceNode*>(pInterfaceObject)->pGetBase();
        KRATOS_ERROR_IF_NOT(p_base_node) << "Base Pointer is nullptr!!!" << std::endl;

        p_base_node->SetValue(MAPPER_NEIGHBOR_INFORMATION, rValue);
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
    // NearestElementMapperMatrix<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(NearestElementMapperMatrix<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    //NearestElementMapperMatrix(NearestElementMapperMatrix const& rOther);

    ///@}

}; // Class NearestElementMapperMatrix

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


// inline std::istream & operator >>(std::istream& rIStream,
//                                   NearestElementMapperMatrix<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
// {
//     return rIStream;
// }

/// output stream function

// inline std::ostream & operator <<(std::ostream& rOStream,
//                                   const NearestElementMapperMatrix<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_NEAREST_ELEMENT_MAPPER_MATRIX_H_INCLUDED  defined