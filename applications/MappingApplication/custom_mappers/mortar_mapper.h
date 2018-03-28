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

#if !defined(KRATOS_MORTAR_MAPPER_H_INCLUDED )
#define  KRATOS_MORTAR_MAPPER_H_INCLUDED

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
template < class TMapperCommunicator,
           class TMappingMatrixBuilder,
           class TBuilderAndSolver
>
class MortarMapper : public MapperMatrixBased< TMapperCommunicator, TMappingMatrixBuilder >
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of MortarMapper
    KRATOS_CLASS_POINTER_DEFINITION(MortarMapper);

    typedef TMappingMatrixBuilder TMappingMatrixBuilderType;

    typedef typename TMappingMatrixBuilderType::TDataType TDataType;

    typedef typename TMappingMatrixBuilderType::TSystemMatrixType TSystemMatrixType;

    typedef typename TMappingMatrixBuilderType::TSystemVectorType TSystemVectorType;

    typedef typename TMappingMatrixBuilderType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename TMappingMatrixBuilderType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename TMappingMatrixBuilderType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename TMappingMatrixBuilderType::LocalSystemVectorType LocalSystemVectorType;

    typedef std::size_t SizeType;

    typedef Node<3> NodeType;

    typedef Geometry<NodeType> GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    MortarMapper(ModelPart& i_model_part_origin, ModelPart& i_model_part_destination,
                          Parameters rJsonParameters) : MapperMatrixBased<TMappingMatrixBuilder, TLinearSolver>(
                                  i_model_part_origin, i_model_part_destination, rJsonParameters)
    {






        // @Jordi the "Condition_Gauss_Point" can be of the type that Riccardo suggested. Then we can get directly the ShapeFunctionValues that we need
        // mpMapperCommunicator->InitializeOrigin(MapperUtilities::Condition_Gauss_Point);
        // mpMapperCommunicator->InitializeDestination(MapperUtilities::Condition_Gauss_Point);
        // mpMapperCommunicator->Initialize();
        // these three steps correspond to "MapperCommunicator::BuildInterface()"

        // //mpMapperCommunicator->GetBuilderAndMultiplier()->BuildLHS(scheme, modelpart, Mdo); // could be moved to baseclass...?
        // FillMappingMatrix();
        // mpMapperCommunicator->GetBuilderAndSolver()->BuildLHS(scheme, modelpart, mM_dd); // this would also initialize this "BuilderAndSolver" => same as below, should be member of this class?
    }

    /// Destructor.
    virtual ~MortarMapper() { }

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
        return "MortarMapper";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MortarMapper";
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

    TSystemVectorPointerType mpQtmp; // for Mortar (needed bcs Trilinos cannot multiply in place)
    TSystemMatrixPointerType mpMdd; // for Mortar

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    void ExecuteMappingStep(const Kratos::Flags& MappingOptions) override // TODO move to private?
    {
        if (MappingOptions.Is(MapperFlags::CONSERVATIVE))
        {
            // mpBuilderAndSolver->SystemSolve(mpMdd, mpQtmp, this->mpQd); // Jordi the trilinos call also wants a modelpart!
            const bool transpose_flag = true;
            this->mpMappingMatrixBuilder->Multiply(*(this->mpMdo), *(this->mpQd), *(this->mpQo), transpose_flag);
        }
        else
        {
            this->mpMappingMatrixBuilder->Multiply(*(this->mpMdo), *(this->mpQo), *(this->mpQd));
            // mpBuilderAndSolver->SystemSolve(mpMdd, this->mpQd, mpQtmp); // Jordi the trilinos call also wants a modelpart!
        }

            //     if (InverseOperation) // for conservative mapping
    //     {
    //        if (mIsMortar)
    //         {
    //             mpBuilderAndSolver->SystemSolve(mpMdd, mpQtmp, mpQd); // Jordi the trilinos call also wants a modelpart!
    //             TSparseSpace::TransposeMult(mpMdo, mpQtmp, mpQo);
    //         }
    //         else
    //         {
    //             TSparseSpace::TransposeMult(mpMdo, mpQd, mpQo);
    //         }
    //     }
    //     else
    //     {
    //         if (mIsMortar)
    //         {
    //             TSparseSpace::Mult(mpMdo, mpQo, mpQtmp);
    //             mpBuilderAndSolver->SystemSolve(mpMdd, mpQd, mpQtmp); // Jordi the trilinos call also wants a modelpart!
    //             // if this turns out to be a problem we can pass the BuilderAndSolver to the MappingMatrixBuilder to call
    //             // the correct function, bcs this is different for serial and Trilinos
    //             // I would try to avoid duplicating the MapperStrategy at all cost!!!
    //         }
    //         else
    //         {
    //             TSparseSpace::Mult(mpMdo, mpQo, mpQd);
    //         }
    //     }
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

    // // @Jordi same question as for the base class
    // TSystemMatrixType mM_dd;
    // TSystemVectorType mQ_tmp;

    // strategy

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void InitializeInverseMapper() override
    {
        this->mpInverseMapper = Mapper::Pointer( new MortarMapper(this->mrModelPartDestination,
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

        Vector geom_information(4*num_points + 1);
        geom_information[0] = r_geom.GetGeometryFamily(); // Save the enum for the later reconstruction

        // "Serializing" the Geometry
        for (SizeType i=0; i<num_points; ++i)
        {
            NodeType& r_node = r_geom.GetPoint(i);

            geom_information[4*i+1] = r_node.GetValue(MAPPING_MATRIX_EQUATION_ID);
            geom_information[4*i+2] = r_node.X();
            geom_information[4*i+3] = r_node.Y();
            geom_information[4*i+4] = r_node.Z();
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
    // MortarMapper<TSparseSpace, TDenseSpace, TLinearSolver>& operator=(MortarMapper<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

    /// Copy constructor.
    //MortarMapper(MortarMapper const& rOther);

    ///@}

}; // Class MortarMapper

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


// inline std::istream & operator >>(std::istream& rIStream,
//                                   MortarMapper<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
// {
//     return rIStream;
// }

// /// output stream function

// inline std::ostream & operator <<(std::ostream& rOStream,
//                                   const MortarMapper<TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);

//     return rOStream;
// }
///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_MORTAR_MAPPER_H_INCLUDED  defined