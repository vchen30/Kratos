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
#include "custom_utilities/interface_preprocess.h"
#include "custom_conditions/base_mapper_condition.h"
// #include "custom_strategies/builders/mapping_matrix_builder.h"
// #ifdef KRATOS_USING_MPI // mpi-parallel compilation
// #include "custom_strategies/builders/trilinos_mapping_matrix_builder.h"
// #endif

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
class MapperMatrixBased : public Mapper
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of MapperMatrixBased
    KRATOS_CLASS_POINTER_DEFINITION(MapperMatrixBased);

    typedef TMappingMatrixBuilder TMappingMatrixBuilderType;

    typedef typename TMappingMatrixBuilderType::Pointer TMappingMatrixBuilderPointerType;

    // Jordi can I use the same names without problems?
    typedef typename TMappingMatrixBuilderType::TDataType TDataType;

    typedef typename TMappingMatrixBuilderType::TSystemMatrixType TSystemMatrixType;

    typedef typename TMappingMatrixBuilderType::TSystemVectorType TSystemVectorType;

    typedef typename TMappingMatrixBuilderType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename TMappingMatrixBuilderType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename TMappingMatrixBuilderType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename TMappingMatrixBuilderType::LocalSystemVectorType LocalSystemVectorType;

    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > VectorComponentType;

    ///@}
    ///@name Life Cycle
    ///@{

    MapperMatrixBased(ModelPart &rModelPartOrigin, ModelPart &rModelpartDestination,
                      Parameters rJsonParameters) : Mapper(rModelPartOrigin, rModelpartDestination, rJsonParameters)
    {
        mpMappingMatrixBuilder = TMappingMatrixBuilderPointerType(new TMappingMatrixBuilderType(this->mEchoLevel));

        mpInterfaceModelPart = ModelPart::Pointer( new ModelPart("Mapper-Interface") );
        this->mpMapperCommunicator->pSetModelpartDestination(&*(this->mpInterfaceModelPart));        

        mpInterfacePreprocessor = InterfacePreprocess::Pointer( new InterfacePreprocess(this->mrModelPartDestination,
                                                                                        this->mpInterfaceModelPart) );

        InitializeMappingMatrix();
    }

    /// Destructor.
    virtual ~MapperMatrixBased() { }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /* This function maps from Origin to Destination */
    void Map(const Variable<double>& rOriginVariable,
             const Variable<double>& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::CONSERVATIVE))
        {
            MappingOptions.Reset(MapperFlags::CONSERVATIVE); // TODO test this!!!
            MappingOptions.Set(MapperFlags::USE_TRANSPOSE); // TODO test this!!!

            InverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else
        {
            // TODO move these three steps to a function
            InitializeMappingStep<>(rOriginVariable, 
                                    rDestinationVariable, 
                                    MappingOptions);
            
            ExecuteMappingStep(MappingOptions);

            FinalizeMappingStep<>(rOriginVariable, 
                                rDestinationVariable, 
                                MappingOptions);
        }
    }

    /* This function maps from Origin to Destination */
    void Map(const Variable< array_1d<double, 3> >& rOriginVariable,
             const Variable< array_1d<double, 3> >& rDestinationVariable,
             Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::CONSERVATIVE))
        {
            MappingOptions.Reset(MapperFlags::CONSERVATIVE); // TODO test this!!!
            MappingOptions.Set(MapperFlags::USE_TRANSPOSE); // TODO test this!!!

            InverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else
        {
            VectorComponentType var_component_x_origin = KratosComponents< VectorComponentType >::Get(rOriginVariable.Name()+std::string("_X"));
            VectorComponentType var_component_y_origin = KratosComponents< VectorComponentType >::Get(rOriginVariable.Name()+std::string("_Y"));
            VectorComponentType var_component_z_origin = KratosComponents< VectorComponentType >::Get(rOriginVariable.Name()+std::string("_Z"));

            VectorComponentType var_component_x_destination = KratosComponents< VectorComponentType >::Get(rDestinationVariable.Name()+std::string("_X"));
            VectorComponentType var_component_y_destination = KratosComponents< VectorComponentType >::Get(rDestinationVariable.Name()+std::string("_Y"));
            VectorComponentType var_component_z_destination = KratosComponents< VectorComponentType >::Get(rDestinationVariable.Name()+std::string("_Z"));

            // X-Component
            InitializeMappingStep<VectorComponentType>(var_component_x_origin, 
                                                    var_component_x_destination,
                                                    MappingOptions);

            ExecuteMappingStep(MappingOptions);

            FinalizeMappingStep<VectorComponentType>(var_component_x_origin, 
                                                    var_component_x_destination,  
                                                    MappingOptions);
            
            // Y-Component
            InitializeMappingStep<VectorComponentType>(var_component_y_origin, 
                                                    var_component_y_destination,
                                                    MappingOptions);

            ExecuteMappingStep(MappingOptions);
            
            FinalizeMappingStep<VectorComponentType>(var_component_y_origin, 
                                                    var_component_y_destination,  
                                                    MappingOptions);

            // Z-Component
            InitializeMappingStep<VectorComponentType>(var_component_z_origin, 
                                                    var_component_z_destination,
                                                    MappingOptions);

            ExecuteMappingStep(MappingOptions);
            
            FinalizeMappingStep<VectorComponentType>(var_component_z_origin, 
                                                    var_component_z_destination,  
                                                    MappingOptions);
        }
    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable<double>& rOriginVariable,
                    const Variable<double>& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::CONSERVATIVE))
        {
            MappingOptions.Reset(MapperFlags::CONSERVATIVE); // TODO test this!!!
            MappingOptions.Set(MapperFlags::USE_TRANSPOSE); // TODO test this!!!

            Map(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else
        {
            // Construct the inverse mapper if it hasn't been done before
            if (!mpInverseMapper) InitializeInverseMapper();
            
            mpInverseMapper->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }
    }

    /* This function maps from Destination to Origin */
    void InverseMap(const Variable< array_1d<double, 3> >& rOriginVariable,
                    const Variable< array_1d<double, 3> >& rDestinationVariable,
                    Kratos::Flags MappingOptions) override
    {
        if (MappingOptions.Is(MapperFlags::CONSERVATIVE))
        {
            MappingOptions.Reset(MapperFlags::CONSERVATIVE); // TODO test this!!! => seems to work
            MappingOptions.Set(MapperFlags::USE_TRANSPOSE); // TODO test this!!!

            Map(rOriginVariable, rDestinationVariable, MappingOptions);
        }
        else
        {
            // Construct the inverse mapper if it hasn't been done before
            if (!mpInverseMapper) InitializeInverseMapper();
            
            mpInverseMapper->Map(rDestinationVariable, rOriginVariable, MappingOptions);
        }        
    }

    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mMdo = *mpMdo;

        return mMdo;
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

    TMappingMatrixBuilderPointerType mpMappingMatrixBuilder;

    TSystemVectorPointerType mpQo;
    TSystemVectorPointerType mpQd;
    TSystemMatrixPointerType mpMdo;

    ModelPart::Pointer mpInterfaceModelPart;
    InterfacePreprocess::Pointer mpInterfacePreprocessor;
    Parameters mInterfaceParameters = Parameters(R"({})");

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void UpdateInterfaceSpecific(Kratos::Flags MappingOptions) override {
        // TODO do the same for the InverseMapper(?)
        if (MappingOptions.Is(MapperFlags::REMESHED))
        {
            GenerateInterfaceModelPart();
            
            mpMappingMatrixBuilder->ClearData(mpMdo);
            mpMappingMatrixBuilder->ClearData(mpQo);
            mpMappingMatrixBuilder->ClearData(mpQd);

            InitializeMappingMatrix();
            ComputeMappingMatrix();
        }
        else
        {
            mpMappingMatrixBuilder->ClearData(mpMdo);
            mpMappingMatrixBuilder->ClearData(mpQo);
            mpMappingMatrixBuilder->ClearData(mpQd);

            InitializeMappingMatrix(false);
            ComputeMappingMatrix();
        }
    }

    void InitializeMappingMatrix(const bool UpdateSystem = true)
    {
        if (UpdateSystem)
        {
            mpMappingMatrixBuilder->SetUpSystem(mrModelPartOrigin);
            mpMappingMatrixBuilder->SetUpSystem(mrModelPartDestination);
        }

        const unsigned int num_local_nodes_origin = mrModelPartOrigin.GetCommunicator().LocalMesh().NumberOfNodes();
        const unsigned int num_local_nodes_destination = mrModelPartDestination.GetCommunicator().LocalMesh().NumberOfNodes();

        mpMappingMatrixBuilder->ResizeAndInitializeVectors(mpMdo, mpQo, mpQd, num_local_nodes_origin, num_local_nodes_destination);
    }

    void ComputeMappingMatrix()
    {
        this->ExchangeInterfaceGeometryData();
        
        mpMappingMatrixBuilder->BuildMappingMatrix(*mpInterfaceModelPart, *mpMdo);
    }

    /**
    This function creates the Interface-ModelPart
    */
    void GenerateInterfaceModelPart()
    {   
        this->mpInterfacePreprocessor->GenerateInterfacePart(mInterfaceParameters);
    }

    template< class TVarType>
    void InitializeMappingStep(const TVarType& rVarOrigin,
                               const TVarType& rVarDestination,
                               const Kratos::Flags& MappingOptions)
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE)) // TODO replace this with USE_TRANSPOSE!
        {
            mpMappingMatrixBuilder->UpdateSystemVector(mrModelPartDestination, *mpQd, rVarDestination);
        }
        else
        {
            mpMappingMatrixBuilder->UpdateSystemVector(mrModelPartOrigin, *mpQo, rVarOrigin);
        }
    }

    virtual void ExecuteMappingStep(const Kratos::Flags& MappingOptions) // Override this class in Mortar
    {
        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE))
        {
            const bool transpose_flag = true; // constexpr?
            mpMappingMatrixBuilder->Multiply(*mpMdo, *mpQd, *mpQo, transpose_flag);
        }
        else
        {
            mpMappingMatrixBuilder->Multiply(*mpMdo, *mpQo, *mpQd);
        }
    }

    template< class TVarType>
    void FinalizeMappingStep(const TVarType& rVarOrigin,
                             const TVarType& rVarDestination,
                             const Kratos::Flags& MappingOptions)
    {
        double factor = 1.0f;
        ProcessMappingOptions(MappingOptions, factor);

        if (MappingOptions.Is(MapperFlags::USE_TRANSPOSE))
        {

            mpMappingMatrixBuilder->Update(mrModelPartOrigin, *mpQo, rVarOrigin, MappingOptions, factor);
        }
        else
        {
            mpMappingMatrixBuilder->Update(mrModelPartDestination, *mpQd, rVarDestination, MappingOptions, factor);
        }
    }

    virtual void ExchangeInterfaceGeometryData() = 0; // TODO make private and only virtual

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