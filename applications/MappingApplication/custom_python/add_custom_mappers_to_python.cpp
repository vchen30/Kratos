// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
//         -        Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//         -        Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
//                  in the documentation and/or other materials provided with the distribution.
//         -        All advertising materials mentioning features or use of this software must display the following acknowledgement:
//                         This product includes Kratos Multi-Physics technology.
//         -        Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"

#include "spaces/ublas_space.h" // Always needed, for "LocalSpaceType"

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "trilinos_space.h"
#include "Epetra_FEVector.h"
#endif

#include "includes/kratos_parameters.h"
#include "custom_utilities/mapper_flags.h"

#include "custom_python/add_custom_mappers_to_python.h"

#include "custom_mappers/mapper.h"
#include "custom_mappers/nearest_neighbor_mapper.h"
#include "custom_mappers/nearest_element_mapper.h"

#include "custom_mappers/nearest_neighbor_mapper_matrix.h"
#include "custom_mappers/mortar_mapper.h"

#include "custom_utilities/mapper_factory_new.h"

#include "custom_strategies/builders/ublas_mapping_matrix_builder.h"
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_strategies/builders/trilinos_mapping_matrix_builder.h"
#endif


namespace Kratos
{

namespace Python
{

// Wrapper functions for taking a default argument for the flags
void UpdateInterface(Mapper& dummy)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    double dummy_search_radius = -1.0f;
    dummy.UpdateInterface(dummy_flags, dummy_search_radius);
}

void UpdateInterface(Mapper& dummy, Kratos::Flags options)
{
    double dummy_search_radius = -1.0f;
    dummy.UpdateInterface(options, dummy_search_radius);
}

void UpdateInterface(Mapper& dummy, double search_radius)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.UpdateInterface(dummy_flags, search_radius);
}

void Map(Mapper& dummy,
         const Variable<double>& origin_variable,
         const Variable<double>& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.Map(origin_variable, destination_variable, dummy_flags);
}

void Map(Mapper& dummy,
         const Variable< array_1d<double, 3> >& origin_variable,
         const Variable< array_1d<double, 3> >& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.Map(origin_variable, destination_variable, dummy_flags);
}

void InverseMap(Mapper& dummy,
                const Variable<double>& origin_variable,
                const Variable<double>& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.InverseMap(origin_variable, destination_variable, dummy_flags);
}

void InverseMap(Mapper& dummy,
                const Variable< array_1d<double, 3> >& origin_variable,
                const Variable< array_1d<double, 3> >& destination_variable)
{
    Kratos::Flags dummy_flags = Kratos::Flags();
    dummy.InverseMap(origin_variable, destination_variable, dummy_flags);
}

void  AddCustomMappersToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, CompressedMatrix, Vector> SerialSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SerialSparseSpaceType, LocalSpaceType> SerialLinearSolverType; // for Mortar
    typedef UblasMappingMatrixBuilder<SerialSparseSpaceType, LocalSpaceType> SerialMappingMatrixBuilderType;

// Overwrite the SparseSpaceType in case of mpi-parallel execution
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef LinearSolver<TrilinosSparseSpaceType, LocalSpaceType> TrilinosLinearSolverType; // for Mortar
    typedef TrilinosMappingMatrixBuilder<TrilinosSparseSpaceType, LocalSpaceType> TrilinosMappingMatrixBuilderType;
#endif

    void (*pUpdateInterface)(Mapper &)
        = &UpdateInterface;

    void (*pUpdateInterfaceOptions)(Mapper &, Kratos::Flags)
        = &UpdateInterface;

    void (*pUpdateInterfaceSearchRadius)(Mapper &, double)
        = &UpdateInterface;

    void (*pMapScalar)(Mapper &,
                       const Variable<double> &,
                       const Variable<double> &)
        = &Map;

    void (*pMapVector)(Mapper &,
                       const Variable< array_1d<double, 3> > &,
                       const Variable< array_1d<double, 3> > &)
        = &Map;

    void (*pInverseMapScalar)(Mapper &,
                              const Variable<double> &,
                              const Variable<double> &)
        = &InverseMap;

    void (*pInverseMapVector)(Mapper &,
                              const Variable< array_1d<double, 3> > &,
                              const Variable< array_1d<double, 3> > &)
        = &InverseMap;


    void (Mapper::*pUpdateInterfaceFull)(Kratos::Flags, double)
        = &Mapper::UpdateInterface;

    void (Mapper::*pMapScalarOptions)(const Variable<double> &,
            const Variable<double> &,
            Kratos::Flags)
        = &Mapper::Map;

    void (Mapper::*pMapVectorOptions)(const Variable< array_1d<double, 3> > &,
            const Variable< array_1d<double, 3> > &,
            Kratos::Flags)
        = &Mapper::Map;

    void (Mapper::*pInverseMapScalarOptions)(const Variable<double> &,
            const Variable<double> &,
            Kratos::Flags)
        = &Mapper::InverseMap;

    void (Mapper::*pInverseMapVectorOptions)(const Variable< array_1d<double, 3> > &,
            const Variable< array_1d<double, 3> > &,
            Kratos::Flags)
        = &Mapper::InverseMap;
    
    // Exposing the base class of the Mappers to Python, but without constructor
    class_< Mapper, boost::noncopyable > mapper 
        = class_< Mapper, boost::noncopyable >("Mapper", no_init)
            .def("UpdateInterface",  pUpdateInterface)
            .def("UpdateInterface",  pUpdateInterfaceOptions)
            .def("UpdateInterface",  pUpdateInterfaceSearchRadius)
            .def("Map",              pMapScalar)
            .def("Map",              pMapVector)
            .def("InverseMap",       pInverseMapScalar)
            .def("InverseMap",       pInverseMapVector)

            .def("UpdateInterface",  pUpdateInterfaceFull)
            .def("Map",              pMapScalarOptions)
            .def("Map",              pMapVectorOptions)
            .def("InverseMap",       pInverseMapScalarOptions)
            .def("InverseMap",       pInverseMapVectorOptions)
            ;
    
    // Adding the flags that can be used while mapping
    mapper.attr("SWAP_SIGN") = MapperFlags::SWAP_SIGN;
    mapper.attr("ADD_VALUES") = MapperFlags::ADD_VALUES;
    mapper.attr("CONSERVATIVE") = MapperFlags::CONSERVATIVE;
    mapper.attr("REMESHED") = MapperFlags::REMESHED;

    // Jordi is it possible to expose the mappers without a constructor and use them only through the factory?
    // This would circumvent problems with the wrong space being selected

    // Exposing the Mappers
        // Matrix-free Mappers
        class_< NearestNeighborMapper, bases<Mapper>, boost::noncopyable>
        ("NearestNeighborMapper", init<ModelPart&, ModelPart&, Parameters>());
        class_< NearestElementMapper, bases<Mapper>, boost::noncopyable>
        ("NearestElementMapper", init<ModelPart&, ModelPart&, Parameters>());
        // Matrix-based Mappers
        class_<NearestNeighborMapperMatrix<SerialMappingMatrixBuilderType, SerialLinearSolverType>,
               bases<Mapper>, boost::noncopyable>("NearestNeighborMapperMatrix", init<ModelPart &, ModelPart &, Parameters>());
        class_<MortarMapper<SerialMappingMatrixBuilderType, SerialLinearSolverType>,
               bases<Mapper>, boost::noncopyable>("MortarMapper", init<ModelPart &, ModelPart &, Parameters>());
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
        class_<NearestNeighborMapperMatrix<TrilinosMappingMatrixBuilderType, TrilinosLinearSolverType>,
            bases<Mapper>, boost::noncopyable>("NearestNeighborMapperMatrix", init<ModelPart &, ModelPart &, Parameters>());
        class_<MortarMapper<TrilinosMappingMatrixBuilderType, TrilinosLinearSolverType>,
            bases<Mapper>, boost::noncopyable>("MortarMapper", init<ModelPart &, ModelPart &, Parameters>());
#endif

    
    // Exposing the MapperFactory
    class_< MapperFactoryNew, boost::noncopyable>("MapperFactoryNew", no_init)
    .def("CreateMapper", &MapperFactoryNew::CreateMapper, return_value_policy<manage_new_object>())
    .staticmethod("CreateMapper")
    ;
    /*
    Jordi
    This should work according to "https://wiki.python.org/moin/boost.python/HowTo", search for "manage_new_object"
    see also "http://www.boost.org/doc/libs/1_61_0/libs/python/doc/html/tutorial/tutorial/exposing.html", search for "ownership"
    and "http://www.boost.org/doc/libs/1_61_0/libs/python/doc/html/tutorial/tutorial/functions.html#tutorial.functions.call_policies", search for "manage_new_object"
    UPDATE it works, but only if the MapperFactory returns a raw pointer to the mapper (i.e. Mapper::Pointer as return type does not work)
    */
    

}

}  // namespace Python.

} // Namespace Kratos
