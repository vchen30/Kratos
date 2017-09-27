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

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "includes/kratos_parameters.h"
#include "custom_utilities/mapper_flags.h"
// #include "custom_utilities/mapper_factory.h"

#include "custom_utilities/mapper_communicator.h"


#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_utilities/mapper_mpi_communicator.h"
#endif


namespace Kratos
{

namespace Python
{

// // Wrapper functions for taking a default argument for the flags
// void UpdateInterface(MapperFactory& dummy)
// {
//     Kratos::Flags dummy_flags = Kratos::Flags();
//     double dummy_search_radius = -1.0f;
//     dummy.UpdateInterface(dummy_flags, dummy_search_radius);
// }

// void UpdateInterface(MapperFactory& dummy, Kratos::Flags& options)
// {
//     double dummy_search_radius = -1.0f;
//     dummy.UpdateInterface(options, dummy_search_radius);
// }

// void UpdateInterface(MapperFactory& dummy, double search_radius)
// {
//     Kratos::Flags dummy_flags = Kratos::Flags();
//     dummy.UpdateInterface(dummy_flags, search_radius);
// }


// void Map(MapperFactory& dummy,
//          const Variable<double>& origin_variable,
//          const Variable<double>& destination_variable)
// {
//     Kratos::Flags dummy_flags = Kratos::Flags();
//     dummy.Map(origin_variable, destination_variable, dummy_flags);
// }

// void Map(MapperFactory& dummy,
//          const Variable< array_1d<double, 3> >& origin_variable,
//          const Variable< array_1d<double, 3> >& destination_variable)
// {
//     Kratos::Flags dummy_flags = Kratos::Flags();
//     dummy.Map(origin_variable, destination_variable, dummy_flags);
// }

// void InverseMap(MapperFactory& dummy,
//                 const Variable<double>& origin_variable,
//                 const Variable<double>& destination_variable)
// {
//     Kratos::Flags dummy_flags = Kratos::Flags();
//     dummy.InverseMap(origin_variable, destination_variable, dummy_flags);
// }

// void InverseMap(MapperFactory& dummy,
//                 const Variable< array_1d<double, 3> >& origin_variable,
//                 const Variable< array_1d<double, 3> >& destination_variable)
// {
//     Kratos::Flags dummy_flags = Kratos::Flags();
//     dummy.InverseMap(origin_variable, destination_variable, dummy_flags);
// }

void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    // void (*pUpdateInterface)(MapperFactory &)
    //     = &UpdateInterface;

    // void (*pUpdateInterfaceOptions)(MapperFactory &, Kratos::Flags &)
    //     = &UpdateInterface;

    // void (*pUpdateInterfaceSearchRadius)(MapperFactory &, double)
    //     = &UpdateInterface;

    // void (*pMapScalar)(MapperFactory &,
    //                    const Variable<double> &,
    //                    const Variable<double> &)
    //     = &Map;

    // void (*pMapVector)(MapperFactory &,
    //                    const Variable< array_1d<double, 3> > &,
    //                    const Variable< array_1d<double, 3> > &)
    //     = &Map;

    // void (*pInverseMapScalar)(MapperFactory &,
    //                           const Variable<double> &,
    //                           const Variable<double> &)
    //     = &InverseMap;

    // void (*pInverseMapVector)(MapperFactory &,
    //                           const Variable< array_1d<double, 3> > &,
    //                           const Variable< array_1d<double, 3> > &)
    //     = &InverseMap;


    // void (MapperFactory::*pUpdateInterfaceFull)(Kratos::Flags &, double)
    //     = &MapperFactory::UpdateInterface;

    // void (MapperFactory::*pMapScalarOptions)(const Variable<double> &,
    //         const Variable<double> &,
    //         Kratos::Flags &)
    //     = &MapperFactory::Map;

    // void (MapperFactory::*pMapVectorOptions)(const Variable< array_1d<double, 3> > &,
    //         const Variable< array_1d<double, 3> > &,
    //         Kratos::Flags &)
    //     = &MapperFactory::Map;

    // void (MapperFactory::*pInverseMapScalarOptions)(const Variable<double> &,
    //         const Variable<double> &,
    //         Kratos::Flags &)
    //     = &MapperFactory::InverseMap;

    // void (MapperFactory::*pInverseMapVectorOptions)(const Variable< array_1d<double, 3> > &,
    //         const Variable< array_1d<double, 3> > &,
    //         Kratos::Flags &)
    //     = &MapperFactory::InverseMap;


    // class_< MapperFactory > mapper_factory = class_<MapperFactory>("MapperFactory", init<ModelPart&, ModelPart&, Parameters>())
    //         .def("UpdateInterface",  pUpdateInterface)
    //         .def("UpdateInterface",  pUpdateInterfaceOptions)
    //         .def("UpdateInterface",  pUpdateInterfaceSearchRadius)
    //         .def("Map",              pMapScalar)
    //         .def("Map",              pMapVector)
    //         .def("InverseMap",       pInverseMapScalar)
    //         .def("InverseMap",       pInverseMapVector)

    //         .def("UpdateInterface",  pUpdateInterfaceFull)
    //         .def("Map",              pMapScalarOptions)
    //         .def("Map",              pMapVectorOptions)
    //         .def("InverseMap",       pInverseMapScalarOptions)
    //         .def("InverseMap",       pInverseMapVectorOptions)
    //         ;

    // mapper_factory.attr("SWAP_SIGN") = MapperFlags::SWAP_SIGN;
    // mapper_factory.attr("ADD_VALUES") = MapperFlags::ADD_VALUES;
    // mapper_factory.attr("CONSERVATIVE") = MapperFlags::CONSERVATIVE;
    // mapper_factory.attr("REMESHED") = MapperFlags::REMESHED;

// stuff for the new MappingApplication design



    class_< MapperCommunicator, boost::noncopyable >
        ("MapperCommunicator", init<ModelPart&, ModelPart&, Parameters>())
        // TODO add the functions
        ;


#ifdef KRATOS_USING_MPI // mpi-parallel compilation
    class_< MapperMPICommunicator, bases< MapperCommunicator >,  boost::noncopyable >
        ("MapperMPICommunicator", init<ModelPart&, ModelPart&, Parameters>())
        ;
#endif
}

}  // namespace Python.

} // Namespace Kratos
