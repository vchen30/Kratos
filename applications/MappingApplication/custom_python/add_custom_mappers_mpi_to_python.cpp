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
#include "custom_python/add_custom_mappers_to_python.h"
#include "trilinos_space.h"
#include "Epetra_FEVector.h"
#include "custom_strategies/builders/trilinos_mapping_matrix_builder.h"


namespace Kratos
{

namespace Python
{

void  AddCustomMappersMPIToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;

    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef LinearSolver<TrilinosSparseSpaceType, LocalSpaceType> TrilinosLinearSolverType; // for Mortar
    typedef TrilinosMappingMatrixBuilder<TrilinosSparseSpaceType, LocalSpaceType> TrilinosMappingMatrixBuilderType;

    
    class_<NearestNeighborMapperMatrix<TrilinosMappingMatrixBuilderType, TrilinosLinearSolverType>,
        bases<Mapper>, boost::noncopyable>("NearestNeighborMapperMatrix", init<ModelPart &, ModelPart &, Parameters>());
    class_<NearestElementMapperMatrix<TrilinosMappingMatrixBuilderType, TrilinosLinearSolverType>,
        bases<Mapper>, boost::noncopyable>("NearestElementMapperMatrix", init<ModelPart &, ModelPart &, Parameters>());
    class_<MortarMapper<TrilinosMappingMatrixBuilderType, TrilinosLinearSolverType>,
        bases<Mapper>, boost::noncopyable>("MortarMapper", init<ModelPart &, ModelPart &, Parameters>());

}

}  // namespace Python.

} // Namespace Kratos
