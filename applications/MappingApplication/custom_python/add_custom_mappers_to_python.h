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

#if !defined(KRATOS_ADD_MAPPERS_TO_PYTHON_H_INCLUDED )
#define  KRATOS_ADD_MAPPERS_TO_PYTHON_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "spaces/ublas_space.h" // Always needed, for "LocalSpaceType"
#include "linear_solvers/linear_solver.h"
#include "includes/kratos_parameters.h"

// Mapper Aux stuff
#include "custom_utilities/mapper_flags.h"
#include "custom_strategies/builders/ublas_mapping_matrix_builder.h"

// Mapper base class
#include "custom_mappers/mapper.h"

// Matrix-free Mappers
#include "custom_mappers/nearest_neighbor_mapper.h"
#include "custom_mappers/nearest_element_mapper.h"

// Matrix-based Mappers
#include "custom_mappers/nearest_neighbor_mapper_matrix.h"
#include "custom_mappers/nearest_element_mapper_matrix.h"
#include "custom_mappers/mortar_mapper.h"


namespace Kratos
{   

namespace Python
{

void AddCustomMappersToPython();
void AddCustomMappersMPIToPython();


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_ADD_MAPPERS_TO_PYTHON_H_INCLUDED  defined
