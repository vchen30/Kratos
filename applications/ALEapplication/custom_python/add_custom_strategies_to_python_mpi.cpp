//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Andreas Winterstein (a.winterstein@tum.de)
//

// System includes

// External includes
#include "spaces/ublas_space.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

//Trilinos includes
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"

// Project includes
#include "trilinos_space.h"
#include "custom_python/add_custom_strategies_to_python.h"

// strategies
#include "custom_strategies/strategies/trilinos_laplacian_meshmoving_strategy.h"
#include "custom_strategies/strategies/trilinos_structural_meshmoving_strategy.h"

// linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos {

namespace Python {
using namespace boost::python;

void AddCustomStrategiesMPIToPython() {


    // Mesh Moving ********************************************************************************************
    typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
    typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;
    typedef SolvingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > TrilinosBaseSolvingStrategyType;

    using TrilinosLaplacianMeshMovingStrategyType = TrilinosLaplacianMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >;
    class_< TrilinosLaplacianMeshMovingStrategyType, bases<TrilinosBaseSolvingStrategyType>, boost::noncopyable >
    ("TrilinosLaplacianMeshMovingStrategy", init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, int, bool, bool, int >() )
    ;

    using TrilinosStructuralMeshMovingStrategyType = TrilinosStructuralMeshMovingStrategy< TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType >;
    class_< TrilinosStructuralMeshMovingStrategyType, bases< TrilinosBaseSolvingStrategyType >, boost::noncopyable >
    ("TrilinosStructuralMeshMovingStrategy", init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, int, bool, bool, int >() )
    ;
}

} // namespace Python.

} // Namespace Kratos