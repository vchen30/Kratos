// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "spaces/ublas_space.h"

// Strategies
#include "custom_strategies/custom_strategies/residualbased_newton_raphson_contact_strategy.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Schemes
#include "solving_strategies/schemes/scheme.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "custom_strategies/custom_convergencecriterias/error_mesh_criteria.h"
#include "custom_Strategies/custom_strategies/residualbased_newton_raphson_contact_strategy_remeshing.h"

// Builders and solvers
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"

// Linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{
using namespace boost::python;

void  AddCustomStrategiesToPython()
{
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
    typedef boost::shared_ptr<ProcessFactoryUtility> ProcessesListType;

    // Base types
    typedef ConvergenceCriteria< SparseSpaceType, LocalSpaceType > ConvergenceCriteriaType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;
    typedef ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointer;
    typedef BuilderAndSolver< SparseSpaceType, LocalSpaceType, LinearSolverType > BuilderAndSolverType;
        
    // Custom strategy types
    typedef ResidualBasedNewtonRaphsonContactStrategyRemeshing< SparseSpaceType, LocalSpaceType , LinearSolverType >  ResidualBasedNewtonRaphsonContactStrategyRemeshingType;
    typedef ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType , LinearSolverType > BaseSolvingStrategyType;
    // Custom scheme types

    // Custom convergence criterion types
    typedef ErrorMeshCriteria< SparseSpaceType,  LocalSpaceType > ErrorMeshCriteriaType;
    
    // Custom builder and solvers types
    
    //********************************************************************
    //*************************STRATEGY CLASSES***************************
    //********************************************************************

    // Residual Based Newton Raphson Contact Strategy      
    class_< ResidualBasedNewtonRaphsonContactStrategyRemeshingType, bases< BaseSolvingStrategyType >, boost::noncopyable >
            ("ResidualBasedNewtonRaphsonContactStrategyRemeshing", init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters, Parameters >())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, unsigned int, bool, bool, bool, Parameters, Parameters, ProcessesListType>())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters, Parameters >())
            .def(init < ModelPart&, BaseSchemeType::Pointer, LinearSolverType::Pointer, ConvergenceCriteriaType::Pointer, BuilderAndSolverType::Pointer, unsigned int, bool, bool, bool, Parameters, Parameters, ProcessesListType>())
            //.def("SetMaxIterationNumber", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetMaxIterationNumber)
            //.def("GetMaxIterationNumber", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetMaxIterationNumber)
            //.def("SetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::SetKeepSystemConstantDuringIterations)
            //.def("GetKeepSystemConstantDuringIterations", &ResidualBasedNewtonRaphsonContactStrategy< SparseSpaceType, LocalSpaceType, LinearSolverType >::GetKeepSystemConstantDuringIterations)
            ;
             
    //********************************************************************
    //*************************SCHEME CLASSES*****************************
    //********************************************************************

    //********************************************************************
    //*******************CONVERGENCE CRITERIA CLASSES*********************
    //********************************************************************

    // Displacement Convergence Criterion
    class_< ErrorMeshCriteriaType,
            bases< ConvergenceCriteriaType >, boost::noncopyable >
            (
            "ErrorMeshCriteria", 
            init<ModelPart&, Parameters>())
            .def("SetEchoLevel", &ErrorMeshCriteriaType::SetEchoLevel)
            ;
            
    //********************************************************************
    //*************************BUILDER AND SOLVER*************************
    //********************************************************************

}

}  // namespace Python.

} // Namespace Kratos

