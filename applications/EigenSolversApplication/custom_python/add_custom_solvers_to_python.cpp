/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Thomas Oberbichler
*/

// System includes
#include <complex>

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_solvers_to_python.h"

#include "includes/ublas_complex_interface.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_solvers/eigen_direct_solver.h"

namespace Kratos
{

namespace Python
{

void AddCustomSolversToPython()
{
	using namespace boost::python;

	typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
	typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
	typedef LinearSolver<SparseSpaceType, LocalSpaceType> LinearSolverType;
	typedef DirectSolver<SparseSpaceType, LocalSpaceType> DirectSolverType;

	typedef std::complex<double> Complex;
    typedef UblasSpace<Complex, ComplexCompressedMatrix, ComplexVector> ComplexSparseSpaceType;
    typedef UblasSpace<Complex, ComplexMatrix, ComplexVector> ComplexLocalSpaceType;
    typedef LinearSolver<ComplexSparseSpaceType, ComplexLocalSpaceType> ComplexLinearSolverType;
	typedef DirectSolver<ComplexSparseSpaceType, ComplexLocalSpaceType> ComplexDirectSolverType;

	using SparseLUSolver = EigenDirectSolver<SparseLU<double>, SparseSpaceType, LocalSpaceType>;
	class_<SparseLUSolver, bases<DirectSolverType>, boost::noncopyable>
		("SparseLUSolver", init<>())
		.def(init<Parameters>());

	using ComplexSparseLUSolver = EigenDirectSolver<SparseLU<Complex>, ComplexSparseSpaceType, ComplexLocalSpaceType>;
	class_<ComplexSparseLUSolver, bases<ComplexDirectSolverType>, boost::noncopyable>
		("ComplexSparseLUSolver", init<>())
		.def(init<Parameters>());
	
	#if defined EIGEN_USE_MKL_ALL
	using PardisoLLTSolver = EigenDirectSolver<PardisoLLT<double>, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLLTSolver, bases<DirectSolverType>, boost::noncopyable>
		("PardisoLLTSolver", init<>())
		.def(init<Parameters>());

	using PardisoLDLTSolver = EigenDirectSolver<PardisoLDLT<double>, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLDLTSolver, bases<DirectSolverType>, boost::noncopyable>
		("PardisoLDLTSolver", init<>())
		.def(init<Parameters>());

	using PardisoLUSolver = EigenDirectSolver<PardisoLU<double>, SparseSpaceType, LocalSpaceType>;
	class_<PardisoLUSolver, bases<DirectSolverType>, boost::noncopyable>
		("PardisoLUSolver", init<>())
		.def(init<Parameters>());
	#endif
}

} // namespace Python

} // namespace Kratos
