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
