//
//   Project Name:        Kratos
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2008-10-13 08:56:42 $
//   Revision:            $Revision: 1.5 $
//
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

// #if !defined(EXECUTE_REMESHING )
// #define  EXECUTE_REMESHING

// /* External includes */
// #include "boost/smart_ptr.hpp"

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

// #include "meshing_application.h"
#ifdef INCLUDE_MMG
    #include "custom_processes/mmg_process.h"
#endif


namespace Kratos
{
    #if !defined(REMESHING_UTILITIES)
    #define REMESHING_UTILITIES
        enum RemeshingUtilities {MMG = 0};
    #endif

template< unsigned int TDim >
class ExecuteRemeshing
{
public:

ExecuteRemeshing(ModelPart& mp, Parameters params = Parameters(R"({})")): mThisModelPart (mp), mThisParameters(params)
{
}

void Execute()
{
#ifdef INCLUDE_MMG
    auto MmgRemesh = MmgProcess<TDim>(mThisModelPart, mThisParameters); 
    MmgRemesh.Execute();
#endif
}
private:
ModelPart mThisModelPart;
Parameters mThisParameters;
};



}  // namespace Kratos.

// #endif // KRATOS_PROJECTION  defined 


