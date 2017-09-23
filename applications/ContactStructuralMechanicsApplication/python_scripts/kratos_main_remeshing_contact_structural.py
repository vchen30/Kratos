from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication  import *
from KratosMultiphysics.StructuralMechanicsApplication  import *
from KratosMultiphysics.ContactStructuralMechanicsApplication  import *
from KratosMultiphysics.MeshingApplication import *
## Import define_output
parameter_file = open("ProjectParameters.json",'r')
ProjectParameters = Parameters( parameter_file.read())

## Get echo level and parallel type
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
parallel_type = ProjectParameters["problem_data"]["parallel_type"].GetString()

## Import parallel modules if needed
if (parallel_type == "MPI"):
    from KratosMultiphysics.mpi import *
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.TrilinosApplication import *

## Structure model part definition
main_model_part = ModelPart(ProjectParameters["problem_data"]["model_part_name"].GetString())
main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, ProjectParameters["problem_data"]["domain_size"].GetInt())

## Solver construction
import python_solvers_wrapper_contact_structural
solver = python_solvers_wrapper_contact_structural.CreateSolver(main_model_part, ProjectParameters)

solver.AddVariables()

# Adding variables related to the remeshing
main_model_part.AddNodalSolutionStepVariable(NODAL_H) 
main_model_part.AddNodalSolutionStepVariable(NODAL_AREA)
main_model_part.AddNodalSolutionStepVariable(ERROR_ESTIMATE)

## Read the model - note that SetBufferSize is done here
solver.ImportModelPart()

## Add AddDofs
solver.AddDofs()

## Initialize GiD  I/O
output_post  = ProjectParameters.Has("output_configuration")
if (output_post == True):
    if (parallel_type == "OpenMP"):
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(solver.GetComputingModelPart(),
                                      ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                      ProjectParameters["output_configuration"])
    elif (parallel_type == "MPI"):
        from gid_output_process_mpi import GiDOutputProcessMPI
        gid_output = GiDOutputProcessMPI(solver.GetComputingModelPart(),
                                         ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                         ProjectParameters["output_configuration"])

    gid_output.ExecuteInitialize()

## Creation of the Kratos model (build sub_model_parts or submeshes)
StructureModel = {ProjectParameters["problem_data"]["model_part_name"].GetString(): main_model_part}

## Get the list of the sub_model_parts in where the processes are to be applied
for i in range(ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
    part_name = ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
    StructureModel.update({part_name: main_model_part.GetSubModelPart(part_name)})

## Print model_part and properties
if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 1):
    print("")
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

## Processes construction
import process_factory

## Remeshing processes construction
if (ProjectParameters.Has("initial_remeshing_process") == True):
    remeshing_processes = process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["initial_remeshing_process"])
    if (ProjectParameters.Has("list_other_processes") == True):
        remeshing_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["list_other_processes"])
    
    ## Remeshing processes initialization
    print("STARTING ADAPTATIVE LOOP")
    if (ProjectParameters.Has("adaptative_loop") == True):
        adaptative_loop = ProjectParameters["adaptative_loop"].GetInt()
    else:
        adaptative_loop = 1
    for n in range(adaptative_loop):
        print("ADAPTATIVE INTERATION: ", n + 1)
        for process in reversed(remeshing_processes):
            process.ExecuteInitialize()

list_of_processes = process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["constraints_process_list"])
list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["loads_process_list"])
if (ProjectParameters.Has("list_other_processes") == True):
    list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["list_other_processes"])
if (ProjectParameters.Has("contact_process_list") == True):
    list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["contact_process_list"])
if (ProjectParameters.Has("json_output_process") == True):
    list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["json_output_process"])
if (ProjectParameters.Has("recursive_remeshing_process") == True):
        list_of_processes += process_factory.KratosProcessFactory(StructureModel).ConstructListOfProcesses(ProjectParameters["recursive_remeshing_process"])

if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 1):
    for process in list_of_processes:
        print(process)

## Processes initialization
for process in list_of_processes:
    process.ExecuteInitialize()

## Add the processes to the solver
solver.AddProcessesList(list_of_processes)

## Solver initialization
solver.Initialize()
solver.SetEchoLevel(echo_level)

if (output_post == True):
    gid_output.ExecuteBeforeSolutionLoop()

for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()

## Writing the full ProjectParameters file before solving
if ((parallel_type == "OpenMP") or (mpi.rank == 0)) and (echo_level > 0):
    f = open("ProjectParametersOutput.json", 'w')
    f.write(ProjectParameters.PrettyPrintJsonString())
    f.close()

## Stepping and time settings
delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
start_time = ProjectParameters["problem_data"]["start_time"].GetDouble()
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

time = start_time
main_model_part.ProcessInfo[TIME_STEPS] = 0

# Solving the problem (time integration)
while(time <= end_time):

    time = time + delta_time
    main_model_part.ProcessInfo[TIME_STEPS] += 1
    main_model_part.CloneTimeStep(time)

    if (parallel_type == "OpenMP") or (mpi.rank == 0):
        print("")
        print("STEP = ", main_model_part.ProcessInfo[TIME_STEPS])
        print("TIME = ", time)

        if (main_model_part.Is(MODIFIED) == True):
            # WE INITIALIZE THE SOLVER
            solver.Initialize()
            # WE RECOMPUTE THE PROCESSES AGAIN
            ## Processes initialization
            for process in list_of_processes:
                process.ExecuteInitialize()
            ## Processes before the loop
            for process in list_of_processes:
                process.ExecuteBeforeSolutionLoop()
            ## Processes of initialize the solution step
            for process in list_of_processes:
                process.ExecuteInitializeSolutionStep()
        else:
            for process in list_of_processes:
                process.ExecuteInitializeSolutionStep()

    if (output_post == True):
        gid_output.ExecuteInitializeSolutionStep()

    solver.Solve()

    for process in list_of_processes:
        process.ExecuteFinalizeSolutionStep()

    if (output_post == True):
        gid_output.ExecuteFinalizeSolutionStep()

    for process in list_of_processes:
        process.ExecuteBeforeOutputStep()

    if (output_post == True) and (gid_output.IsOutputStep()):
        gid_output.PrintOutput()

    for process in list_of_processes:
        process.ExecuteAfterOutputStep()

for process in list_of_processes:
    process.ExecuteFinalize()

if (output_post == True):
    gid_output.ExecuteFinalize()

if (parallel_type == "OpenMP") or (mpi.rank == 0):
    print(" ")
    print("::[KSM Simulation]:: Analysis -END- ")
    print(" ")
