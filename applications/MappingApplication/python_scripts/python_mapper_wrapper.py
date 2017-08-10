from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CreateMapper(model_part_origin, model_part_destination, custom_settings, mapper_index):

    if (type(model_part_origin) != KratosMultiphysics.ModelPart):
        raise Exception("input is expected to be provided as a Kratos ModelPart object")

    if (type(model_part_destination) != KratosMultiphysics.ModelPart):
        raise Exception("input is expected to be provided as a Kratos ModelPart object")
    
    if (type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("input is expected to be provided as a Kratos Parameters object")

    if (type(mapper_index) != int):
        raise Exception("input is expected to be provided as an integer")

    parallel_type = "OpenMP" # this is the default
    if custom_settings["problem_data"].Has("parallel_type"):
        parallel_type = custom_settings["problem_data"]["parallel_type"].GetString()

    # save the type of parallelsim in the mapper settings
    mapper_settings = custom_settings["mapper_settings"][mapper_index]
    mapper_settings["parallel_type"].SetString(parallel_type)

    mapper_type = mapper_settings["mapper_type"]

    # Mapper for OpenMP parallelism
    if (parallel_type == "OpenMP"):
        if (mapper_type == "NearestNeighbor"):
            mapper_module_name = "nearest_neighbor_mapper"

        elif (mapper_type == "NearestElement"):
            mapper_module_name = "nearest_element_mapper"

        elif (mapper_type == "NearestNeighborMatrixBased"):
            mapper_module_name = "nearest_neighbor_mapper_matrix_based"

        elif (mapper_type == "NearestElementMatrixBased"):
            mapper_module_name = "nearest_element_mapper_matrix_based"

        elif (mapper_type == "Mortar"):
            mapper_module_name = "mortar_mapper"

        else:
            raise Exception("the requested mapper type is not in the python mapper wrapper")

    # Mapper for MPI parallelism
    elif (parallel_type == "MPI"):
        if (mapper_type == "NearestNeighbor"):
            mapper_module_name = "nearest_neighbor_mapper_mpi"

        elif (mapper_type == "NearestElement"):
            mapper_module_name = "nearest_element_mapper_mpi"

        elif (mapper_type == "NearestNeighborMatrixBased"):
            mapper_module_name = "trilinos_nearest_neighbor_mapper_matrix_based"

        elif (mapper_type == "NearestElementMatrixBased"):
            mapper_module_name = "trilinos_nearest_element_mapper_matrix_based"

        elif (mapper_type == "Mortar"):
            mapper_module_name = "trilinos_mortar_mapper"

        else:
            raise Exception("the requested mapper type is not in the python mapper wrapper")
    else:
        raise Exception("parallel_type is neither OpenMP nor MPI")

    mapper_module = __import__(mapper_module_name)
    mapper = mapper_module.CreateMapper(model_part_origin, model_part_destination, mapper_settings)

    return mapper
