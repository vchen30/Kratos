set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMappingApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
# include_directories( ${CMAKE_SOURCE_DIR}/applications/mpi_search_application ) # TODO can be removed for now


# generate variables with the sources
set( KRATOS_MAPPING_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/mapping_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mapping_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/mapping_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_mappers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mapper_flags.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/base_mapper_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/nearest_neighbor_mapper_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/nearest_element_mapper_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mortar_mapper_condition.cpp
)

set(KRATOS_MAPPING_APPLICATION_EXTRA_LIBS "")
## If we have mpi support, use it
if(${MPI_NEEDED} MATCHES ON )
    # This is only true if Metis and Trilinos are compiled, right?
    # Jordi is this the correct way of including trilinos?

    # With this line I include the trilinos app. Or do I even need that? 
    # Since I don't use features of the TrilinosApplication but only the matrix, which comes directly from Trilinos
    include_directories( ${CMAKE_SOURCE_DIR}/applications/trilinos_application ) 

    # These three lines I copied from the CMakeLists.txt of the TrilinosApplication
    # Do we need all of them?
    find_package(TRILINOS REQUIRED) # move to core CMakeLists
    include_directories(${TRILINOS_INCLUDE_DIR}) # seems to work without this, too
    #link_directories(${TRILINOS_LIBRARIES})

    # This line was there orignally:
	set(KRATOS_MAPPING_APPLICATION_EXTRA_LIBS ${KRATOS_MAPPING_APPLICATION_EXTRA_LIBS} ${MPI_LIBRARIES}) # also is this needed? I seem to be the only one doing it...

    # How abt this line, do be need to link explicitly?
	set(KRATOS_MAPPING_APPLICATION_EXTRA_LIBS ${KRATOS_MAPPING_APPLICATION_EXTRA_LIBS} ${TRILINOS_LIBRARIES}) # Jordi why does it work without linking?
endif(${MPI_NEEDED} MATCHES ON )

# define library Kratos which defines the basic python interface
add_library(KratosMappingApplication SHARED ${KRATOS_MAPPING_APPLICATION_SOURCES})
target_link_libraries(KratosMappingApplication KratosCore ${KRATOS_MAPPING_APPLICATION_EXTRA_LIBS})
set_target_properties(KratosMappingApplication PROPERTIES PREFIX "")
install(TARGETS KratosMappingApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMappingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMappingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MappingApplication.py" DESTINATION KratosMultiphysics )
