set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosALEApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_ALE_APPLICATION_SOURCES
	      ${CMAKE_CURRENT_SOURCE_DIR}/ale_application.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/laplacian_meshmoving_element.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/structural_meshmoving_element.cpp
	      ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	      ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	      ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/ale_python_application.cpp
		  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/move_mesh_utilities.cpp
)

set(KRATOS_ALE_APPLICATION_EXTRA_LIBS "")

if(${MPI_NEEDED} MATCHES ON )
    set( KRATOS_ALE_APPLICATION_SOURCES ${KRATOS_ALE_APPLICATION_SOURCES}
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python_mpi.cpp
    )

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
	#set(KRATOS_ALE_APPLICATION_EXTRA_LIBS ${KRATOS_ALE_APPLICATION_EXTRA_LIBS} ${MPI_LIBRARIES}) # also is this needed? I seem to be the only one doing it...

    # How abt this line, do be need to link explicitly?
	set(KRATOS_ALE_APPLICATION_EXTRA_LIBS ${KRATOS_ALE_APPLICATION_EXTRA_LIBS} ${TRILINOS_LIBRARIES} ${MPI_LIBRARIES}) # Jordi why does it work without linking?
endif(${MPI_NEEDED} MATCHES ON )


###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosALEApplication SHARED ${KRATOS_ALE_APPLICATION_SOURCES})
target_link_libraries(KratosALEApplication KratosCore ${KRATOS_ALE_APPLICATION_EXTRA_LIBS})
set_target_properties(KratosALEApplication PROPERTIES PREFIX "")
install(TARGETS KratosALEApplication DESTINATION libs )

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosALEApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosALEApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add Cotire
if(USE_COTIRE MATCHES ON)
  cotire(KratosALEApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ALEApplication.py" DESTINATION KratosMultiphysics )
