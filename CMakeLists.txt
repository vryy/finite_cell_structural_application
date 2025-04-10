message("**** configuring KratosFiniteCellStructuralApplication ****")

######### add compilation units here
set( KRATOS_FINITE_CELL_STRUCTURAL_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/extrapolated_kinematic_linear.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/extrapolated_constant_stress_kinematic_linear.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_force_with_function.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_force_with_function.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/immersed_point_force.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ghost_penalty_displacement_gradient_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ghost_penalty_stress_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/finite_cell_structural_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_conditions_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/finite_cell_structural_python_application.cpp
)
#######################################################################

## define library Kratos which defines the basic python interface
add_library(KratosFiniteCellStructuralApplication SHARED ${KRATOS_FINITE_CELL_STRUCTURAL_APPLICATION_SOURCES})
if(MSVC)
    target_compile_options(KratosFiniteCellStructuralApplication PRIVATE $<$<CONFIG:Debug>:/bigobj>)
endif()
target_link_libraries(KratosFiniteCellStructuralApplication PUBLIC KratosCore)
target_link_libraries(KratosFiniteCellStructuralApplication PUBLIC KratosBRepApplication)
target_link_libraries(KratosFiniteCellStructuralApplication PUBLIC KratosStructuralApplication)
target_link_libraries(KratosFiniteCellStructuralApplication PUBLIC KratosFiniteCellApplication)
target_compile_definitions(KratosFiniteCellStructuralApplication PRIVATE KRATOS_CORE=IMPORT)
target_compile_definitions(KratosFiniteCellStructuralApplication PRIVATE FINITE_CELL_STRUCTURAL_APPLICATION=EXPORT)
target_compile_definitions(KratosFiniteCellStructuralApplication PRIVATE FINITE_CELL_APPLICATION=IMPORT)
target_compile_definitions(KratosFiniteCellStructuralApplication PRIVATE STRUCTURAL_APPLICATION=IMPORT)
target_compile_definitions(KratosFiniteCellStructuralApplication PRIVATE BREP_APPLICATION=IMPORT)
set_target_properties(KratosFiniteCellStructuralApplication PROPERTIES PREFIX "")
install(TARGETS KratosFiniteCellStructuralApplication DESTINATION libs )
#######################################################################

## changing the .dll suffix to .pyd for Windows compatibility
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosFiniteCellStructuralApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#######################################################################

## install application's python script if required
if(${INSTALL_PYTHON_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
#######################################################################

## Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FiniteCellStructuralApplication.py" DESTINATION KratosMultiphysics )
#######################################################################
