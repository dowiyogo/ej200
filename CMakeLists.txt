cmake_minimum_required(VERSION 3.16)
project(ej200_bar_sim VERSION 2.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

add_subdirectory(src/external/OPSimTool)
add_subdirectory(src/external/SSLG4)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

function(configure_ej200_target target)
    target_include_directories(${target} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    target_compile_definitions(${target} PRIVATE EJ200_DATA_DIR="${PROJECT_SOURCE_DIR}/data")
    target_link_libraries(${target} ${Geant4_LIBRARIES} OPSim sslg4)
endfunction()

add_executable(ej200_bar_sim main.cc ${sources} ${headers})
configure_ej200_target(ej200_bar_sim)

file(COPY macros DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()
add_test(NAME smoke_test COMMAND ej200_bar_sim -m macros/test.mac)
set_tests_properties(smoke_test PROPERTIES
    PASS_REGULAR_EXPRESSION "=== EJ Scintillator Bar Run Summary ==="
    TIMEOUT 120)

add_test(NAME edge_scan_smoke COMMAND ej200_bar_sim -m macros/edge_scan_smoke.mac)
set_tests_properties(edge_scan_smoke PROPERTIES
    PASS_REGULAR_EXPRESSION "=== EJ Scintillator Bar Run Summary ==="
    TIMEOUT 1200)

add_executable(physics_baseline_check tests/physics_baseline_check.cc ${sources} ${headers})
configure_ej200_target(physics_baseline_check)

add_test(NAME sslg4_properties_check COMMAND physics_baseline_check)

add_executable(readout_config_check tests/readout_config_check.cc ${sources} ${headers})
configure_ej200_target(readout_config_check)

add_test(NAME readout_config_check COMMAND readout_config_check)

add_executable(export_endtop_gdml tests/export_endtop_gdml.cc ${sources} ${headers})
configure_ej200_target(export_endtop_gdml)

add_test(NAME export_endtop_gdml
         COMMAND export_endtop_gdml ${CMAKE_CURRENT_BINARY_DIR}/endtop_geometry.gdml)

find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME endtop_gdml_check
         COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/check_endtop_gdml.py
                 ${CMAKE_CURRENT_BINARY_DIR}/endtop_geometry.gdml)
set_tests_properties(endtop_gdml_check PROPERTIES DEPENDS export_endtop_gdml)

add_test(NAME endtop_balance_smoke
         COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/check_endtop_balance.py
                 ${CMAKE_CURRENT_BINARY_DIR}/ej200_bar_sim
                 ${CMAKE_CURRENT_BINARY_DIR}/macros/endtop_smoke_center.mac)
set_tests_properties(endtop_balance_smoke PROPERTIES TIMEOUT 300)
