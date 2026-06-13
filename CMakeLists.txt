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

find_package(Python3 REQUIRED COMPONENTS Interpreter)

option(USE_GDML "Build optional GDML geometry exporter" OFF)

if(USE_GDML AND Geant4_gdml_FOUND)
    add_executable(export_endonly_gdml tests/export_endonly_gdml.cc ${sources} ${headers})
    configure_ej200_target(export_endonly_gdml)
    target_compile_definitions(export_endonly_gdml PRIVATE USE_GDML)
    add_test(NAME export_endonly_gdml
             COMMAND export_endonly_gdml ${CMAKE_CURRENT_BINARY_DIR}/endonly_geometry.gdml)
elseif(USE_GDML AND NOT Geant4_gdml_FOUND)
    message(FATAL_ERROR "USE_GDML=ON but this Geant4 build has no GDML support.")
endif()
