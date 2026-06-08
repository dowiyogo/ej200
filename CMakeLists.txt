cmake_minimum_required(VERSION 3.16)
project(ej200_bar_sim VERSION 2.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(ej200_bar_sim main.cc ${sources} ${headers})
target_include_directories(ej200_bar_sim PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ej200_bar_sim ${Geant4_LIBRARIES})

file(COPY macros DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()
add_test(NAME smoke_test COMMAND ej200_bar_sim -m macros/test.mac)
set_tests_properties(smoke_test PROPERTIES
    PASS_REGULAR_EXPRESSION "=== EJ Scintillator Bar Run Summary ==="
    TIMEOUT 120)

add_test(NAME edge_scan_smoke COMMAND ej200_bar_sim -m macros/scan_edge.mac)
set_tests_properties(edge_scan_smoke PROPERTIES
    PASS_REGULAR_EXPRESSION "=== EJ Scintillator Bar Run Summary ==="
    TIMEOUT 1200)

add_executable(physics_baseline_check tests/physics_baseline_check.cc ${sources} ${headers})
target_include_directories(physics_baseline_check PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(physics_baseline_check ${Geant4_LIBRARIES})

add_test(NAME physics_baseline_check COMMAND physics_baseline_check)

add_executable(readout_config_check tests/readout_config_check.cc ${sources} ${headers})
target_include_directories(readout_config_check PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(readout_config_check ${Geant4_LIBRARIES})

add_test(NAME readout_config_check COMMAND readout_config_check)
