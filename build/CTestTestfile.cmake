# CMake generated Testfile for 
# Source directory: /mnt/d/SHiP/ej200_v3
# Build directory: /mnt/d/SHiP/ej200_v3/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(smoke_test "/mnt/d/SHiP/ej200_v3/build/ej200_bar_sim" "-m" "macros/test.mac")
set_tests_properties(smoke_test PROPERTIES  PASS_REGULAR_EXPRESSION "=== EJ-200 Bar Run Summary ===" TIMEOUT "120" _BACKTRACE_TRIPLES "/mnt/d/SHiP/ej200_v3/CMakeLists.txt;20;add_test;/mnt/d/SHiP/ej200_v3/CMakeLists.txt;0;")
