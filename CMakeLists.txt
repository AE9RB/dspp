cmake_minimum_required (VERSION 2.8)

project (dspp)
enable_testing()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
    if (CMAKE_BUILD_TYPE MATCHES "Release")
        set(CMAKE_C_FLAGS "-flto ${CMAKE_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "-flto ${CMAKE_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "-flto ${CMAKE_EXE_LINKER_FLAGS}")
        set(CMAKE_SHARED_LINKER_FLAGS "-flto ${CMAKE_SHARED_LINKER_FLAGS}")
        set(CMAKE_STATIC_LINKER_FLAGS "-flto ${CMAKE_STATIC_LINKER_FLAGS}")
        set(CMAKE_MODULE_LINKER_FLAGS "-flto ${CMAKE_MODULE_LINKER_FLAGS}")
    endif()
endif()

### set up header groups for xcode and visual studio

file(GLOB_RECURSE DSPP_HEADERS "include/*.hpp")
file(GLOB_RECURSE BENCHTEST_HEADERS "benchtest/*.hpp")
source_group("dspp" FILES ${DSPP_HEADERS})
source_group("benchtest" FILES ${BENCHTEST_HEADERS})

### selftest benchtest

add_executable(test_benchtest "benchtest/selftest.cpp" ${BENCHTEST_HEADERS})
add_test(NAME run_test_benchtest COMMAND test_benchtest)

### selftest dspp

file(GLOB_RECURSE TEST_DSPP_SRC "test/*.cpp")
add_executable(test_dspp ${TEST_DSPP_SRC} ${DSPP_HEADERS})
target_include_directories(test_dspp PRIVATE include benchtest)
add_test(NAME run_test_dspp COMMAND test_dspp)
