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

file(GLOB_RECURSE INCLUDE_DSPP "include/dspp/*.hpp")
source_group("Header Files\\dspp" FILES ${INCLUDE_DSPP})
LIST(APPEND INCLUDE_DSPP "include/dspp.hpp")

file(GLOB_RECURSE INCLUDE_BENCHTEST "include/benchtest/*.hpp")
source_group("Header Files\\benchtest" FILES ${INCLUDE_BENCHTEST})
LIST(APPEND INCLUDE_BENCHTEST "include/benchtest.hpp")

include_directories(include)

### selftest benchtest

file(GLOB_RECURSE SRC_TEST_BENCHTEST "src/test_benchtest/*.cpp")
add_executable(test_benchtest ${SRC_TEST_BENCHTEST} ${INCLUDE_BENCHTEST})
add_test(NAME run_test_benchtest COMMAND test_benchtest)

### selftest dspp

file(GLOB_RECURSE SRC_TEST_DSPP "src/test_dspp/*.cpp")
add_executable(test_dspp ${SRC_TEST_DSPP} ${INCLUDE_DSPP} ${INCLUDE_BENCHTEST})
add_test(NAME run_test_dspp COMMAND test_dspp)

### documentation

find_package(Doxygen)
if(DOXYGEN_FOUND)
add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/Doxyfile
    COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)
