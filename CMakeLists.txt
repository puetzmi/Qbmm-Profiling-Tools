#############################################################################
# BASIC PROJECT INFORMATION
#############################################################################
cmake_minimum_required(VERSION 3.21.0)
project(QbmmProfilingTools)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
#############################################################################
#############################################################################


#############################################################################
# BUILD TYPE CONFIGURATION
#############################################################################
set(DEFAULT_BUILD_TYPE "Release")
set(VALID_BUILD_TYPES "Release" "Debug" "RelWithDebInfo" "MinSizeRelease")

# Check if build type is empty
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    # Set to default
    message(STATUS "Setting build type to default '${DEFAULT_BUILD_TYPE}'")
    set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE
      STRING "Choose the type of build." FORCE)
# Make sure build type is valid
elseif (NOT CMAKE_BUILD_TYPE IN_LIST VALID_BUILD_TYPES)
    message(FATAL_ERROR "Invalid build type: CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
endif()
# Set options for CMake Curses
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${VALID_BUILD_TYPES})
#############################################################################
#############################################################################


#############################################################################
# COMPILER CONFIGURATION
#############################################################################
find_package(OpenMP REQUIRED)
find_package(MPI)

add_compile_options(${OpenMP_CXX_FLAGS} "-Wall" "-Wextra")
add_link_options(${OpenMP_CXX_FLAGS} "-Wall" "-Wextra")

# AVX2 support
option( USE_AVX2 "Use Advanced Vector Extensions (AVX2)" ON )
if (USE_AVX2)
    add_compile_options("-mavx2")
endif()

if (NOT GDB_DEBUG_LEVEL OR NOT CMAKE_BUILD_TYPE MATCHES "Debug")
    set(GDB_DEBUG_LEVEL "0" CACHE STRING
        "Level of debugging information for GDB (only for build type 'Debug')" FORCE)
endif()
if (CMAKE_BUILD_TYPE MATCHES "Debug")
    set_property(CACHE GDB_DEBUG_LEVEL PROPERTY STRINGS "0" "1" "2" "3")
else()
    set_property(CACHE GDB_DEBUG_LEVEL PROPERTY STRINGS "0")
endif()
if (NOT GDB_DEBUG_LEVEL MATCHES "0")
    #set(CMAKE_CXX_FLAGS_DEBUG "-ggdb${GDB_DEBUG_LEVEL}")
    add_compile_options("-ggdb${GDB_DEBUG_LEVEL}")
    #set(CMAKE_C_FLAGS_DEBUG "-ggdb${GDB_DEBUG_LEVEL}")
endif()

#message(STATUS "${CMAKE_CXX_COMPILER}-${CMAKE_CXX_COMPILER_VERSION}")
# Check computer architecture and add definition accordingly
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    add_compile_definitions(ARCH_X86_64)
else()
    add_compile_definitions(ARCH_X86)
endif()
#############################################################################
#############################################################################


#############################################################################
# ADDITIONAL DEFINITIONS
#############################################################################
if(USE_AVX2)
    set(MALLOC_ALIGN_SIZE 32 CACHE STRING
        "Alignment of dynamically allocated memory in bytes" FORCE)
else()
    set(MALLOC_ALIGN_SIZE 16 CACHE STRING
        "Alignment of dynamically allocated memory in bytes" FORCE)
endif()
set_property(CACHE MALLOC_ALIGN_SIZE PROPERTY STRINGS "16" "32" "64" "128")
add_compile_definitions(MALLOC_ALIGN_SIZE=${MALLOC_ALIGN_SIZE})
#############################################################################
#############################################################################

#############################################################################
# BUILD ID ASSIGNMENT USING GIT TAG/COMMIT ##################################
#############################################################################
find_package(Git REQUIRED)
add_custom_target(version ALL
  ${CMAKE_COMMAND} -D SRC=${CMAKE_CURRENT_SOURCE_DIR}/include/build_id.hpp.in
                   -D DST=${CMAKE_CURRENT_BINARY_DIR}/include/build_id.hpp
                   -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
                   -D COMPILER_INFO=${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}-${CMAKE_BUILD_TYPE}
                   -P ${CMAKE_CURRENT_SOURCE_DIR}/GenerateBuildIdHeader.cmake
)



#############################################################################
# BUILD
#############################################################################
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel" OR CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
    set(MKL_MPI intelmpi)
    set(MKL_THREADING intel_thread)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(MKL_MPI openmpi)
    set(MKL_THREADING gnu_thread)
endif()
set(MKL_INTERFACE lp64)
find_package(MKL CONFIG REQUIRED)

# Set input/output directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src
    $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/include
)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Build shared libraries
add_library(LinAlg SHARED src/cholesky.cpp src/linsolve.cpp src/eigensolve.cpp)
add_library(Qbmm SHARED src/qmom.cpp src/core_inversion.cpp)
add_library(QbmmProfiling SHARED src/qmom_profiler.cpp)
add_library(PhysicalModels SHARED src/physical_models.cpp)

# Build and link executables
add_executable(core_inversion_benchmark applications/core_inversion_benchmark.cpp)
target_link_libraries(core_inversion_benchmark Qbmm LinAlg $<LINK_ONLY:MKL::MKL>)

add_executable(eigensolver_benchmark_general applications/eigensolver_benchmark_general.cpp)
target_link_libraries(eigensolver_benchmark_general LinAlg $<LINK_ONLY:MKL::MKL>)

add_executable(linsolver_benchmark applications/linsolver_benchmark.cpp)
target_link_libraries(linsolver_benchmark LinAlg $<LINK_ONLY:MKL::MKL>)

add_executable(eigensolver_benchmark_moments applications/eigensolver_benchmark_moments.cpp)
target_link_libraries(eigensolver_benchmark_moments Qbmm LinAlg $<LINK_ONLY:MKL::MKL>)

#############################################################################
#############################################################################


#############################################################################
# TESTING
#############################################################################
if (NOT DEFINED WITH_TESTING)
    #set(WITH_TESTING ON CACHE
    #    STRING "Select whether or not test application and targets shall be built" FORCE)
    option( WITH_TESTING "Select whether or not test application and targets shall be built" ON)
endif()
if (WITH_TESTING)
    enable_testing()

    add_executable(test_build_id tests/test_build_id.cpp)
    add_test(NAME test_build_id COMMAND test_build_id)

    add_executable(test_cholesky tests/test_cholesky.cpp)
    target_link_libraries(test_cholesky LinAlg $<LINK_ONLY:MKL::MKL>)
    add_test(NAME test_cholesky COMMAND test_cholesky)

    add_executable(test_eigensolve tests/test_eigensolve.cpp)
    target_link_libraries(test_eigensolve LinAlg $<LINK_ONLY:MKL::MKL>)
    add_test(NAME test_eigensolve COMMAND test_eigensolve)

    add_executable(test_linsolve tests/test_linsolve.cpp)
    target_link_libraries(test_linsolve LinAlg $<LINK_ONLY:MKL::MKL>)
    add_test(NAME test_linsolve COMMAND test_linsolve)

    add_executable(test_core_inversion tests/test_core_inversion.cpp)
    target_link_libraries(test_core_inversion Qbmm LinAlg $<LINK_ONLY:MKL::MKL>)
    add_test(NAME test_core_inversion COMMAND test_core_inversion)

    add_executable(test_qmom tests/test_qmom.cpp)
    target_link_libraries(test_qmom Qbmm LinAlg $<LINK_ONLY:MKL::MKL>)
    add_test(NAME test_qmom COMMAND test_qmom)

    # Enable test for memory leaks
    include(CTest)
    add_custom_target(test_memcheck
        COMMAND ${CMAKE_CTEST_COMMAND}
            --force-new-ctest-process --test-action memcheck
        COMMAND cat "${CMAKE_BINARY_DIR}/Testing/Temporary/MemoryChecker.*.log"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
endif()
#############################################################################
#############################################################################


#############################################################################
# MAKE DOXYGEN DOCUMENTATION IF SELECTED
#############################################################################
option(BUILD_DOC "Build documentation" OFF)
if (BUILD_DOC)
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
        set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
        configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

        # note the option ALL which allows to build the docs together with the application
        add_custom_target( doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
        message("Doxygen not found.")
    endif (DOXYGEN_FOUND)
endif (BUILD_DOC)
#############################################################################
#############################################################################

# This is useful to print all variable when needed
#if (TRUE)
if (FALSE)
get_cmake_property(_variableNames VARIABLES)
list (SORT _variableNames)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
endif()
