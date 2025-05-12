cmake_minimum_required(VERSION 3.24)
include(CheckCXXCompilerFlag)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_program(CCACHE_PROGRAM ccache)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

if(CCACHE_PROGRAM)
  message(STATUS "Using compiler cache")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif()
project(KaHIP C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) 

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# if no build mode is specified build in release mode
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# tweak compiler flags
CHECK_CXX_COMPILER_FLAG(-funroll-loops COMPILER_SUPPORTS_FUNROLL_LOOPS)
if(COMPILER_SUPPORTS_FUNROLL_LOOPS)
  add_definitions(-funroll-loops)
endif()
CHECK_CXX_COMPILER_FLAG(-fno-stack-limit COMPILER_SUPPORTS_FNOSTACKLIMITS)
if(COMPILER_SUPPORTS_FNOSTACKLIMITS)
  add_definitions(-fno-stack-limit)
endif()
CHECK_CXX_COMPILER_FLAG(-Wall COMPILER_SUPPORTS_WALL)
if(COMPILER_SUPPORTS_WALL)
  add_definitions(-Wall)
endif()
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
  add_definitions(-march=native)
endif()
CHECK_CXX_COMPILER_FLAG(-fpermissive COMPILER_SUPPORTS_FPERMISSIVE)
if(COMPILER_SUPPORTS_FPERMISSIVE)
  add_definitions(-fpermissive)
endif()

# check dependencies
find_package(MPI REQUIRED)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP support detected")
  add_definitions(${OpenMP_CXX_FLAGS})
else()
  message(WARNING "OpenMP not available, activating workaround")
  add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
  set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS "")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc)
endif()

find_package(MPI REQUIRED COMPONENTS CXX)

# 64 Bit option
option(64BITMODE "64 bit mode" ON)
if(64BITMODE)
  add_definitions("-DMODE64BITEDGES")
  add_definitions("-DPOINTER64=1")
endif()

option(64BITVERTEXMODE "64 bit mode" ON)
if(64BITVERTEXMODE)
  add_definitions("-DMODE64BITVERTICES")
endif()

# optimized output
option(OPTIMIZED_OUTPUT "optimized output" OFF)
if(OPTIMIZED_OUTPUT)
  add_definitions("-DKAFFPAOUTPUT")
endif()

include(FetchContent)

FetchContent_Declare(
        cpi
        GIT_REPOSITORY https://github.com/kurpicz/cpi.git
        GIT_TAG "origin/main"
)

FetchContent_MakeAvailable(cpi)

include_directories(${cpi_SOURCE_DIR})
include_directories(${MPI_CXX_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/app)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.0.3)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/io)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/partition)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/tools)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib/data_structure/compression_vectors)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/VieClus)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/VieClus/extern/KaHIP/lib)
include_directories(${PROJECT_SOURCE_DIR}/extern/robin-hood-hashing/src/include)

if (NOT TARGET FlatBuffers)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/flatbuffers
        ${CMAKE_CURRENT_BINARY_DIR}/flatbuffers-build
        EXCLUDE_FROM_ALL)
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/stxxl)

set(KAGEN_BUILD_APPS OFF)
set(KAGEN_BUILD_EXAMPLES OFF)
set(KAGEN_BUILD_TESTS OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/KaGen)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/VieClus)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/abseil)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc")

# add STXXL includes path
include_directories(${STXXL_INCLUDE_DIRS})

add_library(robin_hood INTERFACE)
target_include_directories(robin_hood INTERFACE
    ${PROJECT_SOURCE_DIR}/extern/robin-hood-hashing/src/include
)

set(LIBKAFFPA_SOURCE_FILES
  lib/tools/quality_metrics.cpp
  lib/tools/random_functions.cpp)
add_library(libkaffpa OBJECT ${LIBKAFFPA_SOURCE_FILES})
TARGET_LINK_LIBRARIES(libkaffpa PRIVATE cpi)
target_link_libraries(libkaffpa PRIVATE absl::flat_hash_map robin_hood)

set(ARGTABLE_SOURCE_FILES
extern/argtable3-3.0.3/argtable3.c)
add_library(argtab OBJECT ${ARGTABLE_SOURCE_FILES})

set(GRAPH_IO_SOURCE_FILES
        lib/io/graph_io_stream.h
        lib/io/graph_io_stream.cpp)
add_library(graphio OBJECT ${GRAPH_IO_SOURCE_FILES})
target_link_libraries(graphio PRIVATE cpi)
target_link_libraries(graphio PRIVATE absl::flat_hash_map robin_hood)

set(LIBSTREAMPARTITIONING_SOURCE_FILES
  lib/partition/onepass_partitioning/vertex_partitioning.cpp
  lib/partition/onepass_partitioning/floating_block.cpp
  lib/partition/onepass_partitioning/leiden.cpp
  lib/partition/onepass_partitioning/modularity.cpp
  lib/extclustering/extclustering_vieclus.cpp
  lib/extclustering/extclustering.cpp )
add_library(libstreampartitioning OBJECT ${LIBSTREAMPARTITIONING_SOURCE_FILES})
target_link_libraries(libstreampartitioning PRIVATE cpi )
target_link_libraries(libstreampartitioning PRIVATE absl::flat_hash_map robin_hood)

# generate targets for each binary

add_library(clustrelib STATIC
    app/clustering.cpp
    $<TARGET_OBJECTS:graphio>
    $<TARGET_OBJECTS:libstreampartitioning>
    $<TARGET_OBJECTS:libkaffpa>
    $<TARGET_OBJECTS:argtab>
)

target_include_directories(clustrelib PRIVATE extern/VieClus ${MPI_CXX_INCLUDE_DIRS})
target_compile_definitions(clustrelib PRIVATE "-DMODE_FREIGHT_GRAPHS")
target_link_libraries(clustrelib
        PRIVATE cpi
        graphio
        libkaffpa
        argtab
        VieClus
        absl::flat_hash_map
        robin_hood
        libstreampartitioning
        flatbuffers
        ${STXXL_LIBRARIES}
        ${OpenMP_CXX_LIBRARIES}
        ${MPI_CXX_LIBRARIES})
install(TARGETS clustrelib DESTINATION bin)

add_executable(clustre app/clustering.cpp)
target_include_directories(clustre PRIVATE extern/VieClus ${MPI_CXX_INCLUDE_DIRS})
target_compile_definitions(clustre PRIVATE "-DMODE_FREIGHT_GRAPHS")
target_link_libraries(clustre
        PRIVATE cpi
        graphio
        libkaffpa
        argtab
        VieClus
        absl::flat_hash_map
        robin_hood
        libstreampartitioning
        flatbuffers
        ${STXXL_LIBRARIES}
        ${OpenMP_CXX_LIBRARIES}
        ${MPI_CXX_LIBRARIES})
install(TARGETS clustre DESTINATION bin)