cmake_minimum_required(VERSION 3.20)

project(spark-benchmark)

set(CMAKE_CXX_STANDARD 20)
set(FETCHCONTENT_QUIET OFF CACHE BOOL "" FORCE)

# Add CPM
set(CPM_DOWNLOAD_VERSION 0.40.2)
set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
if(NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
    message(STATUS "Downloading CPM.cmake...")
    file(DOWNLOAD https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake ${CPM_DOWNLOAD_LOCATION})
endif()
include(${CPM_DOWNLOAD_LOCATION})

# Add external libs
CPMAddPackage("gh:d99kris/rapidcsv@8.83")
CPMAddPackage("gh:p-ranav/argparse@3.1")
add_subdirectory(spark)

# Source files
add_executable(spark-benchmark
        src/parameters.cpp
        src/parameters.h
        src/reactions.cpp
)

target_link_libraries(spark-benchmark PUBLIC spark::spark rapidcsv argparse)
target_include_directories(spark-benchmark PRIVATE ${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/external)
