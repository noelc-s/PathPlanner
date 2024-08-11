cmake_minimum_required(VERSION 3.10)

project(graph)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -pthread")

find_package(osqp)
if (NOT TARGET osqpstatic)  # This is what osqp-cpp looks for
    message(STATUS "osqpstatic target not found, checking osqp::osqpstatic...")
    if (NOT TARGET osqp::osqpstatic)
        message(STATUS "Could not find osqp::osqpstatic, osqp-cpp will try to fetch it...")
    else ()
        message(STATUS "Found osqp::osqpstatic. Using it in osqp-cpp.")
        add_library(osqpstatic ALIAS osqp::osqpstatic)
    endif ()
endif ()

include(FetchContent)

set(OSQP-CPP_BUILD_TESTS OFF CACHE BOOL "Build OSQP-CPP without tests")
set(ABSL_PROPAGATE_CXX_STD ON CACHE BOOL "Variable for OSQP-CPP ABSL dep")

FetchContent_Declare(
        osqp-cpp
        GIT_REPOSITORY https://github.com/google/osqp-cpp.git
        GIT_TAG 43433736334d6b515ea4b0247156fea9e56c0d3f
)

FetchContent_MakeAvailable(osqp-cpp)
# find_package(matplotlib_cpp HINTS /home/noelcs/repos/matplotlib-cpp_install REQUIRED)

add_executable(graph src/feasibility_test.cpp)
# target_include_directories(graph PUBLIC /home/noelcs/repos/matplotlib-cpp_install/include)
# target_include_directories(graph PUBLIC /usr/include/python3.8/)
target_link_libraries(graph PUBLIC osqp-cpp)