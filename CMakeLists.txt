cmake_minimum_required(VERSION 3.10)

project(path_planning)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -pthread -pg -fopenmp")

set(bezier_DIR /home/noelcs/repos/BezierTubes/cpp CACHE PATH "Path to BezierTubes dir")

### OSQP
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
###

### Bezier Tubes
add_library(bezier SHARED IMPORTED)
set_property(TARGET bezier PROPERTY IMPORTED_LOCATION ${bezier_DIR}/build/libBEZIER.so)
include_directories(bezier INTERFACE ${bezier_DIR}/inc)
####

### CUDA
find_package(Eigen3 REQUIRED)
MESSAGE( STATUS ">> --------------  USING CUDA --------------" )
# set( CUDA_TOOLKIT_ROOT_DIR "/opt/cuda" )
# set(CMAKE_C_COMPILER   /opt/cuda/bin/gcc)
# set(CMAKE_CXX_COMPILER /opt/cuda/bin/g++)

find_package(CUDA REQUIRED)
add_definitions(-DUSE_CUDA)


find_package(CUDA REQUIRED)
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_89")

set(CUDA_SOURCES
    src/kernel.cu
)

cuda_add_library(cuda_kernel STATIC ${CUDA_SOURCES})
target_include_directories(cuda_kernel PUBLIC inc/)
target_link_libraries(cuda_kernel Eigen3::Eigen)

# add_executable(cuda src/main.cpp ${SOURCE_FILES})
# cuda_add_executable(cuda src/main.cpp ${CUDA_SOURCES} src/obstacle.cpp)
# target_link_libraries(cuda ${CUDA_LIBRARIES})

# target_include_directories( cuda PUBLIC inc/)
# target_link_libraries(cuda Eigen3::Eigen)
###

add_executable(test src/main.cpp src/utils.cpp src/obstacle.cpp)
target_include_directories(test PUBLIC inc/)
target_link_libraries(test Eigen3::Eigen yaml-cpp)

find_package(yaml-cpp REQUIRED)

add_library(graph STATIC src/obstacle.cpp src/pathPlanner.cpp src/graph.cpp src/mpc.cpp src/utils.cpp)
target_link_libraries(graph PUBLIC osqp-cpp bezier yaml-cpp cuda_kernel)

add_executable(path_planning src/path_planning.cpp)
target_link_libraries(path_planning PUBLIC graph)

add_executable(integrator src/integrator_test.cpp)
target_link_libraries(integrator PUBLIC graph)

# For some reason cuda cannot compile with gcc > 8
# export CC=/usr/bin/gcc-8
# export CXX=/usr/bin/g++-8