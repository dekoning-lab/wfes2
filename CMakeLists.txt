cmake_minimum_required(VERSION 2.6)

project(WFES)

set(CMAKE_CXX_STANDARD 11)

if(APPLE)
    cmake_policy(SET CMP0042 NEW)
endif()

# CONDA installation
if(NOT DEFINED CONDA_ROOT)
    set(CONDA_ROOT "$ENV{HOME}/miniconda3")
endif()

if(CONDA)
    message(STATUS "Using ${CONDA_ROOT}")

    set(MKL_LIB_DIR "${CONDA_ROOT}/lib")
    set(MKL_INCLUDE_DIR "${CONDA_ROOT}/include")
    set(INTEL_OMP_LIB_DIRS "${CONDA_ROOT}/lib")
endif()

# INTEL installation
if(NOT DEFINED INTEL_ROOT)
    set(INTEL_ROOT "/opt/intel")
endif()

if(INTEL)
    message(STATUS "Using ${INTEL_ROOT}")
    if(UNIX AND NOT APPLE) #linux
        set(LIB_SUFFIX "/intel64")
    endif()

    set(MKL_ROOT_DIR "${INTEL_ROOT}/mkl")
    set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib${LIB_SUFFIX}")
    set(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include")

    set(INTEL_OMP_LIB_DIRS "${INTEL_ROOT}/lib${LIB_SUFFIX}")
endif()

set(MKL_LIBS "-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core")
set(INTEL_OMP_LIBS "-liomp5 -lpthread")
set(OTHER_LIBS "-lm -ldl")

set(CMAKE_CXX_FLAGS "-std=c++11 -DMKL_ILP64 -m64 -Wall -Wformat -Wno-deprecated-declarations")

if(OMP)
    message(STATUS "Using OpenMP")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DOMP")
endif()

file(GLOB LIB_SOURCES "src/lib/*.cpp")

link_directories(${INTEL_OMP_LIB_DIRS} ${MKL_LIB_DIR})
include_directories(${MKL_INCLUDE_DIR})
include_directories("include")

add_library(wfes SHARED ${LIB_SOURCES})
target_link_libraries(wfes ${MKL_LIBS} ${INTEL_OMP_LIBS} ${OTHER_LIBS})

add_executable(wfafs_deterministic "src/wfafle.cpp")
target_link_libraries(wfafs_deterministic wfes)

add_executable(wfes_sweep "src/wfes_sweep.cpp")
target_link_libraries(wfes_sweep wfes)

add_executable(wfes_single "src/wfes_single.cpp")
target_link_libraries(wfes_single wfes)

add_executable(wfes_switching "src/wfes_switching.cpp")
target_link_libraries(wfes_switching wfes)

add_executable(wfes_sequential "src/wfes_sequential.cpp")
target_link_libraries(wfes_sequential wfes)

add_executable(test_wfes "src/test.cpp")
target_link_libraries(test_wfes wfes)

add_executable(phase_type_dist "src/phase_type_dist.cpp")
target_link_libraries(phase_type_dist wfes)

add_executable(phase_type_moments "src/phase_type_moments.cpp")
target_link_libraries(phase_type_moments wfes)

add_executable(wfafs_stochastic "src/wfas.cpp")
target_link_libraries(wfafs_stochastic wfes)

# install(TARGETS wfafle RUNTIME DESTINATION bin)
