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

set(CMAKE_CXX_FLAGS "-DMKL_ILP64 -m64 -Wall -Wformat -Wno-deprecated-declarations")

if(OMP)
    message(STATUS "Using OpenMP")
    string(APPEND CMAKE_CXX_FLAGS " -fopenmp -DOMP")
endif()

file(GLOB LIB_SOURCES "src/*.cpp")
list(FILTER LIB_SOURCES EXCLUDE REGEX "src/.*_exe.cpp")

link_directories(${INTEL_OMP_LIB_DIRS} ${MKL_LIB_DIR})
include_directories(${MKL_INCLUDE_DIR})
include_directories("include")

add_library(wfes SHARED ${LIB_SOURCES})
target_link_libraries(wfes ${MKL_LIBS} ${INTEL_OMP_LIBS} ${OTHER_LIBS})

add_executable(wfafle "src/wfafle_exe.cpp")
target_link_libraries(wfafle wfes)

add_executable(wfes_sweep "src/wfes_sweep_exe.cpp")
target_link_libraries(wfes_sweep wfes)

add_executable(wfes_single "src/wfes_single_exe.cpp")
target_link_libraries(wfes_single wfes)

add_executable(wfes_switching "src/wfes_switching_exe.cpp")
target_link_libraries(wfes_switching wfes)

add_executable(test_wfes "src/test_exe.cpp")
target_link_libraries(test_wfes wfes)

#add_executable(dummy "src/dummy_exe.cpp")
#target_link_libraries(dummy wfes)

# install(TARGETS wfafle RUNTIME DESTINATION bin)
