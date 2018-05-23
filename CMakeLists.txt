cmake_minimum_required(VERSION 2.6)

project(WFAFLE)

set(CMAKE_CXX_STANDARD 11)

set(INTEL_ROOT "$ENV{HOME}/intel/")

if(UNIX AND NOT APPLE) #linux
    set(LIB_SUFFIX "/intel64")
endif()
if(APPLE)
    cmake_policy(SET CMP0042 NEW)
endif()

set(MKL_ROOT_DIR "${INTEL_ROOT}/mkl")
set(MKL_LIB_DIR "${MKL_ROOT_DIR}/lib${LIB_SUFFIX}")
set(MKL_INCLUDE_DIR "${MKL_ROOT_DIR}/include")
set(MKL_LIBS "-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core")

set(INTEL_OMP_LIB_DIRS "${INTEL_ROOT}/lib${LIB_SUFFIX}")
set(INTEL_OMP_LIBS "-liomp5 -lpthread")

set(OTHER_LIBS "-lm -ldl")

# add -fopenmp
set(CMAKE_CXX_FLAGS "-DMKL_ILP64 -m64 -Wall -Wformat -Wno-deprecated-declarations")

# find_package(Eigen3 3.3 REQUIRED)

file(GLOB LIB_SOURCES "src/*.cpp")
list(FILTER LIB_SOURCES EXCLUDE REGEX "src/.*_exe.cpp")

link_directories(${INTEL_OMP_LIB_DIRS} ${MKL_LIB_DIR})
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${MKL_INCLUDE_DIR})
include_directories("include")

add_library(wfes SHARED ${LIB_SOURCES})
target_link_libraries(wfes ${MKL_LIBS} ${INTEL_OMP_LIBS} ${OTHER_LIBS})

add_executable(wfafle "src/wfafle_exe.cpp")
target_link_libraries(wfafle wfes)

add_executable(wfes_soft "src/wfes_soft_exe.cpp")
target_link_libraries(wfes_soft wfes)

add_executable(wfes_single "src/wfes_single_exe.cpp")
target_link_libraries(wfes_single wfes)

add_executable(test_wfes "src/test_exe.cpp")
target_link_libraries(test_wfes wfes)

# install(TARGETS wfafle RUNTIME DESTINATION bin)
