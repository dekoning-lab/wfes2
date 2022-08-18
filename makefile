PLAT:=$(shell uname)

LIB:=-L${CONDA_PREFIX}/lib
LIB+=-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core 
LIB+=-liomp5 -lpthread
LIB+=-lm -ldl

INC:=-Iinclude
INC+=-I${CONDA_PREFIX}/include
LIB+=-I${CONDA_PREFIX}/include/eigen3

CXXFLAGS:=-std=c++11 -Wall -Wformat -Wno-deprecated-declarations -fpermissive
CXXFLAGS+=-DMKL_ILP64 -m64

# Platform-specific tweaks
ifeq (${PLAT},Linux)
	CXXFLAGS+=-DOMP -fopenmp
else ifeq (${PLAT},Darwin)
	LIB+=-lmkl_avx -lmkl_vml_avx
	LIB+=-rpath ${CONDA_PREFIX}/lib
	CXXFLAGS+=-target x86_64-apple-macos10.14
endif

LIBFLAGS:=-fPIC -shared

TARGETS:=wfes_single wfes_switching wfes_sweep wfes_sequential wfafs_deterministic wfafs_stochastic phase_type_dist phase_type_moments time_dist time_dist_sgv wfes_test

.PHONY: all PRE clean

all: PRE $(addprefix bin/,${TARGETS})

PRE:
	mkdir -p bin


bin/libwfes.so: src/lib/*
	mkdir -p bin
	${CXX} ${CXXFLAGS} ${INC} ${LIB} ${LIBFLAGS} $^ -o $@ ${EXTRAFLAGS}

bin/%: src/lib/* src/%.cpp
	mkdir -p bin
	@echo "# Building $(@F)"
	${CXX} ${CXXFLAGS} ${INC} ${LIB} $^ -o $@ ${EXTRAFLAGS}

clean:
	rm -r bin
