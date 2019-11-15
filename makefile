PLAT:=$(shell uname)

LIB:=-L${CONDA_PREFIX}/lib
LIB+=-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lmkl_avx -lmkl_vml_avx
LIB+=-liomp5 -lpthread
LIB+=-lm -ldl

INC:=-Iinclude
INC+=-I${CONDA_PREFIX}/include
LIB+=-I${CONDA_PREFIX}/include/eigen3

CXXFLAGS:=-std=c++11 -Wall -Wformat -Wno-deprecated-declarations
CXXFLAGS+=-DMKL_ILP64 -m64

# Platform-specific tweaks
ifeq (${PLAT},Linux)
	CXXFLAGS+=-DOMP -fopenmp
else ifeq (${PLAT},Darwin)
	LIB+=-rpath ${CONDA_PREFIX}/lib
endif

LIBFLAGS:=-fPIC -shared

TARGETS:=wfes_single wfes_switching wfafle wfes_sweep wfes_sequential wfas phase_type_dist phase_type_moments wfes_test

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
