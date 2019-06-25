LIB:=-L${CONDA_PREFIX}/lib
LIB+=-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -lmkl_avx -lmkl_vml_avx
LIB+=-liomp5 -lpthread
LIB+=-lm -ldl

INC:=-Iinclude
INC+=-I${CONDA_PREFIX}/include
LIB+=-I${CONDA_PREFIX}/include/eigen3

CXXFLAGS:=-std=c++11 -Wall -Wformat -Wno-deprecated-declarations
CXXFLAGS+=-DMKL_ILP64 -m64
CXXFLAGS+=-DOMP -fopenmp

LIBFLAGS:=-fPIC -shared

TARGETS:=wfes_single wfes_switching wfafle wfes_sweep wfes_sequential wfas phase_type_dest phase_type_moments wfes_test

.PHONY: all PRE clean

PRE:
	mkdir -p bin

all: $(addprefix bin/,${TARGETS})

bin/libwfes.so: src/lib/*
	${CXX} ${CXXFLAGS} ${INC} ${LIB} ${LIBFLAGS} $^ -o $@

bin/%: src/lib/* src/%.cpp
	@echo "# Building $(@F)"
	${CXX} ${CXXFLAGS} ${INC} ${LIB} $^ -o $@

clean:
	rm -r bin
