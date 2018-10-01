MKL_HOME := /cm/shared/apps/intel/compilers_and_libraries_2018.2.199/linux/mkl
WSMP_HOME := $(HOME)/wsmp/lib

LIBS := -lm -lpthread -lgfortran -ldl
MKL_LIBS := $(addprefix $(MKL_HOME)/lib/intel64/libmkl_, intel_ilp64.a sequential.a core.a)
MKL_INC := $(MKL_HOME)/include

 #blas95_lp64.a?

SRC_PATH := ../src
INC_PATH := ../include
DEP := PardisoSolver.cpp SparseMatrix.cpp WrightFisher.cpp parsing.cpp rdist.cpp util.cpp
SRC := $(addprefix $(SRC_PATH)/, $(DEP))

TARGETS := wfes_single wfes_switching wfes_sweep wfafle

all: $(TARGETS)

%: $(SRC) $(SRC_PATH)/%_exe.cpp
	g++ -o $@ $^ --std=c++11 -I$(INC_PATH) -I$(MKL_INC) -Wl,--start-group $(MKL_LIBS) $(WSMP_HOME)/libwsmp8_8.a -Wl,--end-group $(LIBS) -m64 -DMKL_ILP64

clean: 
	rm -rf wfes_single
