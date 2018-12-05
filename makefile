LIB_DIR:=lib
BIN_DIR:=bin

INCD:=include
INTEL_ROOTD:=/opt/intel
MKL_ROOTD:=${INTEL_ROOTD}/mkl
MKL_LIBD:=${MKL_ROOTD}/lib
MKL_INCD:=${MKL_ROOTD}/include
INTEL_LIBD:=${INTEL_ROOTD}/lib
ALL_INCD:=$(addprefix -I, ${INCD} ${MKL_INCD})
ALL_LIBD:=$(addprefix -L, ${LIB_DIR} ${MKL_LIBD} ${INTEL_LIBD})

RPATH:=-Wl,-rpath,${MKL_LIBD},-rpath,${INTEL_LIBD}

MKL_LIBS:=mkl_intel_ilp64 mkl_intel_thread mkl_core
INTEL_LIBS:=iomp5 pthread
OTHER_LIBS:=m dl
ALL_LIBS:=$(addprefix -l, ${MKL_LIBS} ${INTEL_LIBS} ${OTHER_LIBS})

MKL_FLAGS:=-DMKL_ILP64 -m64
OPT_FLAGS:=-std=c++11
WARN_FLAGS:=-Wall -Wformat -Wno-deprecated-declarations
ALL_FLAGS:=${MKL_FLAGS} ${OPT_FLAGS} ${WARN_FLAGS}

EXE:=$(addprefix ${BIN_DIR}/, wfes_single wfes_switching wfes_sweep)

.PHONY: all prep clean
all: prep ${EXE}
prep:
	mkdir -p ${LIB_DIR} ${BIN_DIR}
clean:
	rm -rf ${LIB_DIR} ${BIN_DIR}

${LIB_DIR}/libwfes.so: src/lib/*.cpp
	${CXX} -shared $^ ${ALL_INCD} ${ALL_LIBD} ${ALL_LIBS} ${ALL_FLAGS} -o $@

${BIN_DIR}/%: src/%_exe.cpp ${LIB_DIR}/libwfes.so
	${CXX} $< ${ALL_INCD} ${ALL_LIBD} -lwfes ${RPATH} ${ALL_LIBS} ${ALL_FLAGS} -o $@
