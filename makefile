LIBD:=lib
BIND:=bin
INCD:=include

uname:=$(shell uname)
ifeq (${uname},Linux)
INTEL_LIB_SUFFIX:=/intel64
endif

pwd:=$(shell pwd)

INTEL_ROOTD:=/opt/intel
MKL_ROOTD:=${INTEL_ROOTD}/mkl
MKL_LIBD:=${MKL_ROOTD}/lib${INTEL_LIB_SUFFIX}
MKL_INCD:=${MKL_ROOTD}/include
INTEL_LIBD:=${INTEL_ROOTD}/lib${INTEL_LIB_SUFFIX}
ALL_INCD:=$(addprefix -I, ${INCD} ${MKL_INCD})
ALL_LIBD:=$(addprefix -L, ${LIBD} ${MKL_LIBD} ${INTEL_LIBD})

RPATH:=-Wl,-rpath,${MKL_LIBD},-rpath,${INTEL_LIBD},-rpath,${pwd}/${LIBD}

MKL_LIBS:=mkl_intel_ilp64 mkl_intel_thread mkl_core
INTEL_LIBS:=iomp5 pthread
OTHER_LIBS:=m dl
ALL_LIBS:=$(addprefix -l, ${MKL_LIBS} ${INTEL_LIBS} ${OTHER_LIBS})

MKL_FLAGS:=-DMKL_ILP64 -m64
OPT_FLAGS:=-std=c++11
WARN_FLAGS:=-Wall -Wformat -Wno-deprecated-declarations
EXTRA_FLAGS:=

ifeq (${OMP},1)
EXTRA_FLAGS+=-fopenmp -DOMP
endif
ifeq (${DEBUG},1)
EXTRA_FLAGS+=-g
endif

ALL_FLAGS:=${MKL_FLAGS} ${OPT_FLAGS} ${WARN_FLAGS} ${EXTRA_FLAGS}

EXE:=$(addprefix ${BIND}/, wfes_single wfes_switching wfes_allele_freq wfes_sweep wfes_phase_type_dist wfes_test)

.PHONY: all prep clean install
all: prep ${EXE}
prep:
	mkdir -p ${LIBD} ${BIND}
clean:
	rm -rf ${LIBD} ${BIND}

${LIBD}/libwfes.so: src/lib/*.cpp
	${CXX} -shared -fPIC $^ ${ALL_INCD} ${ALL_LIBD} ${ALL_LIBS} ${ALL_FLAGS} -o $@

${BIND}/%: src/exe/%.cpp ${LIBD}/libwfes.so
	${CXX} $< ${ALL_INCD} ${ALL_LIBD} -lwfes ${RPATH} ${ALL_LIBS} ${ALL_FLAGS} -o $@

install: all
	cp bin/wfes_* /usr/local/bin
	cp lib/libwfes.so /usr/local/lib

uninstall:
	rm -f /usr/local/bin/wfes_*
	rm -f /usr/local/lib/libwfes.so
