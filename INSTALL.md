# WFES2

# Installation

Currently, `wfes2` is supported on `linux` and `macos` systems.

## Dependencies

### Compilation

You will need `cmake` in order to compile the code.

#### Setting the compiler

On some systems, it might be necessary to specify the compiler before running `cmake`. The preferred way to do this is with `CC` and `CXX` variables:

```bash
export CC=$(which gcc)
export CXX=$(which g++)
cmake ...
```

### Runtime

`WFES` depends on the INTEL `MKL` and `omp` libraries. They can be isntalled with the free installer from INTEL, or through the `conda` package manager.

#### `INTEL` installation

The required libraries can be found on the [intel website](https://software.intel.com/en-us/mkl). After download and installation, build the executables:

```bash
mkdir build && cd build
cmake -DINTEL=ON ..
make
```

By default, the `INTEL` installed puts the libraries and headers into `/opt/intel`. If the installation is located elsewhere, point to it with:

```bash
cmake -DINTEL=ON -DINTEL_ROOT=<INTEL PATH> ..
```

#### `conda` isntallation

The `conda` package manager provides the required packages. Note that `miniconda` will suffice for the installation. First, use `conda` to isntall the required packages:

```bash
conda install mkl mkl-include
```

Then, compile the executables with:

```bash
mkdir build && cd build
cmake -DCONDA=ON ..
make
```

By default, the build scripts will look for the libraries and headers in `$HOME/miniconda3/`. This setting can be changed:

```bash
cmake -DCONDA=ON -DCONDA_ROOT=<CONDA PATH> ..
```

#### `OpenMP`

Note that `OpenMP` will always be used for the `MKL` routines, which take care of the bulk of the computation. To use `OpenMP` for the matrix construction routines, it needs to be enabled before compilation:

``` bash
cmake <...> -DOMP=ON ..
```

This option is disabled by default.
