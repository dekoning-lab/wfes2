# WFES2

# Installation

Currently, `wfes2` is supported on `linux` and `macos` systems. The
installation is done trhough `conda`.

# Conda setup

If you do no have `conda` installed, go
[here](https://docs.conda.io/en/latest/miniconda.html) and follow the
installation instructions for your platform. Otherwise, proceed to the next
step.

# Conda environment setup

Create a new isolated `conda` environment by using platform-appropriate file:

```
conda env create -f conda-env-linux.yaml # linux

conda env create -f conda-env-macos.yaml # macos
```

This should install the necessary libraries. Now, activate the environment with:

```
conda activate wfes
```

You would have to repeat this step for every new shell session you open.

## Install all the packages directly into `base`

Alternatively, if you do not want to use a separate conda environment, you can
install the dependencies into the `base` env:

```
#linux
conda install mkl mkl-include eigen gxx_linux-64
#macos
conda install mkl mkl-include eigen clangxx_osx-64
```

## MacOS caveats

Since MacOS thinks different, you woulb need to install the compiler yourself. In general, this should be sufficient:

```
xcode-select --install
```

# Compiling

With the dependencies in place, simply:

```
make all
```

This creates a `bin` directory, with multiple executables ready to use.

# Run the tests

The test target will be built by `make all`:

```
bin/wfes_test
```

If all is well, you should see this:

```
===============================================================================
All tests passed (XX assertions in xx test cases)
```

