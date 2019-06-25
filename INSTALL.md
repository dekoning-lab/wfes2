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

# Compiling

With the dependencies in place, simply:

```
make all
```

This creates a `bin` directory, with multiple executables ready to use.

# Run the test

The test target will be built by `make all`:

```
bin/wfes_test
```

If all is well, you should see this:

```
===============================================================================
All tests passed (XX assertions in xx test cases)
```

