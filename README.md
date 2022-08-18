# WFES2

Wright-Fisher Exact Solver 2

Wright-Fisher Exact Solver (`WFES`) implements a variety of exact calculations with the Wright-Fisher model. Unlike other approaches, `WFES` does not use simulations or strong simplifying assumptions. `WFES` benefits from high-performance linear algebra techniques, making it possible to compute exact quantities for biologically realistic population sizes. The following document details the usage of the `WFES` code.

# Installation

*Note:* We are in the process of developing a version of WFES2 that fully takes advantage of Apple Silicon. In the meantime, the make system has been updated to compile on Apple Silicon using Rosetta2. The binary is automatically generated for Mac OS 10.4 x86_64

Please see [INSTALL.md](INSTALL.md) for installation instructions.

# Usage

Please consult the [manual](https://github.com/dekoning-lab/wfes2/blob/master/doc/manual.pdf) for information on usage and models.

Also note that for clarity, the allele frequency spectrum calculators have been renamed. `wfafle` is now known as `wfafs_deterministic` and `wfas` has been renamed `wfafs_stochastic`. The manual needs to be updated to reflect this change.

