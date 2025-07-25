# Data-driven optimal approximation on Hardy spaces in simply connected domains
This repository contains code for the manuscript 
> A. Borghi, T. Breiten, *Data-driven optimal approximation on Hardy spaces in simply connected domains*, [ArXiv](https://arxiv.org/abs/2507.15837) (2025)

For the experiments to work [chebfun](https://www.chebfun.org) needs to be installed.

The algorithm described in the manuscript can be found in `code/functions/algorithm1.m`. The functions `fdm_2d_matrix.m` and `fdm_2d_vector.m` can be found in the [LYAPACK toolbox](https://www.netlib.org/lyapack/) (see also [here](https://morwiki.mpi-magdeburg.mpg.de/morwiki/Convection-Diffusion#cite_note-lyapack-1)).

To get the results showed in the paper, the following files can be executed in the `code` folder:

### (Plots of the conformal mappings and reflections)
- `mappings.m`  

### (Heat equation example)
- `mobius_E2D_heat.m`  

### (Backward differentiation methods examples)
- `BDF2_test_disc.m` 
- `BDF4_test_disc.m` 

