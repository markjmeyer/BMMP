# BMMP

Code for estimating Bayesian methods for multivariate matched proportions (MMP) as described by Meyer, Cheng, and Knutson (2022), currently under review, but available at https://arxiv.org/abs/2108.04096.

## Folder Details

The folder titled Code contains the file bmmp.R which has all of the relevant functions for estimating Bayesian MMPs as well as non-Bayesian methods from the existing literature. This file also contains the data used in the reanalysis described in Meyer, Cheng, and Knutson (2022). The folder titled Application contains the file soc_reanalysis.R which can be used to replicate the reanalysis. The files in the folder Simulation were used to for the Simulation Study. Files ending in _mvp can take several days to run.

### Brief Function and Object Description

All functions two matrices, X1 and X2, as arguments. X1 is an n x K matrix of the first set of binary multivariate outcomes. X2 is an n x K matrix of the second or paired set of binary multivariate outcomes. Each function has additional arguments for prior specification and output controls. Priors default to those described in the manuscript.


