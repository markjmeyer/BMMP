# BMMP

Code for estimating Bayesian methods for multivariate matched proportions (MMP) as described by Meyer, Cheng, and Knutson (2022), currently under review, but available at https://arxiv.org/abs/2108.04096.

## Folder Details

The folder titled Code contains the file bmmp.R which has all of the relevant functions for estimating Bayesian MMPs as well as non-Bayesian methods from the existing literature. This file also contains the data used in the reanalysis described in Meyer, Cheng, and Knutson (2022). The folder titled Application contains the file soc_reanalysis.R which can be used to replicate the reanalysis. The files in the folder Simulation were used to for the Simulation Study. Files ending in _mvp can take several days to run.

## Brief Function and Object Description

All functions two matrices, X1 and X2, as arguments. X1 is an n x K matrix of the first set of binary multivariate outcomes. X2 is an n x K matrix of the second or paired set of binary multivariate outcomes. Each function has additional arguments for prior specification and output controls. Priors default to those described in the manuscript.

### bpm: Bayesian Probit Model

> bpm(X1, X2, B = 10000, burnin = B, 
>     verbose = TRUE, up = 2000, dots = 200)

Arguments:
X1, X2    n x K data matrices split by type
B         total number of retained posterior samples
burnin    the burnin or warmup
verbose   if TRUE, keeps track of how many samples have been taken
up        if verbose = TRUE, prints updates ever up samples
dots      if verbose = TRUE, prints a dot every dots samples

### pbpm: Penalized Bayesian Probit Model

> pbpm(X1, X2, random = FALSE, B = 10000, burnin = B, 
>      Al = 1, Ae = 1, hc = TRUE, al = 0.01, bl = 0.01,
>      verbose = TRUE, up = 2000, dots = 200)

### bmvp: Bayesian Multivariate Probit

> bmvp(X1, X2, B = 10000, burnin = B, pen = TRUE, cholesky = FALSE,
>      Al = 1, Ae = 1, hc = FALSE, al = 0.01, bl = 0.01,
>      Kp = 2, alpha = 0.001, As = 1, Bs = 1,
>      verbose = TRUE, up = 2000, dots = 200)

### soc: data from the manuscript

> head(soc)
