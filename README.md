# BMMP

Code for estimating Bayesian methods for multivariate matched proportions (MMP) as described by Meyer, Cheng, and Knutson (2023), recently accepted for publication in Statistics in Biosciences. A link to the pre-print is forth-coming.

## Folder Details

The folder titled Code contains the file bmmp.R which has all of the relevant functions to estimate our multivariate probit-based approach for Bayesian MMPs as well as non-Bayesian methods from the existing literature. This file also contains the data used in the reanalysis described in Meyer, Cheng, and Knutson (2023). The folder titled Application contains the file soc_reanalysis.R which can be used to replicate the reanalysis. The files in the folder Simulation were used to for the Simulation Study. These files can take several days to run.

## Brief Function and Object Description

All functions two matrices, X1 and X2, as arguments. X1 is an n x K matrix of the first set of binary multivariate outcomes. X2 is an n x K matrix of the second or paired set of binary multivariate outcomes. The method described in the paper is implemented by the function bmvp(). It imple Priors default to those described in the manuscript. 

### bmvp: Bayesian Multivariate Probit

> bmvp(X1, X2, B = 10000, burnin = B, pen = TRUE, cholesky = FALSE,
>      Al = 1, Ae = 1, hc = FALSE, al = 0.01, bl = 0.01,
>      Kp = 2, alpha = 0.001, As = 1, Bs = 1,
>      verbose = TRUE, up = 2000, dots = 200)

### soc: data from the manuscript

> head(soc)

Additional functions include those to implement the Bootstrap- and GEE-based methods from the exisiting literature.

### geemmp: GEE-based method, Klingenberg and Agresti (2006)

> geemmp(X1, X2, family = gaussian(link = 'identity'), corstr = 'independence')

### bootmmp: Bootstrap-based method, Westfall, Troendle, and Pennello (2010)

> bootmmp(X1, X2, B, int.hyp = FALSE)

All methods are classed with corresponding print(), summary(), and coef() functions to print out estimated differences in marginal probabilities for each set of matched proportions. summary() further provides interval estimates, hypothesis testing features, and model fit criteria.
