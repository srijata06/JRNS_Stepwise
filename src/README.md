# Simultaneous Regression and Covariance Selection (JRNS/Stepwise)

### Rcpp files

The 3 main Rcpp files are `JRNS`, `StepB` and `StepOmega` and each contains a single function described as follows:  
- `JRNS` is the function that implements the JRNS algorithm and it takes the data matrices `X`, `Y` and the initial values of regression coefficient matrix, `B` and inverse covariance matrix, `$\Omega$` as inputs and returns the estimates of `B` and `$\Omega$`.  
- `StepB` is the function which implements the first step of the Stepwise algorithm.  
- `StepOmega` is the function which implements the second step of the Stepwise algorithm.  

The JRNS and StepOmega functions have an optional argument called `posdef` which when set to 1 will allow to check whether the resulting estimate of `$\Omega$` is positive definite or not and in case it is not will induce a tranformation on the estimate. The default value is 0. These functions along with the StepB function also have another argument called `hyp` which when set to 1 allows the user to specify a Beta hyper-prior for the prior mixture weights `$q_1$` and `$q_2$`. Default is set to 0.

A few variants to each of these functions based on hyperparameter selection:  
 - `JRNS_fixedtau1.cpp` contains the `JointRegfixed` function which implements the JRNS algorithm but sets the slab variance of the prior of each entry of `B` fixed at some value.  
- `StepBfixed.cpp` contains the `StepBfixed` function implements the first step of the Stepwise algorithm by fixing the the slab variance of the prior of each entry of `B` to some value.


