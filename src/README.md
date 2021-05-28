# Simultaneous Regression and Covariance Selection (JRNS/Stepwise)

### Rcpp files

The 3 main Rcpp files are `JRNS`, `StepB` and `StepOmega` and their respective functions in them are described as follows:  
- `JRNS` is the function that implements the JRNS algorithm and it takes the data matrices `X`, `Y` and the initial values of regression coefficient matrix, `B` and inverse covariance matrix, `$\Omega$` as inputs and returns the estimates of `B` and `$\Omega$`.  
- `StepB` is the function which implements the first step of the Stepwise algorithm.  
- `StepOmega` ids the function which implements the second step of the Stepwise algorithm.  

A few variants to each of these functions based on hyperparameter selection:  
- `JRNS_q1q2hp.cpp` contains the `JRNShp` funtion which implements the JRNS algorithm by considering uniform hyperperiors on the prior mixture weights for `B` and `$\Omega$`.  
- `StepB_q1hp.cpp` contains the `StepBhp` function which implements the first step of the Stepwise algorithm by considering a uniform hyperperior on the prior mixture weight.  
- `StepOmega_q2hp.cpp` contains the `StepOmegahp` function which implements the second step of the Stepwise algorithm by considering a uniform hyperperior on the prior mixture weight.  
- `JRNS_fixedtau1.cpp` contains the `JointRegfixed` function which implements the JRNS algorithm but sets the slab variance of the prior of each entry of `B` fixed at some value.  
- `StepBfixed.cpp` contains the `StepBfixed` function implements the first step of the Stepwise algorithm by fixing the the slab variance of the prior of each entry of `B` to some value.


