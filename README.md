# Simultaneous Regression and Covariance Selection (JRNS/Stepwise)

R codes for implementing the methods proposed in the paper, **"A generalized likelihood based Bayesian approach for scalable joint regression and covariance estimation in high 
dimensions"** by Samanta, Khare and Michailidis for joint regression and covariance selection are presented here.

The functions listed here are:  
1. `JRNS` is the function that implements the JRNS algorithm and it takes the data matrices `X`, `Y` and the initial values of regression coefficient matrix, `B` and inverse covariance matrix, `$\Omega$` as inputs and returns the estimates of `B` and `$\Omega$`.  
2. `StepB` is the function which implements the first step of the Stepwise algorithm.  
3. `StepOmega` ids the function which implements the second step of the Stepwise algorithm.  

### Illustration to use the function

Let us first see how one can generate input data for the functions:

```R
#### True value of the parameters are:
n = 100
p = 200
q = 200
rho = 0.7

sigmax=matrix(NA,p,p)
for(i in 1:p) {
  for(j in 1:p){
    sigmax[i,j]=rho**(abs(i-j))
  }
}

B0 <- matrix(sample(c(runif(floor(p*0.5), 1,2), rep(0, times = p*q - floor(p*0.5)))), 
                 nrow = p, ncol = q)

ofdiag=rep(0,q*(q-1)/2)
index=sample(seq(1,qq),(q*0.2/2),replace=FALSE)
for(i in index){
  ofdiag[i]=runif(1,0.5,1)*((-1)**(rbinom(1,1,0.5)))
}

sigma_inv0[lower.tri(sigma_inv0)]=ofdiag

sigma_inv0[upper.tri(sigma_inv0)] = t(sigma_inv0)[upper.tri(sigma_inv0)]
diag(sigma_inv0) = runif(q,1,2)
var_mod0 = solve(sigma_inv0)

#### Generating the X and Y matrices 

X = mvtnorm::rmvnorm(n, sigma=sigmax)
epsilon = mvtnorm::rmvnorm(n, mean=rep(0,q),sigma=var_mod0)
Y = X %*% B0 + epsilon

```

We then need to set some initial values of the parameters to start the Gibbs Sampler:

```R
#### INITIAL VALUE OF THE PARAMETERS ##### 

  B_lasso=matrix(0,p,q)
  for(i in 1:q){
    cvglm=cv.glmnet(X,Y[,i])
    
    mod_lasso=glmnet(X,Y[,i],alpha=1,family="gaussian",lambda=cvglm$lambda.1se,intercept=F,standardize=F)
    B_lasso[,i]=as.vector(mod_lasso$beta)
  }
  
  
  er=Y-X%*%(B_lasso)
  s0=(1/n)*t(er)%*%er
  mod_glasso=glasso(s=s0,rho=sqrt(log(p)/n))
  sigma_inv_lasso=mod_glasso$wi
```

Finally, we can find the estimates by JRNS and Stepwise algorithms using the following lines of code:

```R
fit_JRNS = JRNS(B_ini=B_lasso,Omega_ini=sigma_inv_lasso,Y=Y,X=X,nmc=1000,burnin=2000)
fit_StepB = StepB(B_lasso,Y,X,ng,nburn)
fit_StepOmega = StepOmega(StepOmega(S=S_hat,temp=sigma_inv_lasso,nmc=1000,burnin=2000,n=n))
```

A few variants to each of these functions based on hyperparameter selection:  
1. `JRNS_q1q2hp.cpp` contains the `JRNShp` funtion which implements the JRNS algorithm by considering uniform hyperperiors on the prior mixture weights for `B` and `$\Omega$`.  
2. `StepB_q1hp.cpp` contains the `StepBhp` function which implements the first step of the Stepwise algorithm by considering a uniform hyperperior on the prior mixture weight.  
3. `StepOmega_q2hp.cpp` contains the `StepOmegahp` function which implements the second step of the Stepwise algorithm by considering a uniform hyperperior on the prior mixture weight.  
4. `JRNS_fixedtau1.cpp` contains the `JointRegfixed` function which implements the JRNS algorithm but sets the slab variance of the prior of each entry of `B` fixed at some value.  
5. `StepBfixed.cpp` contains the `StepBfixed` function implements the first step of the Stepwise algorithm by fixing the the slab variance of the prior of each entry of `B` to some value.


