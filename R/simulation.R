library(gglasso)
library(glasso)
library(glmnet)
library(mvtnorm)
library(invgamma)
library(Rcpp)
library(RcppArmadillo)
library(rms)
library(dlm)

sourceCpp("JRNS.cpp")
sourceCpp("StepB.cpp")
sourceCpp("StepOmega.cpp")

##### The following code has a function which is available at
##### https://github.com/skdeshpande91/multivariate_SSL

source("Sparsity_measures.R")


n=100
nsim=50
nburn=1000
ng=2000
p=30
q=60

B_hat = array(NA,dim = c(p,q,nsim,2))
Omega_hat = array(NA,dim = c(q,q,nsim,2))

Bcount = array(NA,dim = c(p,q,2,nsim))
Omegacount = array(NA,dim = c(q,q,2,nsim))
Omegacount_BANS = array(NA,dim = c(q,q,nsim))
stv_B=matrix(NA,nsim,2)
spc_B=matrix(NA,nsim,2)
mcc_B=matrix(NA,nsim,2)
norm_B=matrix(NA,nsim,2)
cov_B_step = array(rep(NA,p*q*nsim),c(p,q,nsim))
CI_lower_step = array(rep(NA,p*q*nsim),c(p,q,nsim))
CI_upper_step = array(rep(NA,p*q*nsim),c(p,q,nsim))
cov_B_joint = array(rep(NA,p*q*nsim),c(p,q,nsim))
CI_lower_joint = array(rep(NA,p*q*nsim),c(p,q,nsim))
CI_upper_joint = array(rep(NA,p*q*nsim),c(p,q,nsim))

stv_omega=matrix(NA,nsim,2)
spc_omega=matrix(NA,nsim,2)
mcc_omega=matrix(NA,nsim,2)
norm_omega=matrix(NA,nsim,2)
time = matrix(NA,nsim,2)

#### TRUE MODEL PARAMETERS ######    

nz_B=p/5
nz_omega=q/5 
qq=q*(q-1)/2

sigma_inv0=matrix(NA,q,q)
ofdiag=rep(0,qq)
sigmax=matrix(NA,p,p)
rho=0.6

for(i in 1:p) {
  for(j in 1:p){
    sigmax[i,j]=rho**(abs(i-j))
  }
} 

B0=rep(0,p*q)
index1=sample(seq(1,q*p),nz_B,replace=FALSE)
for(i in index1){
  B0[i]=runif(1,1,2)
}
B0=matrix(B0,p,q)

ofdiag=rep(0,qq)
index2=sample(seq(1,qq),(nz_omega/2),replace=FALSE)
for(i in index2){
  ofdiag[i]=runif(1,0.5,1)*((-1)**(rbinom(1,1,0.5)))
}

sigma_inv0[lower.tri(sigma_inv0)]=ofdiag

for(i in 1:(q-1)){
  for(j in (i+1):q){
    sigma_inv0[i,j]=sigma_inv0[j,i]
  }
}

for(i in 1:q)
{
  sigma_inv0[i,i]=runif(1,1,2)
}

##### If the sigma_inv0 matrix does not turn out to be PD we do the following ####

eigen = eigen(sigma_inv0,symmetric=TRUE)$values
if(min(eigen)<0){
  diag(sigma_inv0) = diag(sigma_inv0) + abs(min(eigen)) +1
}


var_mod0 = solve(sigma_inv0)    #### sigma
X = mvtnorm::rmvnorm(n, sigma=sigmax)

#### Checking SNR. Sometimes one might need to inflate B0 so that the SNR is above 2 ###

SNR=norm(B0,type="F")/norm(var_mod0,type="F")
SNR

################ GIbbs Sampler ######

for (it in 1:nsim){
  if (it %% 1 == 0) print(paste0("Sim: ", it))
  
  #### Generating DATA ######
  
  epsilon = mvtnorm::rmvnorm(n, mean=rep(0,q),sigma=var_mod0)
  Y = X %*% B0 + epsilon
  
  full_data = cbind(X,Y)
  xtx = t(X) %*% X
  ytx = t(Y) %*% X
  
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
  
  JointResults1 = JRNS(B_lasso,sigma_inv_lasso,Y,X,ng,nburn)
  B_hat[,,it,1] = JointResults1$Bhat
  Omega_hat[,,it,1] = JointResults1$Omegahat
  Bcount[,,,it] = JointResults1$Bcount
  Omegacount[,,,it] = JointResults1$Omegacount
  time[it,1]= JointResults1$time
  eB = error_B(B_hat[,,it,1],B0)
  tpB[it,1] = eB["TP"]
  tnB[it,1] = eB["TN"]
  fpB[it,1] = eB["FP"]
  fnB[it,1] = eB["FN"]
  stv_B[it,1] = eB["SEN"]
  spc_B[it,1] = eB["SPE"]
  mcc_B[it,1] = eB["MCC"]
  norm_B[it,1] = norm(B_hat[,,it,1]-B0,type="F")/norm(B0,type="F")
  
  B_all = JointResults1$B
  B_all[B_all==0]<-NA
  CI_lower_joint[,,it] = apply(B_all,1:2,quantile,probs = 0.025,na.rm = TRUE)
  CI_upper_joint[,,it] = apply(B_all,1:2,quantile,probs = 0.975,na.rm = TRUE)
  for(ii in 1:p){
    for(jj in 1:q){
      cov_B_joint[ii,jj,it] = 1*(CI_lower_joint[ii,jj,it] <= B0[ii,jj] & CI_upper_joint[ii,jj,it] >= B0[ii,jj])
    }
  }
  eO1 = error_Omega(Omega_hat[,,it,1],sigma_inv0)
  tpS[it,1] = eO1["TP"]
  tnS[it,1] = eO1["TN"]
  fpS[it,1] = eO1["FP"]
  fnS[it,1] = eO1["FN"]
  stv_omega[it,1] = eO1["SEN"]
  spc_omega[it,1] = eO1["SPE"]
  mcc_omega[it,1] = eO1["MCC"]
  norm_omega[it,1] = norm(Omega_hat[,,it,1]-sigma_inv0,type="F")/norm(sigma_inv0,type="F")
  
  StepBResults = StepB(B_lasso,Y,X,ng,nburn)
  B_hat[,,it,2] = StepBResults$Bhat
  B_all = StepBResults$B
  B_all[B_all==0]<-NA
  CI_lower_step[,,it] = apply(B_all,1:2,quantile,probs = 0.025,na.rm = TRUE)
  CI_upper_step[,,it] = apply(B_all,1:2,quantile,probs = 0.975,na.rm = TRUE)
  for(ii in 1:p){
    for(jj in 1:q){
      cov_B1[ii,jj,it] = 1*(CI_lower1[ii,jj,it] <= B0[ii,jj] & CI_upper1[ii,jj,it] >= B0[ii,jj])
      cov_B2[ii,jj,it] = 1*(CI_lower2[ii,jj,it] <= B0[ii,jj] & CI_upper2[ii,jj,it] >= B0[ii,jj])
    }
  }
  
  S_hat = (1/n) * t(Y - X%*%B_hat[,,it,2]) %*% (Y - X %*% B_hat[,,it,2])
  StepOResults = StepOmega(S=S_hat,temp=sigma_inv_lasso,nmc=ng,burnin=nburn,n)
  Omega_hat[,,it,2] = StepOResults$Omegahat
  time[it,2]= StepBResults$time + StepOResults$time
  eB = error_B(B_hat[,,it,2],B0)
  stv_B[it,2] = eB["SEN"]
  spc_B[it,2] = eB["SPE"]
  mcc_B[it,2] = eB["MCC"]
  norm_B[it,2] = norm(B_hat[,,it,2]-B0,type="F")/norm(B0,type="F")
  
  eO2 = error_Omega(Omega_hat[,,it,2],sigma_inv0)
  stv_omega[it,2] = eO2["SEN"]
  spc_omega[it,2] = eO2["SPE"]
  mcc_omega[it,2] = eO2["MCC"]
  norm_omega[it,2] = norm(Omega_hat[,,it,2]-sigma_inv0,type="F")/norm(sigma_inv0,type="F")
  
  
}



