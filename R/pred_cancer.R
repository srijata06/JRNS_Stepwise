## finding prediction errors for different methods using TCGA cancer data 
library(gglasso)
library(glasso)
library(glmnet)
library(mvtnorm)
library(invgamma)
library(Rcpp)
library(RcppArmadillo)
library(rms)
library(dlm)

sourceCpp("JRNS_fixedtau1.cpp")
sourceCpp("StepB_fixedtau.cpp")
sourceCpp("StepOmega.cpp")

##### The following two source codes have been taken from
##### https://github.com/skdeshpande91/multivariate_SSL

sourceCpp("mSSL_dpe_source_code.cpp")
sourceCpp("mSSL_dcpe_source_code.cpp")

##### The following datasets have been taken from the supplementary of the paper
##### Bayesian structured learning in multi-layered genomic networks availabe here at
##### https://www.tandfonline.com/doi/suppl/10.1080/01621459.2020.1775611?scroll=top

##### Load the relevent data. Here the mRNA and RPPA data for just the LUAD cancer 
##### have been loaded.

load("mRNA_LUAD.dat")
load("RPPA_LUAD.dat")

mRNA = scale(mRNA)
RPPA = scale(RPPA)

nburn=1000
ng=2000
n = dim(mRNA)[1]
p = dim(mRNA)[2]
q = dim(RPPA)[2]

joint_error = rep(NA,10)
step_error = rep(NA,10)
DPE_error = rep(NA,10)
DCPE_error = rep(NA,10)
lm_error = rep(NA,10)


#### DATA ######

#Randomly shuffle the data
rand = sample(1:dim(mRNA)[1],n,replace = FALSE)
Y = RPPA[rand,]
Y = as.matrix(Y)
X = mRNA[rand,]


#Create 5 equally size folds
folds <- cut(seq(1,nrow(X)),breaks=5,labels=FALSE)

#Perform 5 fold cross validation
for(i in 1:5){
  if (i %% 1 == 0) print(paste0("Fold no.: ", i))
  #Segement your data by fold using the which() function 
  testIndices <- which(folds==i,arr.ind=TRUE)
  test.X <- X[testIndices, ]
  test.Y <- Y[testIndices, ]
  train.X <- X[-testIndices, ]
  train.Y <- Y[-testIndices, ]
  full_data = cbind(train.X,train.Y)
  #### INITIAL VALUE OF THE PARAMETERS ##### 
  B_lasso=matrix(0,p,q)
  for(ii in 1:q){
    cvglm=cv.glmnet(train.X,train.Y[,ii])
    
    mod_lasso=glmnet(train.X,train.Y[,ii],alpha=1,family="gaussian",lambda=cvglm$lambda.1se,intercept=F,standardize=F)
    B_lasso[,ii]=as.vector(mod_lasso$beta)
  }
  
  
  er=train.Y-train.X%*%(B_lasso)
  s0=(1/n)*t(er)%*%er
  mod_glasso=glasso(s=s0,rho=sqrt(log(p)/n))
  sigma_inv_lasso=mod_glasso$wi
  
  JointResults1 = JRNSfixed(B_lasso,sigma_inv_lasso,train.Y,train.X,ng,nburn)
  B_hat = JointResults1$Bhat2
  Omega_hat = JointResults1$Omegahat
  Y.est = test.X %*% B_hat
  joint_error[i] = norm(Y.est-test.Y, type = "F")
  
  StepBResults = StepBfixed(B_lasso,train.Y,train.X,ng,nburn)
  B_hat = StepBResults$Bhat1
  S_hat = (1/n) * t(train.Y - train.X%*%B_hat) %*% (train.Y - train.X %*% B_hat)
  StepOResults = StepOmega(S=S_hat,temp=sigma_inv_lasso,ng,nburn,n)
  Omega_hat = StepOResults$Omegahat
  Y.est = test.X %*% B_hat
  step_error[i] = norm(Y.est-test.Y, type = "F")
  
  DCPEResults = mSSL_dcpe(train.X, train.Y,
                          lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                          xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                          theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                          eta_hyper_params = c(1, ncol(Y)),
                          diag_penalty = 1,
                          max_iter = 500,
                          eps = 1e-3,
                          verbose = 5)
  
  B_hat = DCPEResults$B
  Omega_hat = DCPEResults$Omega
  Y.est = test.X %*% B_hat
  DCPE_error[i] = norm(Y.est-test.Y, type = "F")
  
  DPEResults = mSSL_dpe(train.X, train.Y,
                        lambdas = list(lambda1 = 1, lambda0 = seq(10, nrow(X), length = 10)),
                        xis = list(xi1 = 0.01 * nrow(X), xi0 = seq(0.1 * nrow(X), nrow(X), length = 10)),
                        theta_hyper_params = c(1, ncol(X) * ncol(Y)),
                        eta_hyper_params = c(1, ncol(Y)),
                        diag_penalty = 1,
                        max_iter = 500,
                        eps = 1e-3,
                        s_max_condition = 10*nrow(X),
                        obj_counter_max = 5,
                        verbose = 0)
  
  B_hat = DPEResults$B
  Omega_hat = DPEResults$Omega
  Y.est = test.X %*% B_hat
  DPE_error[i] = norm(Y.est-test.Y, type = "F")
  
  Y.est.lm = matrix(NA,dim(test.Y)[1],q)
  for(j in 1: q){
    lin.mod = lm(train.Y[,j]~0 + train.X)
    Y.est.lm[,j] = test.X %*% lin.mod$coefficients 
  }
  
  lm_error[i] = norm(Y.est.lm-test.Y, type = "F")
  
  
  
}

joint_error = mean(joint_error,na.rm = T)
step_error = mean(step_error,na.rm = T)
DPE_error = mean(DPE_error,na.rm = T)
DCPE_error = mean(DCPE_error,na.rm = T)
lm_error = mean(lm_error,na.rm = T)
BANS_error = mean(BANS_error,na.rm = T)

pred_error = c(joint_error,step_error,DPE_error,DCPE_error,lm_error)
pred_error = pred_error/lm_error
pred_error

