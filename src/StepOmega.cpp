#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <ctime>

using namespace Rcpp;
using namespace arma;
using namespace std;


// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins("cpp11")]]


// [[Rcpp::export]]
List StepOmega(mat S, mat temp, int nmc, int burnin, int n, double r = 1e-4, double s = 1e-8){
	
  int p = S.n_cols;
  mat Omega = zeros(p,p);
  Omega = temp;
  mat Omegahat = zeros(p,p);
  cube out = zeros(p,p,nmc+burnin+1); // out is an array of estimated precision matrices in the sampling stage (i.e after burning stage)
  cube OmegaCount = zeros(p,p,2); // OmegaCount is an array of the votes given to zero and non-zero edges
  vec P = zeros(2);
  vec myvec = zeros(2);
  for(int t=0; t<2;t++)
    myvec(t) = t;
  
  double a = 0;
  double b = 0;
  double bd = 0;
  double lambda = 0;
  double q2 = 1-(1/(double)p);
  int rand = 1;
  mat OmegaSum = zeros(p,p);
  clock_t time_start, time_end;
  time_start = clock();
  
  

  for(int it = 0; it < (nmc + burnin+1); it++){

    for(int j=0; j<(p-1); j++){
      for(int k=(j+1); k<p; k++){
        
        a = S(j,j) + S(k,k);
        b = - sum(Omega.col(j)%S.col(k) + Omega.col(k)%S.col(j)) + Omega(j,k)*a;
		
        if((Omega(j,k) == 0)){
		 lambda = R::rgamma(r, s);   
	    }
		else{
		lambda = R::rgamma(r + 0.5, 1/(0.5*Omega(j,k)*Omega(j,k) + s));   
		}
         
        
		P(1) = sqrt(lambda/(n*a+lambda)) * exp((n*b)*(n*b)/(2*(n*a+lambda))) *(1-q2)/q2;
        P(0) = 1;
		
	
        if(P.has_inf()==true){
        	
          rand = 1;
          Omega(j,k) = Omega(k,j) = R::rnorm(n*b/(n*a+lambda), 1/sqrt(n*a + lambda));
          
          if((it > burnin)){
          	
          	OmegaCount(j,k,rand) += 1;
			OmegaCount(k,j,rand) =  OmegaCount(j,k,rand);
          	
          	OmegaSum(j,k)  += Omega(j,k);
          	OmegaSum(k,j)  += Omega(k,j);
          	
		    }
            
            
        }  else {
          
          
          P = P/sum(P);
          
          rand = (Rcpp::RcppArmadillo::sample(myvec,1,false, P))(0);
          
          
          if(rand == 1){
            Omega(j,k) = Omega(k,j) = R::rnorm(n*b/(n*a+lambda), 1/sqrt(n*a + lambda));
            
            if((it > burnin)){
            	
            OmegaSum(j,k) +=  Omega(j,k);
            OmegaSum(k,j) += Omega(k,j);
			}
            
    	      
    	       
          } else {
            Omega(j,k) = Omega(k,j) = 0;
          }
     
     
          if((it > burnin))
            OmegaCount(j,k,rand) += 1; 
		    OmegaCount(k,j,rand) =  OmegaCount(j,k,rand);
            
        }

        
      }
      
      bd = sum(Omega.col(j)%S.col(j)) - Omega(j,j)*S(j,j) + (R::rgamma(r + 1, 1/(abs(Omega(j,j)) + s)))/n;
      Omega(j,j) = (sqrt(bd*bd + 4*S(j,j)) - bd)/(2*S(j,j));
	  
      if(it > burnin)
        OmegaSum(j,j) += Omega(j,j); 
    }
    bd = sum(Omega.col(p-1)%S.col(p-1)) - Omega(p-1,p-1)*S(p-1,p-1) + R::rgamma(r + 1, 1/(abs(Omega(p-1, p-1)) + s))/n;
    Omega(p-1,p-1) = (sqrt(bd*bd + 4*S(p-1,p-1)) - bd)/(2*S(p-1,p-1));
	if(it > burnin)
      OmegaSum(p-1,p-1) += Omega(p-1,p-1);
    
    
    
   
   if(it - floor(it/1000)*1000 == 0){
      vec iteration = ones(1)*it;
      iteration.print();
    }
    if(it >=burnin)
      out.slice(it - burnin) = Omega; 
    
  }
  
  double dnmc = (double)nmc;
  
  for(int ni = 0; ni < (p-1); ni++){
    for(int nj = (ni+1); nj < p; nj++){
    	
      if((OmegaCount(ni,nj,0) > (0.5*nmc))){
      	
        Omegahat(ni,nj) = 0;
        Omegahat(nj,ni) = 0;
      }else{
        Omegahat(ni,nj)  = OmegaSum(ni,nj)/OmegaCount(ni,nj,1);
        Omegahat(nj,ni) = Omegahat(ni,nj) ;
	  }
    }
    Omegahat(ni,ni) = OmegaSum(ni,ni)/dnmc;
  }
  
  Omegahat(p-1,p-1) = OmegaSum(p-1,p-1)/dnmc; 
  
  time_end = clock();
  double time_taken;
  time_taken = (time_end - time_start) / (double)CLOCKS_PER_SEC;
  return(List::create(Named("OmegaCount") = OmegaCount, Named("OmegaSum") = OmegaSum, Named("Omega") = out, Named("Omegahat") = Omegahat, Named("time") = time_taken));
}
