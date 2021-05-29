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
List StepBhp(mat B_ini, vec sigma2_ini, mat Y, mat X, int nmc, int burnin, double r1 = 1e-4, double s1 = 1e-8,double alpha1 = 1, double beta1 = 1){
	
	int n = Y.n_rows;
	int p = B_ini.n_rows;
	int q = B_ini.n_cols;
	mat XtX = X.t()*X;
	mat YtX = Y.t()*X;
	mat B = zeros(p,q);
	mat Bhat = zeros(p,q);
	B = B_ini ;
	double tau;
	double eta;
	double q1;
	double nb ;
	vec resid =zeros(p);
	double beta ;
	double c1 ;
	double b;
	double c2 ;
	
	int components = 1;
	vec P_B = zeros(2);
	vec var_mod = sigma2_ini;
	cube Bcount = zeros(p,q,2);
        mat CI_lower = zeros(p,q);
	mat CI_upper = zeros(p,q);
	cube CI = zeros(p,q,2);
	mat Bsum = zeros(p,q);
	cube out = zeros(p,q,nmc+1);
	vec myvec = zeros(2);
	
	clock_t time_start, time_end;
    time_start = clock();
    
    for(int t=0; t<2;t++){
		myvec(t) = t;
	}
    
	for(int it = 0; it < (nmc + burnin+1); it++){
	 	 	
	    for(int k=0; k<q; k++){
	    	
	    	
	    	uvec idx = find(B.col(k) != 0);
	    	nb = idx.size();
	    	resid = Y.col(k) - X*B.col(k);
            
            
           for(int j=0; j<p; j++){
           	
			   b = XtX(j,j) ;
			   c2 = YtX(k,j) - sum(XtX.col(j)%B.col(k)) + B(j,k)*XtX(j,j);
			   
			   if((B(j,k) == 0)){
				eta = R::rgamma(r1, s1);   
				q1 = R::rbeta(alpha1+1,beta1);
			   }
			   else{
				eta = R::rgamma(r1 + 0.5, 1/((0.5*B(j,k)*B(j,k)/var_mod(k)) + s1));  
                                q1 = R::rbeta(alpha1,beta1+1);				
			   }
               
			   
			if(eta == 0) {
				   P_B(1) = 0;
			   }
			   else {
				   P_B(1) = sqrt(eta/(b + eta)) * exp(c2*c2/(2*(b + eta))) *(1-q1)/q1; 
			   }
            P_B(0) = 1;
            tau = 1/sqrt(eta);
			c1 = b + eta;
			beta =  s1 + 0.5*sum(resid % resid) + sum(B.col(k) % B.col(k))/(2 * tau*tau);
            var_mod(k) = R:: rgamma( r1 + (nb/2) + (n/2), 1/beta);
            var_mod(k) = 1/var_mod(k);
			   
			   if(P_B.has_inf()==true){
			   	
                  int components = 1;
                  B(j,k) = R::rnorm(c2/c1, sqrt(var_mod(k)/c1));
                  
                  if((it > burnin)){
                  	
				     Bcount(j,k,components) += 1;
                     Bsum(j,k) +=  B(j,k);
                     
                     }
                     
                } else {
                	
                	
                	P_B = P_B/sum(P_B);
		   
			        components = (Rcpp::RcppArmadillo::sample(myvec,1,false, P_B) )(0);
			        
			        if(components == 1){
			        	
    		          B(j,k) = R::rnorm(c2/c1, sqrt(var_mod(k)/c1));
    		          
    		          if((it > burnin))
    		          Bsum(j,k) +=  B(j,k);
    		          
          		    } else {
          		    	
            		B(j,k) = 0;
          		    }
          		    
            	    if((it > burnin))
            		Bcount(j,k,components) += 1;
            	}
            	
            	
       		}
       		
       		
   		} 
   		if(it - floor(it/1000)*1000 == 0){
        vec iteration = ones(1)*it;
        iteration.print();
        }
   		
		if(it >= burnin)
		    out.slice(it - burnin) = B	;
		    
		    
	}
	mat Bcount0 = Bcount.slice(0);
	
    
    for (int mi = 0 ; mi < p ; mi++){
       for(int mj = 0 ; mj < q ; mj++){
       	
          if(Bcount(mi,mj,0)> (nmc/2)){
        
             Bhat(mi,mj) = 0;
        
           }else{
              Bhat(mi,mj) = Bsum(mi,mj)/Bcount(mi,mj,1);
        
           }
        }
    }
	
	time_end = clock();
    double time_taken;
    time_taken = (time_end - time_start) / (double)CLOCKS_PER_SEC;
  
	return(List::create(Named("Bcount") = Bcount,Named("Bsum") = Bsum, Named("B") = out, Named("Bhat") = Bhat, Named("time") = time_taken, Named("lower_limit") = CI_lower, Named("upper_limit") = CI_upper ));
}
	

