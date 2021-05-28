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
List StepBfixed(mat B_ini, mat Y, mat X, int nmc, int burnin,double r1 = 1e-4, double s1 = 1e-8){
	
	int n = Y.n_rows;
	int p = B_ini.n_rows;
	int q = B_ini.n_cols;
	mat XtX = X.t()*X;
	mat YtX = Y.t()*X;
	mat B = zeros(p,q);
	mat Bhat1 = zeros(p,q);
	mat Bhat2 = zeros(p,q);
	B = B_ini ;
	double tau = 1;
	double p1 = (double)p;
    double q1 = 1-(1/p1);
	double nb ;
	vec resid =zeros(n);
	double beta ;
	double c1 ;
	double c2 ;
	
	
	int components = 1;
	vec P_B = zeros(2);
	vec var_mod = zeros(q);
	cube Bcount = zeros(p,q,2);
    mat CI_lower = zeros(p,q);
	mat CI_upper = zeros(p,q);
	cube CI = zeros(p,q,2);
	mat Bsum = zeros(p,q);
	cube out = zeros(p,q,nmc+1);
	vec myvec = zeros(2);
	mat se = zeros(p,q);
	mat invrs = zeros(p,q);
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
            beta = s1 + 0.5*sum(resid % resid) + sum(B.col(k) % B.col(k))/(2 * tau*tau);
            var_mod(k) = R:: rgamma( r1 + (nb/2) + ((double)n/2), 1/beta);
            var_mod(k) = 1/var_mod(k);
            
           for(int j=0; j<p; j++){
           	
			   c1 = XtX(j,j) + 1/(tau*tau)	;
			   c2 = YtX(k,j) - sum(XtX.col(j)%B.col(k)) + B(j,k)*XtX(j,j);
			   P_B(1) = (1 - q1)*exp(c2*c2/(2*var_mod(k)*c1)) / (q1*tau*sqrt(c1));
			   P_B(0) = 1;
			   
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
	vec S1 = zeros(q);
	mat s(1,1);
	mat s_new(1,1);
	uvec gamma1;
	mat temp1;
	vec temp2;
	mat check1;
	mat check2;
	mat check3;
	
	for(int mj = 0 ; mj < q ; mj++) {
	  	
       gamma1 = find( Bcount0.col(mj) <= (nmc/2));
       if(gamma1.size() > 0){
       	mat I = eye(gamma1.size(),gamma1.size());
       	temp1 = inv( (X.cols(gamma1)).t() * X.cols(gamma1) + (1/(tau*tau)) * I);
		temp2 = temp1 * ((X.cols(gamma1)).t() * Y.col(mj));
       	vec Bhat1mj = zeros(p);
       	Bhat1mj(gamma1) =  temp2;
       	Bhat1.col(mj) = Bhat1mj;
		s = (1/(double)n) * ( sum(Y.col(mj) % Y.col(mj)) - (Y.col(mj).t() * (X.cols(gamma1)) * temp2) );
		check1 = sum(Y.col(mj) % Y.col(mj));
		check2 = Y.col(mj).t(); 
		check3 = X.cols(gamma1);
		s_new =  ( sum(Y.col(mj) % Y.col(mj)) - (Y.col(mj).t() * (X.cols(gamma1)) * temp2) );
		S1(mj) = s_new(0,0);
		vec ciu = zeros(p);
		vec cil = zeros(p);
		vec semj = zeros(p);
		vec invrsmj = zeros(p);
		invrsmj(gamma1) = diagvec(temp1);
		invrs.col(mj) = invrsmj;
		semj(gamma1) = sqrt((1 + n*s(0,0))* diagvec(temp1)/(n + 2));
		se.col(mj) =  semj;
		ciu(gamma1) = temp2 + (R :: qt(0.975, n + 2, TRUE, FALSE)) * semj(gamma1);
		cil(gamma1) = temp2 - (R :: qt(0.975, n + 2, TRUE, FALSE)) * semj(gamma1);
        CI_lower.col(mj) = cil;
		CI_upper.col(mj) = ciu;
		
		
       } 
    }
    
    for (int mi = 0 ; mi < p ; mi++){
       for(int mj = 0 ; mj < q ; mj++){
       	
          if(Bcount(mi,mj,0)> (nmc/2)){
        
             Bhat2(mi,mj) = 0;
        
           }else{
              Bhat2(mi,mj) = Bsum(mi,mj)/Bcount(mi,mj,1);
        
           }
        }
    }
	
	time_end = clock();
    double time_taken;
    time_taken = (time_end - time_start) / (double)CLOCKS_PER_SEC;
  
	return(List::create(Named("Bcount") = Bcount,Named("Bsum") = Bsum, Named("B") = out, Named("Bhat1") = Bhat1, Named("Bhat2") = Bhat2, Named("time") = time_taken, Named("lower_limit") = CI_lower, Named("upper_limit") = CI_upper ));
}
	

