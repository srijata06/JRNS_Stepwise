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
List JRNS(mat B_ini, mat Omega_ini, mat Y, mat X, int nmc, int burnin, double r1 = 1e-4, double s1 = 1e-8){
	
	int n = Y.n_rows;
	int p = B_ini.n_rows;
	int q = B_ini.n_cols;
	mat XtX = X.t()*X;
	mat XtY = X.t()*Y;
	mat B = zeros(p,q);
	B = B_ini ;
	mat M2 = B.t()*X.t()*X;
	mat Omega = zeros(q,q);
	Omega = Omega_ini;
	mat Omega2 = Omega*Omega;
	mat Bhat = zeros(p,q);
	mat Omegahat = zeros(q,q);
	
	
	double p1 = (double)p;
    double q1 = 1-(1/p1);
	
	vec resid =zeros(p);
	
	double b ;
	double c1;
	double c2 ;
	double a = 0;
    double d2 = 0;
    double lambda = 0;
	double eta = 0;
    double q2 = 1-(1/(double)q);  
	double new1;
	
	int components = 1;
	int rand = 1;
	double bd = 0;
	vec P_B = zeros(2);
	vec var_mod = zeros(q);
	cube Bcount = zeros(p,q,2);
	mat Bsum = zeros(p,q);
	cube count2 = zeros(q,q,2);
	mat sum2 = zeros(q,q);
	cube out1 = zeros(p,q,nmc+1);
	cube out2 = zeros(q,q,nmc+1);
	cube B_avg = zeros(p,q,nmc+burnin+1);
    vec myvec = zeros(2);
    vec P = zeros(2);
    
	clock_t time_start, time_end;
    time_start = clock();
    
    for(int t=0; t<2;t++){
		myvec(t) = t;
	}
    
	for(int it = 0; it < (nmc + burnin+1); it++){
		
		
		mat M1 = XtY*Omega2;
		
		for(int r=0; r<p; r++){
		
		  for(int s=0; s<q; s++){	
		  
		    b = XtX(r,r)*Omega2(s,s);
		  	c2 = M1(r,s) - sum(Omega2.col(s) % M2.col(r)) + B(r,s)*XtX(r,r)*Omega2(s,s);
		  	if((B(r,s) == 0)){
				eta = R::rgamma(r1, s1);   
			   }
			   else{
				eta = R::rgamma(r1 + 0.5, 1/(0.5*B(r,s)*B(r,s) + s1));   
			   }
               
			   if(eta == 0) {
				   P_B(1) = 0;
			   }
			   else {
				   P_B(1) = sqrt(eta/(b + eta)) * exp(c2*c2/(2*(b + eta))) *(1-q1)/q1; 
			   }
            	P_B(0) = 1;
			
	    	c1 = b + eta;
		if(P_B.has_inf()==true){
			   	
                int components = 1;
                new1 = R::rnorm(c2/c1, sqrt(1/c1));
                  
                if((it > burnin)){  	
			Bcount(r,s,components) += 1;
                Bsum(r,s) +=  new1;
                
                }
				
            } else {
            	
            	
            	P_B = P_B/sum(P_B);
            	components = (Rcpp::RcppArmadillo::sample(myvec,1,false, P_B) )(0);
            	if(components == 1){
			        new1 = R::rnorm(c2/c1, sqrt(1/c1));
    		        if((it > burnin))
    		         Bsum(r,s) +=  new1;
				 
				} else {
          		 
				 new1 = 0;
          		}
          		    
            	if((it > burnin))
            	Bcount(r,s,components) += 1;
            	
            	
                }
		  	M2.row(s) = M2.row(s) + (new1-B(r,s)) * XtX.row(r);
			B(r,s) = new1;
			if(it > 0){
				
				B_avg(r,s,it) =  (it*B_avg(r,s,it-1) + new1)/(it+1);
			}else{
				B_avg(r,s,it) = new1;
				
			}
		    }
    	}
		
    	if(it >= burnin)
		  out1.slice(it - burnin) = B;
	    mat E = Y-X*B;
        mat S =  E.t()*E;
		
		for(int j=0; j<(q-1); j++){
            for(int k=(j+1); k<q; k++){
            	
               a = S(j,j) + S(k,k);
               d2 = - sum(Omega.col(j)%S.col(k) + Omega.col(k)%S.col(j)) + Omega(j,k)*a;
               
			   
			   if((Omega(j,k) == 0)){
				lambda = R::rgamma(r1, s1);   
			   }
			   else{
				lambda = R::rgamma(r1 + 0.5, 1/(0.5*Omega(j,k)*Omega(j,k) + s1));   
			   }
               
			   if(lambda == 0) {
				   P(1) = 0;
			   }
			   else {
				   P(1) = sqrt(lambda/(a+lambda)) * exp(d2*d2/(2*(a+lambda))) *(1-q2)/q2; 
			   }
               P(0) = 1;
               
			   
					
               if(P.has_inf()==true){
                 int rand = 1;
                 Omega(j,k) = Omega(k,j) = R::rnorm(d2/(a+lambda), 1/sqrt(a + lambda));
				 
				 
          
                 if((it > burnin)){
				 
                     count2(j,k,rand) += 1;
                     sum2(j,k) +=  Omega(j,k);
					 
					 
                  }
                  
                  
                }  else {
          
          
                 P = P/sum(P);
          
                 rand = (Rcpp::RcppArmadillo::sample(myvec,1,false, P))(0);
          
          
                 if(rand == 1){
                 	
                     Omega(j,k) = Omega(k,j) = R::rnorm(d2/(a+lambda), 1/sqrt(a + lambda));
					 
					  
					 if(it > burnin){
					 sum2(j,k) +=  Omega(j,k);
					 }
                 } else {
                 	
                     Omega(j,k) = Omega(k,j) = 0;
					 
				    }
     
     
                  if((it > burnin))
                     count2(j,k,rand) += 1; 
            
                }
                
        
           }
		   
		   
           bd = sum(Omega.col(j)%S.col(j)) - Omega(j,j)*S(j,j) + (R::rgamma(r1 + 1, 1/(abs(Omega(j,j)) + s1)));
           Omega(j,j) = (sqrt(bd*bd + 4*S(j,j)*(double)n) - bd)/(2*S(j,j));
		   if(it > burnin)
		    sum2(j,j) += Omega(j,j);
		   
        }
        bd = sum(Omega.col(q-1)%S.col(q-1)) - Omega(q-1,q-1)*S(q-1,q-1) + (R::rgamma(r1 + 1, 1/(abs(Omega(q-1, q-1)) + s1)));
        Omega(q-1,q-1) = (sqrt(bd*bd + 4*S(q-1,q-1)*(double)n) - bd)/(2*S(q-1,q-1));
		if(it > burnin)
          sum2(q-1,q-1) += Omega(q-1,q-1);
		
		Omega2 = Omega*Omega;
		
	if(it - floor(it/1000)*1000 == 0){
      vec iteration = ones(1)*it;
      iteration.print();
    }	
    if(it >=burnin)
    out2.slice(it - burnin) = Omega; 	
            	         	
    }
	
	for (int mi = 0 ; mi < p ; mi++){
       for(int mj = 0 ; mj < q ; mj++){
       	
          if(Bcount(mi,mj,0)> (nmc/2)){
        
             Bhat(mi,mj) = 0;
        
           }else{
              Bhat(mi,mj) = Bsum(mi,mj)/Bcount(mi,mj,1);
        
           }
        }
    }
	double dnmc = (double)nmc;
  
   for(int ni = 0; ni < (q-1); ni++){
    for(int nj = (ni+1); nj < q; nj++){
    	
      if(count2(ni,nj,0) > (nmc/2)){
      	
        Omegahat(ni,nj) = 0;
        Omegahat(nj,ni) = 0;
      }else{
        Omegahat(ni,nj)  = sum2(ni,nj)/count2(ni,nj,1);
        Omegahat(nj,ni) = Omegahat(ni,nj) ;
	  }
    }
    Omegahat(ni,ni) = sum2(ni,ni)/dnmc;
  }
    Omegahat(q-1,q-1) = sum2(q-1,q-1)/dnmc; 
    
	time_end = clock();
    double time_taken;
    time_taken = (time_end - time_start) / (double)CLOCKS_PER_SEC;
    return(List::create(Named("Bcount") = Bcount,Named("Bsum") = Bsum, Named("B") = out1, Named("Bcumavg") = B_avg, Named("Omega") = out2, Named("Omegacount") = count2, Named("Omegasum") = sum2, Named("Bhat") = Bhat, Named("Omegahat") = Omegahat, Named("time") = time_taken));
       
}
	
