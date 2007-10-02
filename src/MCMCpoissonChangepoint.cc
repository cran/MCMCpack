// MCMCpoissonChange.cc is C++ code to estimate a Poisson changepoint model
// with a gamma prior
//
// Jong Hee Park
// Dept. of Political Science
// University of Chicago
// jhp@uchicago.edu
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Copyright (C) 2004 Andrew D. Martin and Kevin M. Quinn

#ifndef MCMCPOISSONCHANGEPOINT_CC
#define MCMCPOISSONCHANGEPOINT_CC

#include "MCMCrng.h"
#include "MCMCfcds.h"
#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "rng.h"
#include "mersenne.h"
#include "lecuyer.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;


template <typename RNGTYPE>
Matrix<double> poisson_state_sampler(rng<RNGTYPE>& stream, const int& m,
					  const Matrix<double>& Y,
					  const Matrix<double>& lambda,
					  const Matrix<double>& P){
	const int ns = m + 1;
	const int n = Y.rows();
	Matrix<> F(n, ns);
	Matrix<> pr1(ns, 1);
	pr1[0] = 1;
	Matrix<> py(ns, 1);
	Matrix<> pstyt1(ns, 1);
	
	//
	// Forward sampling: update F matrix
	//
	for (int t=0; t<n ; ++t){
		int yt = (int) Y[t]; 
		for (int j=0; j<ns ; ++j){
			py[j]  =  dpois(yt, lambda[j]);
		} 
		if (t==0) pstyt1 = pr1;
		else { 
			pstyt1 =  ::t(F(t-1,_)*P); 
		}        
		Matrix<> unnorm_pstyt = pstyt1%py;
		Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt); // pstyt = Pr(st|Yt)
		for (int j=0; j<ns ; ++j){
			F(t,j) = pstyt(j);
		}
	}// end of F matrix filtering
	
	//
	// Backward recursions
	//
	Matrix<int> s(n, 1);                            
	Matrix<> ps(n, ns); 	
	ps(n-1,_) = F(n-1,_);                      
	s(n-1) = ns;                                
	
	Matrix<> pstyn(ns, 1);
	double pone = 0.0;
	int t = n-2;
	while (t >= 0){
		int st = s(t+1);
		Matrix<> Pst_1 = ::t(P(_,st-1)); 
		Matrix<> unnorm_pstyn = F(t,_)%Pst_1; 
		pstyn = unnorm_pstyn/sum(unnorm_pstyn);   
		if (st==1)   s(t) = 1;                  		
		else{
			pone = pstyn(st-2);
			if(stream.runif () < pone) s(t) = st-1;
			else s(t) = st;
		}
		ps(t,_) = pstyn;
		--t;
	}// end of while loop
	
	Matrix<> Sout(n, ns+1); 
	Sout(_, 0) = s(_,0);
	for (int j = 0; j<ns; ++j){
		Sout(_,j+1) = ps(_, j);
		}
		
		return Sout;
		
	} // end of state sampler

//////////////////////////////////////////// 
// MCMCpoissonChangepoint implementation.  
//////////////////////////////////////////// 
template <typename RNGTYPE>
void MCMCpoissonChangepoint_impl(rng<RNGTYPE>& stream, const Matrix<>& Y,
								 Matrix<>& lambda, Matrix<>& P, const Matrix<>& A0,
								 double m, double c0, double d0,
								 unsigned int burnin, unsigned int mcmc, unsigned int thin,
								 unsigned int verbose, bool chib, 
								 Matrix<>& lambda_store, Matrix<>& P_store, 
								 Matrix<>& ps_store, Matrix<>& s_store, 
								 double& logmarglike)
{
   // define constants and form cross-product matrices
   const unsigned int tot_iter = burnin + mcmc; //total iterations
   const unsigned int nstore = mcmc / thin; // number of draws to store
   const int n = Y.rows();
   const int ns = m + 1;                 // number of states

    //MCMC loop
    unsigned int count = 0;    
	for (int iter = 0; iter < tot_iter; ++iter){

    //////////////////////
    // 1. Sample s
    //////////////////////
	Matrix<> Sout = poisson_state_sampler(stream, m, Y, lambda, P);
	Matrix<> s = Sout(_, 0);        
    Matrix<> ps(n, ns); 
    for (int j = 0; j<ns; ++j){
        ps(_,j) = Sout(_,j+1);
        }

    //////////////////////
    // 2. Sample lambda 
    //////////////////////
    Matrix<> addY(ns, 1);
    Matrix<> addN(ns, 1);

    for (int j = 0; j<ns; ++j){
         for (int i = 0; i<n; ++i){
            if (s[i] == (j+1)) { // since j starts from 0
                addN[j] = addN[j] + 1;
                addY[j] = addY[j] + Y[i];
                }// end of if
            }
            double c1 = addY[j] + c0;
            double d1 = addN[j] + 1/ d0;       
            lambda[j] = stream.rgamma(c1, d1);    
        }
        
    //////////////////////
    // 3. Sample P
    //////////////////////        
    double shape1 = 0;
    double shape2 = 0;    
    P(ns-1, ns-1) = 1; //no jump at the last state

    for (int j =0; j< (ns-1); ++j){
        shape1 =  A0(j,j) + addN[j] - 1;  
        shape2 =  A0(j,j+1) + 1; // SS(j,j+1);        
        P(j,j) = stream.rbeta(shape1, shape2);
        P(j,j+1) = 1 - P(j,j);
        }
      
    // load draws into sample array
    if (iter >= burnin && ((iter % thin)==0)){  
        for (int i=0; i<ns; ++i)
            lambda_store(count,i) = lambda[i];
        for (int j=0; j<ns*ns; ++j)    
            P_store(count,j)= P[j];
			s_store(count,_) = s;
        for (int l=0; l<n ; ++l)           
            ps_store(l,_) = ps_store(l,_) + ps(l,_);           // add two matrices
        
        ++count; // count when (iter % *thin)==0
    
        }   // end of if(iter...)
           
    
    // print output to stdout
    if(verbose > 0 && iter % verbose == 0){
        Rprintf("\n\nMCMCpoissonChangepoint iteration %i of %i \n", (iter+1), tot_iter);
        Rprintf("lambda = \n");
        for (int j=0; j<ns; ++j)
        Rprintf("%10.5f\n", lambda[j]);
        }
       
    R_CheckUserInterrupt(); // allow user interrupts   

    }// end MCMC loop
  
if(chib ==1){
        
    Matrix<> lambda_st = meanc(lambda_store); 
    Matrix<> P_vec_st = meanc(P_store);
    const Matrix<> P_st(ns, ns);
    for (int j = 0; j< ns*ns; ++j){  
        P_st[j] = P_vec_st[j]; 
        }    

    //////////////////////
    // 1. pdf.lambda
    //////////////////////  
    Matrix<> density_lambda(nstore, ns);
    for (int iter = 0; iter<nstore; ++iter){

    Matrix<> addY(ns, 1);
    Matrix<> addN(ns, 1);
    
    for (int j = 0; j<ns; ++j){
         for (int i = 0; i<n; ++i){
            if (s_store(iter, i) == (j+1)) { 
                addN[j] = addN[j] + 1;
                addY[j] = addY[j] + Y[i];
                }// end of if
            } // end of for i
        double c1 = addY[j] + c0;
        double d1 = addN[j] + d0;   
        density_lambda(iter, j) = dgamma(lambda_st[j], c1, 1/d1);
        } // end of for j 
             
    }// end of pdf.lambda MCMC run  
    double pdf_lambda = log(prod(meanc(density_lambda)));
	
    //////////////////////
    // 2. pdf.P
    //////////////////////
    Matrix<> density_P(nstore, ns);
    for (int iter = 0; iter < nstore; ++iter){
		Matrix<> Sout = poisson_state_sampler(stream, m, Y, lambda_st, P);
        Matrix <> s = Sout(_, 0);
        Matrix <> ps(n, ns); 
        for (int j = 0; j<ns; ++j){
            ps(_,j) = Sout(_,j+1);
            }

        double shape1 = 0;
        double shape2 = 0;    
        P(ns-1, ns-1) = 1; 
        Matrix <> P_addN(ns, 1); 
        for (int j = 0; j<ns; ++j){
            for (int i = 0; i<n; ++i){
                if (s[i] == (j+1)) { // since j starts from 0
                P_addN[j] = P_addN[j] + 1;            
                    }// end of if
                } // end of for i
            } // end of for j        
        
        for (int j =0; j< (ns-1); ++j){
            shape1 =  A0(j,j) + P_addN[j] - 1;  
            shape2 =  A0(j,j+1) + 1; //         
            P(j,j) = stream.rbeta(shape1, shape2);
            P(j,j+1) = 1 - P(j,j);
            density_P(iter, j) = dbeta(P_st(j,j), shape1, shape2); 
            } // end of for j
            density_P(iter, ns-1) = 1; //
          
    }// end of pdf.P MCMC run            
    double pdf_P = log(prod(meanc(density_P)));
    
    //////////////////////
    // likelihood
    //////////////////////       
    Matrix<> F(n, ns);
    Matrix<> like(n, 1);
    Matrix<> pr1(ns, 1);
    pr1[0] = 1;
    Matrix<> py(ns, 1);
    Matrix<> pstyt1(ns, 1);

    for (int t=0; t<n ; ++t){
       int yt = (int) Y[t]; 
       for (int j=0; j<ns ; ++j){
       py[j]  =  dpois(yt, lambda_st[j]);
            } 
       if (t==0) pstyt1 = pr1;
       else { 
            pstyt1 =  ::t(F(t-1,_)*P_st); // make it an ns by 1 matrix
            }        
       Matrix<double> unnorm_pstyt = pstyt1%py;
       Matrix<double> pstyt = unnorm_pstyt/sum(unnorm_pstyt); // pstyt = Pr(st|Yt)
       for (int j=0; j<ns ; ++j){
        F(t,j) = pstyt(j);
            }
       like[t] = sum(unnorm_pstyt);
        }// end of likelihood computation
    
    double loglike = sum(log(like));
	
	Rprintf("loglike is \n");
	Rprintf("%10.5f\n", loglike); 
    
    //////////////////////
    // log prior ordinate 
    //////////////////////
    Matrix<> density_lambda_prior(ns, 1);
    Matrix<> density_P_prior(ns, 1);
    density_P[ns-1] = 1; //
    
    for (int j=0; j<ns ; ++j){
        density_lambda_prior[j] = log(dgamma(lambda_st[j], c0, d0));    
        }   
        
    for (int j =0; j< (ns-1); ++j){
        density_P_prior[j] = log(dbeta(P_st(j,j), A0(j,j), A0(j,j+1))); 
        }        
    
    // compute marginal likelihood
    double logprior = sum(density_lambda_prior) + sum(density_P_prior);
    logmarglike = (loglike + logprior) - (pdf_lambda + pdf_P);
	Rprintf("logmarglike is \n");
	Rprintf("%10.5f\n", logmarglike); 
	Rprintf("-------------------------------------------\n");
	
        
}// end of marginal likelihood computation

}
//////////////////////////////////////////// 
// Start MCMCpoissonChangepoint function
///////////////////////////////////////////
extern "C"{
    void MCMCpoissonChangepoint(double *lambdaout, double *Pout, double *psout, double *sout,  
           const double *Ydata, const int *Yrow, const int *Ycol, const int *m, 
		   const int *burnin, const int *mcmc, const int *thin, const int *uselecuyer, const int *seedarray,
		   const int *lecuyerstream, const int *verbose, 
           const double *lambdastart, const double *Pstart, 
           const double *a, const double *b, const double *c0, const double *d0,  
           const double *A0data, double *logmarglikeholder, 
           const int *chib){
   
    // pull together Matrix objects
    const Matrix <double> Y(*Yrow, *Ycol, Ydata);  
    const unsigned int tot_iter = *burnin + *mcmc; //total iterations
    const unsigned int nstore = *mcmc / *thin; // number of draws to store
    const int n = Y.rows();
    const int ns = *m + 1;                 // number of states
 	
   // generate starting values
    Matrix <> lambda(ns, 1, lambdastart);
    const Matrix <> A0(ns, ns, A0data);
    Matrix <> P(ns, ns, Pstart);

    double logmarglike;
	// double loglike;

    // storage matrices
    Matrix<> lambda_store(nstore, ns);
    Matrix<> P_store(nstore, ns*ns);
    Matrix<> ps_store(n, ns);
    Matrix<> s_store(nstore, n);

    MCMCPACK_PASSRNG2MODEL(MCMCpoissonChangepoint_impl, Y,
						   lambda, P, A0, *m, *c0, *d0,  *burnin, *mcmc, *thin, *verbose, *chib, 
						   lambda_store, P_store, ps_store, s_store, logmarglike)
        
	logmarglikeholder[0] = logmarglike;
	// loglikeholder[0] = loglike;    
    
    // return output
    for (int i = 0; i<(nstore*ns); ++i){
        lambdaout[i] = lambda_store[i]; 
        }
    for (int i = 0; i<(nstore*ns*ns); ++i){
        Pout[i] = P_store[i]; 
        }
    for (int i = 0; i<(n*ns); ++i){
        psout[i] = ps_store[i]; 
        }
    for (int i = 0; i<(nstore*n); ++i){
        sout[i] = s_store[i];
        }
 
 }// end of MCMCpoissonChange
}//end extern "C"

#endif
