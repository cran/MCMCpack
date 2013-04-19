//////////////////////////////////////////////////////////////////////////
// MCMChierBetaBinom.cc is C++ code to fit a hierarchical beta binomial 
// model
//
//  y_{ij} ~ Binomial(s_{ij}, theta_{ij})
// theta_{ij} ~ Beta(alpha_j, beta_j)
// alpha_j ~ Pareto(1, a)
// beta_j  ~ Pareto(1, b)
//
//
// uses adaptive Metropolis scheme similar to that of Haario et al. (2001)
//   as implemented by Roberts and Rosenthal (2006/2008) 
//   "Examples of Adaptive MCMC"
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
//  5/28/2011 KQ
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////





#ifndef MCMCHIERBETABINOM_CC
#define MCMCHIERBETABINOM_CC

#include <iostream>

#include "matrix.h"
#include "algorithm.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

typedef Matrix<double,Row,View> rmview;

using namespace std;
using namespace scythe;


// used to access Ydata like a 2d array 
#define M(ROW,COL,NROWS) (COL*NROWS+ROW)




// log of the pareto density
double logdpareto(const double& x, const double& xm, const double& a){
  double logfunval;
  if (x > xm && a > 0){
    logfunval = log(a) + a*log(xm) - (a+1)*log(x);
  }
  else{
    logfunval = -numeric_limits<double>::infinity();
  }
  return logfunval;
}




// log of the full conditional density for alpha_j, beta_j
double logABfcd(const double& alpha, const double& beta, 
		   const vector<const double*>& theta, 
		   const double& a, const double& b){
  
  double term1 = 0.0;
  double term2 = 0.0;
  const int len_theta = theta.size();
  if (alpha > 1.0 && beta > 1.0){
    for (int i=0; i<len_theta; ++i){
      term1 += lndbeta1(*theta[i], alpha, beta); 
    }      
  }
  else{
    term1 = -numeric_limits<double>::infinity();
  }

  // a and/or b <= 0 is treated as improper uniform prior 
  if (a > 0){
    term2 += logdpareto(alpha, 1.0, a);
  }
  if (b > 0){
    term2 += logdpareto(beta, 1.0, b);
  }
  double logfcd = term1 + term2;

  return logfcd;

}


// candidate generating function
template <typename RNGTYPE>
Matrix<double> mixcangen(const double& alpha, const double& beta, 
			 const double& mixfrac, const double& base_sigma, 
			 const double& alpha_var_n, const double& beta_var_n, 
			 const double& ab_cov_n, 
			 rng<RNGTYPE>& stream){
  
  Matrix<double> ab_V(2, 2, false);
  ab_V = alpha_var_n, ab_cov_n, ab_cov_n, beta_var_n;
  ab_V = (5.6644 * ab_V) / 2.0;

  Matrix<double> ab_mean(2, 1, false);
  ab_mean = alpha, beta;

  Matrix<double> ab_can(2, 1, false);
  double u = stream.runif();
  if (u < mixfrac){
    double alpha_can = stream.rnorm(alpha, base_sigma);
    double beta_can = stream.rnorm(beta, base_sigma);
    ab_can = alpha_can, beta_can;
  }
  else{
    ab_can = stream.rmvnorm(ab_mean, ab_V);
  }

  return ab_can;

}




template <typename RNGTYPE>
void hierBetaBinom_impl(rng<RNGTYPE>& stream, 
			double* sampledata, const int samplerow, 
			const int samplecol,
			const int* y,
			const int* s, 
			const double* theta_start,
			const double* alpha_start,
			const double* beta_start,
			const double a,
			const double b,
			const int* ilabels,
			const int* jlabels,
			const int* ilabelsunique,
			const int* jlabelsunique,
			const int n,
			const int ni,
			const int nj,
			const int burnin,
			const int mcmc,  
			const int thin,
			const int verbose,
			int * accepts, 
			const double* base_sigma){
  
  const int tot_iter = burnin + mcmc;
// JHP  const int nstore = mcmc/thin;

  // these probably should not be hard coded
  const double mixfrac = 0.05;
  


  double* theta;
  theta = new double[n];
  for (int i=0; i<n; ++i){
    theta[i] = theta_start[i];
  }

  double* alpha;
  alpha = new double[nj];
  for (int i=0; i<nj; ++i){
    alpha[i] = alpha_start[i];
  }

  double* beta;
  beta = new double[nj];
  for (int i=0; i<nj; ++i){
    beta[i] = beta_start[i];
  }

  double* alpha_storage;
  alpha_storage = new double[nj * tot_iter];
  for (int i=0; i<(nj*tot_iter); ++i){
    alpha_storage[i] = -999.0;
  }

  double* beta_storage;
  beta_storage = new double[nj * tot_iter];
  for (int i=0; i<(nj*tot_iter); ++i){
    beta_storage[i] = -999.0;
  }
  
  double* sums_alpha_n;
  sums_alpha_n = new double[nj];
  for (int i=0; i<nj; ++i){
    sums_alpha_n[i] = 0.0;
  }
  
  double* sums_beta_n;
  sums_beta_n = new double[nj];
  for (int i=0; i<nj; ++i){
    sums_beta_n[i] = 0.0;
  }
  

  


  
  vector< vector<const double*> > js_thetas_ptr;
  js_thetas_ptr.reserve(nj);
  for (int j=0; j<nj; ++j){
    vector<const double*> holder;
    holder.reserve(n);
    for (int i=0; i<n; ++i){
      if (jlabels[i] == (j+1)){
	holder.push_back( (const double*)&theta[i] );
      }
    }
    js_thetas_ptr.push_back(holder);
  }



  int count = 0; // counter for output storage
  // MCMC iterations
  for (int iter=0; iter<tot_iter; ++iter){
   
    // sample [theta_{ij} | alpha_j, beta_j, y, s]
    for (int i=0; i<n; ++i){
      int cluster_id = jlabels[i];
      theta[i] = stream.rbeta(static_cast<double>(y[i]) + alpha[cluster_id-1],
			      static_cast<double>(s[i] - y[i]) + 
			      beta[cluster_id - 1]);
    }







    // sample [alpha_j, beta_j | theta, y, s, a, b]
    for (int j=0; j<nj; ++j){
// JHP     const int len_thetaj = js_thetas_ptr[j].size();
// JHP      double logfcd_cur = logABfcd(alpha[j], beta[j], 
// JHP			      js_thetas_ptr[j], a, b);
      if (iter < 10){
	double alpha_can = stream.rnorm(alpha[j], base_sigma[j]);
	double beta_can = stream.rnorm(beta[j], base_sigma[j]);
	double accept_numer;
	double accept_denom;
	accept_numer = logABfcd(alpha_can, beta_can, 
				   js_thetas_ptr[j], a, b);
	accept_denom = logABfcd(alpha[j], beta[j], 
				   js_thetas_ptr[j], a, b);	  
	const double ratio = exp(accept_numer - accept_denom); 
	
	if (stream.runif() < ratio) {
	  alpha[j] = alpha_can;
	  beta[j] = beta_can;
	  ++accepts[j];
	}	
      } // end iter < 10
      else{ // adaptive MCMC after the first 10 iterations
	// use iter and not iter -1 in 2 lines below b/c of start w/ 0  
	double alpha_mean_n = sums_alpha_n[j] / static_cast<double>(iter);
	double beta_mean_n = sums_beta_n[j] / static_cast<double>(iter);
	double alpha_var_n = 0.0;
	double beta_var_n = 0.0;
	double ab_cov_n = 0.0;
	for (int i=0; i<iter; ++i){
	  const double alpha_dev = (alpha_mean_n - alpha_storage[M(i, j, tot_iter)]);
	  const double beta_dev = (beta_mean_n - beta_storage[M(i, j, tot_iter)]);
	  alpha_var_n += std::pow(alpha_dev, 2);
	  beta_var_n += std::pow(beta_dev, 2);
	  ab_cov_n += alpha_dev * beta_dev;
	}
	alpha_var_n = alpha_var_n / static_cast<double>(iter);
	beta_var_n = beta_var_n / static_cast<double>(iter);
	ab_cov_n = ab_cov_n / static_cast<double>(iter);
	if (alpha_var_n <= 0.0){
	  alpha_var_n = std::pow(base_sigma[j], 2);
	  ab_cov_n = 0.0;
	}
	if (beta_var_n <= 0.0){
	  beta_var_n = std::pow(base_sigma[j], 2);
	  ab_cov_n = 0.0;
	}

	Matrix<double> ab_can = mixcangen(alpha[j], beta[j], mixfrac, 
					  base_sigma[j], alpha_var_n, 
					  beta_var_n, ab_cov_n, stream);
	double alpha_can = ab_can(0);
	double beta_can = ab_can(1);
	double accept_numer;
	double accept_denom;
	accept_numer = logABfcd(alpha_can, beta_can, 
				   js_thetas_ptr[j], a, b);
	accept_denom = logABfcd(alpha[j], beta[j], 
				   js_thetas_ptr[j], a, b);	  
	const double ratio = exp(accept_numer - accept_denom); 
	
	if (stream.runif() < ratio) {
	  alpha[j] = alpha_can;
	  beta[j] = beta_can;
	  ++accepts[j];
	}	
	
      } // end else iter not less than 10
      sums_alpha_n[j] += alpha[j];
      sums_beta_n[j] += beta[j];
      alpha_storage[M(iter, j, tot_iter)] = alpha[j];
      beta_storage[M(iter, j, tot_iter)] = beta[j];
    } // end j loop
    








    // store values 
    if (iter >= burnin && (iter % thin == 0)) {
      for (int i=0; i<n; ++i){
	sampledata[M(count, i, samplerow)] = theta[i];
      }
      for (int i=0; i<nj; ++i){
	sampledata[M(count, (i+n), samplerow)] = alpha[i];	
      }
      for (int i=0; i<nj; ++i){
	sampledata[M(count, (i+n+nj), samplerow)] = beta[i];	
      }
      ++count;
    }







    // print output to screen
    if(verbose > 0 && iter % verbose == 0){
      Rprintf("\n\niteration %i of %i \n", (iter+1), tot_iter);  
    }
    







    // allow user interrupts
    R_CheckUserInterrupt(); 
    
  } // end MCMC iterations





  // clear memory
  delete [] theta;
  delete [] alpha;
  delete [] beta;
  delete [] alpha_storage;
  delete [] beta_storage;
  delete [] sums_alpha_n;
  delete [] sums_beta_n;
  

} // end hierBetaBinom_impl













extern "C"{

  // function called by R to fit model
  void
  hierBetaBinom(double* sampledata, 
		const int* samplerow,
		const int* samplecol,
		const int* y,
		const int* s,
		const double* theta_start,
		const double* alpha_start,
		const double* beta_start,
		const double* a,
		const double* b,
		const int* ilabels,
		const int* jlabels,
		const int* ilabelsunique,
		const int* jlabelsunique,
		const int* n,
		const int* ni,
		const int* nj,
		const int* burnin, const int* mcmc,  const int* thin,
		const int *uselecuyer, 
		const int *seedarray,
		const int *lecuyerstream, 
		const int* verbose, 
		int *accepts, 
		const double* base_sigma){


    MCMCPACK_PASSRNG2MODEL(hierBetaBinom_impl, 
			   sampledata, 
			   *samplerow,
			   *samplecol,
			   y,
			   s,
			   theta_start,
			   alpha_start,
			   beta_start,
			   *a,
			   *b,
			   ilabels,
			   jlabels,
			   ilabelsunique,
			   jlabelsunique,
			   *n,
			   *ni,
			   *nj,
			   *burnin, *mcmc,  *thin,
			   *verbose, 
			   accepts, 
			   base_sigma);

  } // end hierBetaBinom
  
} // end extern "C"



#endif
