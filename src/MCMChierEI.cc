// fits Wakefield's hierarchical model for ecological inference using
// Wakefield's normal approximation to the binomial convolution likelihood
// and slice sampling and Gibbs sampling to sample from the posterior
//
// KQ 3/2/2002
// KQ 10/25/2002 [ported to Scythe0.3 and written for an R interface]
// KQ 7/20/2004 [minor changes regarding output and user interrupts]
// ADM 7/24/2004 [updated to new Scythe version]
// KQ 8/10/2004 bug fix and major overhaul of sampling scheme

#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts


using namespace SCYTHE;
using namespace std;

static double Lev1thetaPost(double theta[], const double& r0, const double& r1,
		     const double& c0, const double& mu0, const double& mu1,
		     const double& sigma0, const double& sigma1){
  const double theta0 = theta[0];
  const double theta1 = theta[1];
  const double p0 = 1.0/(1.0 + exp(-1*theta0));
  const double p1 = 1.0/(1.0 + exp(-1*theta1));
  const double logprior = lndnorm(theta0, mu0, sqrt(sigma0)) + 
    lndnorm(theta1, mu1, sqrt(sigma1));
  const double loglike = lndnorm(c0, r0*p0 + r1*p1,
				 sqrt(r0*p0*(1.0-p0) + 
				      r1*p1*(1.0-p1)));
  return(loglike + logprior);  
} 



// eventually all of the slice sampling functions should be made more 
// general and put in MCMCfcds.{h cc}
//
// Radford Neal's (2000) doubling procedure coded for a logdensity
static void doubling(double (*logfun)(double[], const double&, const double&,
			       const double&, const double&, const double&,
			       const double&, const double&), 
	      double theta[], const int& index, const double& z, 
	      const double& w, const int& p, const double& r0, 
	      const double& r1, const double& c0, const double& mu0, 
	      const double& mu1, const double& sigma0, const double& sigma1, 
	      rng* stream, double& L, double& R){

  const double U = stream->runif();
  const double x0 = theta[index];
  double theta_L[2];
  double theta_R[2];
  theta_L[0] = theta_R[0] = theta[0];
  theta_L[1] = theta_R[1] = theta[1];
  L = x0 - w*U;
  theta_L[index] = L;
  R = L + w;
  theta_R[index] = R;
  int K = p;
  while (K > 0 && 
	 (z < logfun(theta_L, r0, r1, c0, mu0, mu1, sigma0, sigma1) | 
	  z < logfun(theta_R, r0, r1, c0, mu0, mu1, sigma0, sigma1))){
    double V = stream->runif();
    if (V < 0.5){
      L = L - (R - L);
      theta_L[index] = L;
    }
    else {
      R = R + (R - L);
      theta_R[index] = R;
    }
    --K;
  }  
}

// Radford Neal's (2000) Accept procedure coded for a logdensity
static const bool Accept(double (*logfun)(double[], const double&, 
					  const double&,
					  const double&, const double&, 
					  const double&,
					  const double&, const double&), 
			 double theta[], const int& index, const double x0, 
			 const double& z, const double& w, const double& r0, 
			 const double& r1, const double& c0, 
			 const double& mu0, const double& mu1, 
			 const double& sigma0, const double& sigma1, 
			 const double& L, const double& R){

  double Lhat = L;
  double Rhat = R;
  bool D = false;
  while ((Rhat - Lhat ) > 1.1 * w){
    double M = (Lhat + Rhat) / 2.0;
    if ( (x0 < M && theta[index] >= M) || (x0 >= M && theta[index] < M)){
      D = true;
    }
    if (theta[index] < M){
      Rhat = M;
    }
    else {
      Lhat = M;
    }
    int ind0;
    if (index==0){
      ind0 = 1;
    }
    else {
      ind0 = 0;
    }
    double theta_L[2];
    double theta_R[2];
    theta_L[ind0] = theta_R[ind0] = theta[ind0];
    theta_L[index] = Lhat;
    theta_R[index] = Rhat;
    if (D && z >= logfun(theta_L, r0, r1, c0, mu0, mu1, sigma0, sigma1) && 
	z >=  logfun(theta_R, r0, r1, c0, mu0, mu1, sigma0, sigma1)){
      return(false);
    }    
  }
  return(true);
}


// Radford Neal's (2000) shrinkage procedure coded for a log density
static double shrinkage(double (*logfun)(double[], const double&, 
					 const double&,
					 const double&, const double&, 
					 const double&,
					 const double&, const double&), 
			double theta[], const int& index, const double& z, 
			const double& w, const double& r0, 
			const double& r1, const double& c0, const double& mu0, 
			const double& mu1, const double& sigma0, 
			const double& sigma1, rng* 
			stream, const double& L, const double& R){

  double Lbar = L;
  double Rbar = R;
  //  int ind0;
  //if (index==0){
  //  ind0 = 1;
  // }
  //else {
  //  ind0 = 0;
  // }
  double theta_x1[2];
  theta_x1[0] = theta[0];
  theta_x1[1] = theta[1];
  const double x0 = theta[index]; 
  for (;;){
    const double U = stream->runif();
    const double x1 = Lbar + U*(Rbar - Lbar);
    theta_x1[index] = x1;
    if (z < logfun(theta_x1, r0, r1, c0, mu0, mu1, sigma0, sigma1) &&
	Accept(logfun, theta_x1, index, x0, z, w, 
	       r0, r1, c0, mu0, mu1, 
	       sigma0, sigma1, L, R)){
      return(x1);
    }
    if (x1 < x0){
      Lbar = x1;
    }
    else {
      Rbar = x1;
    }
  } // end infinite loop
}




extern "C"{
  
  
  void hierEI(double* sample, const int* samrow, const int* samcol,
	      const double* Rr0, const double* Rr1, const double* Rc0,
	      const double* Rc1, const int* Rntables, const int* Rburnin,
	      const int* Rmcmc, const int* Rthin, 
	      const double* Rmu0pm, const double* Rmu0pv,
	      const double* Rmu1pm, const double* Rmu1pv,
	      const double* Rnu0, const double* Rdelta0,
	      const double* Rnu1, const double* Rdelta1,
	      const int* Rverbose, 
	      const int *lecuyer, const int *seedarray,
	      const int *lecuyerstream){
    

    // load data
    // table notation is:
    // --------------------
    //   Y0  |     | r0
    // --------------------
    //   Y1  |     | r1
    // --------------------
    //   c0  | c1  | N
  
     // initialize rng stream
     rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);
    
    const int ntables = *Rntables;
    const int verbose = *Rverbose;

    Matrix<double> r0(ntables, 1, Rr0);
    Matrix<double> r1(ntables, 1, Rr1);
    Matrix<double> c0(ntables, 1, Rc0);
    Matrix<double> c1(ntables, 1, Rc1);
    Matrix<double> N = c0 + c1;

    
    // MCMC-related quantities
    int burnin = *Rburnin;
    int mcmc =   *Rmcmc;
    int thin =   *Rthin;
    int tot_iter = burnin + mcmc;

    
    // prior for mu0 ~ N(mu0_prior_mean, mu0_prior_var)
    double mu0_prior_mean = *Rmu0pm;
    double mu0_prior_var = *Rmu0pv;

    // prior for mu1 ~ N(mu1_prior_mean, mu1_prior_var)
    double mu1_prior_mean = *Rmu1pm;
    double mu1_prior_var = *Rmu1pv;
    
    // prior for sigma0 ~ IG(nu0/2, delta0/2)
    double nu0 = *Rnu0;
    double delta0 = *Rdelta0;

    // prior for sigma1 ~ IG(nu1/2, delta1/2)
    double nu1 = *Rnu1;
    double delta1 = *Rdelta1;
    
    // storage matrices
    Matrix<double> p0mat(mcmc/thin, ntables);
    Matrix<double> p1mat(mcmc/thin, ntables);
    Matrix<double> mu0mat(mcmc/thin, 1);
    Matrix<double> mu1mat(mcmc/thin, 1);
    Matrix<double> sig0mat(mcmc/thin, 1);
    Matrix<double> sig1mat(mcmc/thin, 1);
    int count = 0;
    
    // starting values
    Matrix<double> p0 = stream->runif(ntables,1)*0.5 + 0.25;
    Matrix<double> p1 = stream->runif(ntables,1)*0.5 + 0.25;
    Matrix<double> theta0 = log(p0/(1.0 - p0));
    Matrix<double> theta1 = log(p1/(1.0 - p1));
    double mu0 = 0.0;
    double mu1 = 0.0;
    double sigma0 = 1.0;
    double sigma1 = 1.0;
    double L = -2.0;
    double R = 2.0;
    
    // sampling constants
    const int warmup_iter = 4000;
    const int warmup_burnin = 2000;
    const double w_init = .000000001;
    const int p_init = 50;
    const Matrix<double> widthmat(warmup_iter - warmup_burnin, 2);

    // warm up sampling to choose slice sampling parameters adaptively
    for (int iter=0; iter<warmup_iter; ++iter){

      // loop over tables
      for (int i=0; i<ntables; ++i){

	// sample theta0, theta1 using slice sampling
	for (int index = 0; index<2; ++index){
	  double theta_i[2];
	  theta_i[0] = theta0[i];
	  theta_i[1] = theta1[i];
	  double funval = Lev1thetaPost(theta_i, r0[i], r1[i], c0[i], 
					mu0, mu1, sigma0, sigma1);
	  
	  double z = funval - stream->rexp(1.0);
	  doubling(&Lev1thetaPost, theta_i, index, z, w_init, p_init, r0[i], 
		   r1[i], c0[i], mu0, mu1, sigma0, sigma1, stream, L, R);

	  theta_i[index] = shrinkage(&Lev1thetaPost, theta_i, index, z, 
				     w_init, r0[i], r1[i], c0[i], mu0, mu1, 
				     sigma0, sigma1, stream, L, R);
	  if (iter >= warmup_burnin){
	    widthmat(iter-warmup_burnin, index) =  R - L;
	  }
	  theta0[i] = theta_i[0];
	  theta1[i] = theta_i[1];	  
	} // end index loop	  
	
      } // end tables loop
      
      // sample mu0 and mu1
      // mu0
      double post_var = 1.0/(1.0/mu0_prior_var + ntables*(1.0/sigma0));
      double post_mean = post_var*(sumc(theta0)[0]*(1.0/sigma0) + 
				   (1.0/mu0_prior_var)*mu0_prior_mean);
      mu0 = stream->rnorm(post_mean, sqrt(post_var));
      
      // mu1
      post_var = 1.0/(1.0/mu1_prior_var + ntables*(1.0/sigma1));
      post_mean = post_var*(sumc(theta1)[0]*(1.0/sigma1) + 
			    (1.0/mu1_prior_var)*mu1_prior_mean);
      mu1 = stream->rnorm(post_mean, sqrt(post_var));
           
      // sample sigma0 and sigma1
      // sigma0
      Matrix<double> e = theta0 - mu0;
      Matrix<double> SSE = crossprod(e);
      double nu2 = (nu0 + ntables)*0.5;
      double delta2 = (delta0 + SSE[0])*0.5;
      sigma0 = stream->rigamma(nu2,delta2);
      
      // sigma1
      e = theta1 - mu1;
      SSE = crossprod(e);
      nu2 = (nu1 + ntables)*0.5;
      delta2 = (delta1 + SSE[0])*0.5;
      sigma1 = stream->rigamma(nu2,delta2);           

      // allow user interrupts
      R_CheckUserInterrupt();        
    }
    // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


    // sampling constants
    const double w = mean(widthmat);
    int p_temp = 2;
    while ((w * pow(2.0, p_temp) ) < max(widthmat)){
      ++p_temp;
    } 
    const int p = p_temp + 1;

    // @@@@@@@@@@@@@@  The real sampling @@@@@@@@@@@@@@@    
    for (int iter=0; iter<tot_iter; ++iter){

      // loop over tables
      for (int i=0; i<ntables; ++i){

	// sample theta0, theta1 using slice sampling
	for (int index = 0; index<2; ++index){
	  double theta_i[2];
	  theta_i[0] = theta0[i];
	  theta_i[1] = theta1[i];
	  double funval = Lev1thetaPost(theta_i, r0[i], r1[i], c0[i], 
					mu0, mu1, sigma0, sigma1);
	  
	  double z = funval - stream->rexp(1.0);
	  doubling(&Lev1thetaPost, theta_i, index, z, w, p, r0[i], 
		   r1[i], c0[i], mu0, mu1, sigma0, sigma1, stream, L, R);

	  //Rprintf("L = %10.5f  R = %10.5f\n", L, R);

	  theta_i[index] = shrinkage(&Lev1thetaPost, theta_i, index, z, w, 
				     r0[i], r1[i], c0[i], mu0, mu1, 
				     sigma0, sigma1, stream, L, R);


	  theta0[i] = theta_i[0];
	  theta1[i] = theta_i[1];	  
	} // end index loop	  
	
	// if after burnin store samples
	if ((iter >= burnin) && ((iter%thin)==0)){
	  p0mat(count,i) = 1.0/(1.0 + exp(-1*theta0[i]));
	  p1mat(count,i) = 1.0/(1.0 + exp(-1*theta1[i]));
	  
	}
      } // end tables loop
      

      // sample mu0 and mu1
      // mu0
      double post_var = 1.0/(1.0/mu0_prior_var + ntables*(1.0/sigma0));
      double post_mean = post_var*(sumc(theta0)[0]*(1.0/sigma0) + 
				   (1.0/mu0_prior_var)*mu0_prior_mean);
      mu0 = stream->rnorm(post_mean, sqrt(post_var));
      
      // mu1
      post_var = 1.0/(1.0/mu1_prior_var + ntables*(1.0/sigma1));
      post_mean = post_var*(sumc(theta1)[0]*(1.0/sigma1) + 
			    (1.0/mu1_prior_var)*mu1_prior_mean);
      mu1 = stream->rnorm(post_mean, sqrt(post_var));
      
     
      // sample sigma0 and sigma1
      // sigma0
      Matrix<double> e = theta0 - mu0;
      Matrix<double> SSE = crossprod(e);
      double nu2 = (nu0 + ntables)*0.5;
      double delta2 = (delta0 + SSE[0])*0.5;
      sigma0 = stream->rigamma(nu2,delta2);
      
      // sigma1
      e = theta1 - mu1;
      SSE = crossprod(e);
      nu2 = (nu1 + ntables)*0.5;
      delta2 = (delta1 + SSE[0])*0.5;
      sigma1 = stream->rigamma(nu2,delta2);
      
      
      // if after burnin store samples
      if ((iter >= burnin) && ((iter%thin)==0)){
	mu0mat(count,0) = mu0;
	mu1mat(count,0) =  mu1;
	sig0mat(count,0) = sigma0;
	sig1mat(count,0) = sigma1;
	++ count;
      }
      
      // print output to screen
      if (verbose>0 && (iter%verbose)==0){
	Rprintf("\nMCMChierEI iteration %i of %i \n", (iter+1), 
		tot_iter);
      }
      
      // allow user interrupts
      R_CheckUserInterrupt();        
    }
    
    delete stream; // clean up random number stream   
    
    // return sample
    Matrix<double> storeagem = cbind(p0mat, p1mat);
    storeagem = cbind(storeagem, mu0mat);
    storeagem = cbind(storeagem, mu1mat);
    storeagem = cbind(storeagem, sig0mat);
    storeagem = cbind(storeagem, sig1mat);
    int mat_size = samrow[0] * samcol[0];
    for (int i=0; i<mat_size; ++i)
      sample[i] = storeagem[i];

  }

} // extern "C"

