// fits a model derived from Wakefield's baseline model for 
// ecological inference in which logit(p_i) follows a random walk in time
// a priori. The model is fit using Wakefield's normal approximation 
// to the binomial convolution likelihood and the Metropolis-Hastings 
// algorithm to sample from the posterior
//
//  evolution variances are estimated 
//
// Kevin Quinn
// Dept. of Political Science and CSSS
// Univ. Washington
// quinn@stat.washington.edu
//
// 3/9/2002
// 10/25/2002 -- ported to Scythe0.3 and written for an R interface
//



#include <iostream> 
#include "Scythe_Matrix.h"
#include "Scythe_Simulate.h"
#include "Scythe_Stat.h"
#include "Scythe_Math.h"
#include "Scythe_LA.h"
#include "Scythe_IDE.h"



extern "C"{
 
  using namespace SCYTHE;
	using namespace std;
  
  void dynamicEI(double* sample, const int* samrow, const int* samcol,
		  const double* Rr0, const double* Rr1, const double* Rc0,
		  const double* Rc1, const int* Rntables, const int* Rburnin,
		  const int* Rmcmc, const int* Rthin, 
		  const double* RW, const double* Rnu0,
		  const double* Rdelta0, const double* Rnu1, 
		  const double* Rdelta1, const int* Rverbose, 
		  const double* Rtune, const int* Rseed, int* accepts){


    // load data
    // table notation is:
    // --------------------
    //   Y0  |     | r0
    // --------------------
    //   Y1  |     | r1
    // --------------------
    //   c0  | c1  | N

   
    // initialize seed (mersenne twister / use default seed unless specified)
    if(Rseed[0]==0) set_mersenne_seed(5489UL);
    else set_mersenne_seed(Rseed[0]);

    
    int ntables = *Rntables;
    int verbose = *Rverbose;

    Matrix<double> r0(ntables, 1, Rr0);
    Matrix<double> r1(ntables, 1, Rr1);
    Matrix<double> c0(ntables, 1, Rc0);
    Matrix<double> c1(ntables, 1, Rc1);
    Matrix<double> N = c0 + c1;
    
    Matrix<double> W(ntables, ntables, RW);
    
    // MCMC-related quantities
    int burnin = *Rburnin;
    int mcmc =   *Rmcmc;
    int thin =   *Rthin;
    int tot_iter = burnin + mcmc;
    double tune = *Rtune;
    
    
    // sum of weights across each row
    Matrix<double> W_sum = sumc(t(W));  
    
    // precision matrix (not the weight matrix) for theta0 and theta1
    Matrix<double> P = -1*W;
    for (int i=0; i<ntables; ++i)
      P(i,i) = W_sum[i];
  

    // sigma_theta0 ~ IG(nu0/2, delta0/2)
    double nu0 = *Rnu0;
    double delta0 = *Rdelta0;

    // sigma_theta1 ~ IG(nu1/2, delta1/2)
    double nu1 = *Rnu1;
    double delta1 = *Rdelta1;


    // storage matrices
    Matrix<double> p0mat(mcmc/thin, ntables);
    Matrix<double> p1mat(mcmc/thin, ntables);
    Matrix<double> sig0mat(mcmc/thin, 1);
    Matrix<double> sig1mat(mcmc/thin, 1);

    int count = 0;

    // starting values
    Matrix<double> p0 = ones<double>(ntables,1)*0.5;
    Matrix<double> p1 = ones<double>(ntables,1)*0.5;
    Matrix<double> y0(ntables,1);
    Matrix<double> y1(ntables,1);
    for (int j=0; j<ntables; ++j){
      // sample y0|c1,r0,r1,p0,p1
      double psi = ( p0[j]*(1.0-p1[j]) ) / ( p1[j]*(1.0-p0[j]));
      y0[j] = rnchypgeom(c0[j], r0[j], r1[j], psi);
      y1[j] = c0[j] - y0[j];
    
      // sample (p0,p1)|y0,y1,r0,r1,c0,c1
      p0[j] = rbeta(1+y0[j], 1+(r0[j]-y0[j]));
      p1[j] = rbeta(1+y1[j], 1+(r1[j]-y1[j]));
    }
    Matrix<double> theta0 = log(p0/(1.0 - p0));
    Matrix<double> theta1 = log(p1/(1.0 - p1));
    Matrix<double> logjumpdens_cur = ones<double>(ntables,1)*1e200;
    // evolution variance for theta0
    double sigma_theta0 = ::pow(0.25, 2);
    // evolution variance for theta1
    double sigma_theta1 = ::pow(0.25, 2);

    // tomography line quantities
    // p0 on X axis and p1 on Y axis
    Matrix<double> r0frac = r0/(r0+r1);
    Matrix<double> r1frac = r1/(r0+r1);
    Matrix<double> c0frac = c0/(c0+c1);
    Matrix<double> c1frac = c1/(c0+c1);
    Matrix<double> intercept = c0/r1;   // intercept of tomography line 
    Matrix<double> slope = -1*(r0/r1);  // slope of tomography line
    Matrix<double> orthoSD = tune/sqrt(N); // sd for MH sampline
    Matrix<double> orthoVar = pow(orthoSD, 2);
    Matrix<double> p0min(ntables,1);
    Matrix<double> p0max(ntables,1);

    // calculate min and max possible values of p0
    for (int i=0; i<ntables; ++i){
      p0min[i] = SCYTHE::max((1.0 - intercept[i])/slope[i], 0.0);
      p0max[i] = SCYTHE::min(-1*intercept[i]/slope[i], 1.0);
    }


    // adjust p0min and p0max so all of [0,1]^2 is covered by the 
    // proposal density in the MH step
    for (int i=0; i<ntables; ++i){
      double run = ::sqrt(4.0 * ::pow(orthoSD[i], 2) * 
			  (1.0 + 1.0/::pow(slope[i], 2))) /
	(2.0*(1.0 + 1.0 / ::pow(slope[i], 2)));
      p0min[i] = p0min[i] - run;
      p0max[i] = p0max[i] + run;
    }


    for (int iter=0; iter<tot_iter; ++iter){
      for (int i=0; i<ntables; ++i){
      
	// sample (p0,p1)|r0,r1,c0,c1
      
	// sample candidate values of p0 and p1
	double u = runif()*(p0max[i]-p0min[i]) + p0min[i]; 
	double length = rnorm(0.0, orthoSD[i]);
	double s = sgn(length);
	length = fabs(length);
	double run = s * ::sqrt(4.0 * ::pow(length, 2) * 
				(1.0 + 1.0/::pow(slope[i], 2))) /
	  (2.0*(1.0 + 1.0 / ::pow(slope[i], 2)));
	double rise = -1.0 * run/slope[i];
	// the candidate values
	double p0_can =  u + run;
	double p1_can = intercept[i] + slope[i]*(u) + rise;
	double theta0_can = ::log(p0_can/(1.0 - p0_can));
	double theta1_can = ::log(p1_can/(1.0 - p1_can));
      
	// log density ordinates
	double logjumpdens_can = lndnorm(length, 0.0, orthoSD[i]);
	double logprior_can, loglike_can, logpost_can, prior_mean0, 
	  prior_mean1, prior_var0, prior_var1;
	if ((p0_can < 1.0) && (p0_can > 0.0) && (p1_can < 1.0) && 
	    (p1_can > 0.0)){	
	  prior_mean0 = ((W(i,_) * theta0) / W_sum[i])[0];
	  prior_mean1 = ((W(i,_) * theta1) / W_sum[i])[0];
	  prior_var0 = sigma_theta0/W_sum[i];
	  prior_var1 = sigma_theta1/W_sum[i];
	  logprior_can = lndnorm(theta0_can, prior_mean0, sqrt(prior_var0)) +
	    lndnorm(theta1_can, prior_mean1, sqrt(prior_var1));

	  loglike_can = lndnorm(c0[i], r0[i]*p0_can + r1[i]*p1_can,
				 sqrt(r0[i]*p0_can*(1.0-p0_can) + 
				      r1[i]*p1_can*(1.0-p1_can)));
	  logpost_can = loglike_can + logprior_can;
	}
	else{
	  logpost_can = ::log(0);
	}
      
	prior_mean0 = ((W(i,_) * theta0) / W_sum[i])[0];
	prior_mean1 = ((W(i,_) * theta1) / W_sum[i])[0];
	prior_var0 = sigma_theta0/W_sum[i];
	prior_var1 = sigma_theta1/W_sum[i];
	double logprior_cur = lndnorm(theta0[i], prior_mean0, sqrt(prior_var0)) +
	  lndnorm(theta1[i], prior_mean1, sqrt(prior_var1));
      

	double loglike_cur = lndnorm(c0[i], r0[i]*p0[i] + r1[i]*p1[i],
				     sqrt(r0[i]*p0[i]*(1.0-p0[i]) + 
					  r1[i]*p1[i]*(1.0-p1[i])));
      
	double logpost_cur = loglike_cur + logprior_cur;
      

	double alpha = ::exp(logpost_can - logpost_cur + logjumpdens_cur[i] - 
			     logjumpdens_can);
      
      
	if (runif() < alpha){
	  p0[i] = p0_can;
	  p1[i] = p1_can;
	  theta0[i] = theta0_can;
	  theta1[i] = theta1_can;
	  logjumpdens_cur[i] = logjumpdens_can;
	  ++accepts[0];
	}


	// sample sigma_theta0 and sigma_theta1
	Matrix<double> SSE = t(theta0-meanc(theta0)) * P * 
	  (theta0 - meanc(theta0));
	double nu2 = (nu0 + ntables)*0.5;
	double delta2 = (delta0 + SSE[0])*0.5;
	sigma_theta0 = rigamma(nu2, delta2);
      
	SSE = t(theta1-meanc(theta1)) * P * (theta1 - meanc(theta1));
	nu2 = (nu1 + ntables)*0.5;
	delta2 = (delta1 + SSE[0])*0.5;
	sigma_theta1 = rigamma(nu2, delta2);
      

	// if after burnin store samples
	if ((iter >= burnin) && ((iter%thin)==0)){
	  p0mat(count,i) = p0[i];
	  p1mat(count,i) = p1[i];
	}
      }
    

      if ((iter >= burnin) && ((iter%thin)==0)){
	sig0mat(count,0) = sigma_theta0;
	sig1mat(count,0) = sigma_theta1;
	++count;
      }  

    
      // print output to screen
      if (verbose==1 && (iter%10000)==0){
	cout << "MCMCdynamicEI iteration = " << iter <<  endl;
	cout << " sigma_theta0 = " << sigma_theta0 << endl;
	cout << " sigma_theta1 = " << sigma_theta1 << endl;
	cout << " MH acceptance rate = " << 
	  static_cast<double>(accepts[0]) / static_cast<double>(iter) 
	  / static_cast<double>(ntables) << endl << endl;
      }
    
    }


    // return sample
    Matrix<double> storeagem = cbind(p0mat, p1mat);
    storeagem = cbind(storeagem, sig0mat);
    storeagem = cbind(storeagem, sig1mat);
    int mat_size = samrow[0] * samcol[0];
    for (int i=0; i<mat_size; ++i)
      sample[i] = storeagem[i];
  

  }

} // extern "C"

