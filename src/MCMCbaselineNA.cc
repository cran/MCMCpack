// fits Wakefield's baseline model for ecological inference using
// Wakefield's normal approximation to the binomial convolution likelihood
// and Metropolis-Hastings to sample from the posterior
//
// KQ 3/2/2002
// KQ 10/22/2002 [ported to Scythe0.3 and written for an R interface]
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
  
  void baselineNA(double* sample, const int* samrow, const int* samcol,
		  const double* Rr0, const double* Rr1, const double* Rc0,
		  const double* Rc1, const int* Rntables, const int* Rburnin,
		  const int* Rmcmc, const int* Rthin, const double* Ralpha0,
		  const double* Rbeta0, const double* Ralpha1, 
		  const double* Rbeta1, const int* Rverbose, 
		  const double* Rtune, const int* Rseed, int* accepts){


    // load data
    // table notation is:
    // --------------------
    //   Y0    |     | r0
    // --------------------
    //   Y1    |     | r1
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

    
    // MCMC-related quantities
    int burnin = *Rburnin;
    int mcmc =   *Rmcmc;
    int thin =   *Rthin;
    int tot_iter = burnin + mcmc;
    double tune = *Rtune;

    // prior for p0
    // p0 ~ beta(alpha0, beta0)
    double alpha0 = *Ralpha0;
    double beta0  = *Rbeta0;
    
    // prior for p1
    // p1 ~ beta(alpha1, beta1)
    double alpha1 = *Ralpha1;
    double beta1  = *Rbeta1;

    // storage matrices
    Matrix<double> p0mat(mcmc/thin, ntables);
    Matrix<double> p1mat(mcmc/thin, ntables);
    int count = 0;
    
    // starting values
    Matrix<double> p0 = ones<double>(ntables,1)*0.5;
    Matrix<double> p1 = ones<double>(ntables,1)*0.5;
    Matrix<double> y0(ntables,1);
    Matrix<double> y1(ntables,1);
    Matrix<double> logjumpdens_cur = ones<double>(ntables,1)*1e20;
    
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
	double p0_can = u + run;
	double p1_can = intercept[i] + slope[i]*(u) + rise;
      
	// log density ordinates
	double logjumpdens_can = lndnorm(length, 0.0, sqrt(orthoVar[i]));
	double logprior_can, loglike_can, logpost_can;
	if ((p0_can < 1.0) && (p0_can > 0.0) && (p1_can < 1.0) && 
	    (p1_can > 0.0)){
	  logprior_can = lndbeta1(p0_can, alpha0, beta0) + 
	    lndbeta1(p1_can, alpha1, beta1);
	  loglike_can = lndnorm(c0[i], r0[i]*p0_can + r1[i]*p1_can,
				sqrt(r0[i]*p0_can*(1.0-p0_can) + 
				     r1[i]*p1_can*(1.0-p1_can) ) );
	  logpost_can = loglike_can + logprior_can;
	}
	else{
	  logpost_can = ::log(0);
	}

	double logprior_cur = lndbeta1(p0[i], alpha0, beta0) + 
	  lndbeta1(p1[i], alpha1, beta1);
      
	double loglike_cur = lndnorm(c0[i], r0[i]*p0[i] + r1[i]*p1[i],
				     sqrt( r0[i]*p0[i]*(1.0-p0[i]) + 
					   r1[i]*p1[i]*(1.0-p1[i]) ) );
      
	double logpost_cur = loglike_cur + logprior_cur;
      

	double alpha = ::exp(logpost_can - logpost_cur + logjumpdens_cur[i] - 
			     logjumpdens_can);
	
	if (runif() < alpha){
	  p0[i] = p0_can;
	  p1[i] = p1_can;
	  logjumpdens_cur[i] = logjumpdens_can;
	  ++accepts[0];
	}
      

	// if after burnin store samples
	if ((iter >= burnin) && ((iter%thin)==0)){
	  p0mat(count,i) = p0[i];
	  p1mat(count,i) = p1[i];
	}
      }

      if ((iter >= burnin) && ((iter%thin)==0)) ++count;

      // print output to screen
      if (verbose==1 && (iter%25000)==0){
	cout << "MCMCbaselineEI iteration = " << iter <<  endl;
	cout << " Metropolis-Hastings acceptance rate = " << 
	  static_cast<double>(accepts[0]) / static_cast<double>(iter) / 
	  static_cast<double>(ntables) << endl << endl;
      }
    
    }



    // return sample
    Matrix<double> zeros(mcmc/thin, ntables);
    Matrix<double> storeagem = cbind(p0mat, p1mat);
    storeagem = cbind(storeagem, zeros);
    storeagem = cbind(storeagem, zeros);
    int mat_size = samrow[0] * samcol[0];
    for (int i=0; i<mat_size; ++i)
      sample[i] = storeagem[i];

  }

} //extern "C"

