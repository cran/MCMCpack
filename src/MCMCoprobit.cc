// MCMCprobit.cc is a program that simulates draws from the posterior
// density of a probit regression model.  It is written
// modularly, with the goal of creating an R object that will reference
// this code.
// 
// Kevin Quinn
// Dept. of Political Science and CSSS
// University of Washington
// quinn@stat.washington.edu
// 
// Andrew D. Martin
// Dept. of Political Science
// 
// 1/24/2003 
// 

#include <iostream>
#include "Scythe_Matrix.h"
#include "Scythe_Simulate.h"
#include "Scythe_Stat.h"
#include "Scythe_LA.h"
#include "Scythe_IDE.h"


using namespace SCYTHE;
using namespace std;


extern "C"{

  using namespace SCYTHE;


  // simulate from posterior density and return a Gibbs by parameters matrix 
  // of the posterior density sample
  void
  oprobitpost (double* sample, const int* samrow, const int* samcol,
	       const int* Y, const int* Yrow, const int* Ycol,
	       const double* X, const int* Xrow, const int* Xcol,
	       const int* burnin, const int* mcmc,  const int* thin,
	       const double* tune, const int* seed, const int* verbose, 
	       const double* bstart, const int* bstartrow, 
	       const int* bstartcol,
	       const double* gammastart, const int* gammastartrow, 
	       const int* gammastartcol,
	       const double* b0, const int* b0row, const int* b0col,
	       const double* B0, const int* B0row, const int* B0col,
	       int* accepts){

    // put together matrices
    Matrix<double> Msample(samcol[0], samrow[0], sample);
    Msample = t(Msample);
    Matrix<double> MX(Xcol[0], Xrow[0], X);
    MX = t(MX);
    Matrix<int> MY(Ycol[0], Yrow[0], Y);
    MY = t(MY);
    Matrix<double> Mbetastart(bstartcol[0],bstartrow[0], bstart);
    Mbetastart = t(Mbetastart);
    Matrix<double> Mgammastart(gammastartcol[0], gammastartrow[0], gammastart);
    Mgammastart = t(Mgammastart);
    Matrix<double> Mb0(b0col[0], b0row[0], b0);
    Mb0 = t(Mb0);
    Matrix<double> MB0(B0col[0], B0row[0], B0);
    MB0 = t(MB0);

    // define constants and from cross-product matrices
    int k = MX.cols();
    int ncat = Mgammastart.rows() - 1;
    int N = MX.rows();
    int tot_iter = burnin[0] + mcmc[0];
    Matrix<double> XpX = crossprod (MX);


    // holding matrices
    Matrix<double> betam(mcmc[0]/thin[0], k);
    Matrix<double> gammam(mcmc[0]/thin[0], ncat+1);
    
    // initialize seed (mersenne twister / use default seed unless specified)
    if(seed==0) set_mersenne_seed(5489UL);
    else set_mersenne_seed(seed[0]);
     
    // starting values
    Matrix<double> beta = Mbetastart;
    Matrix<double> gamma = Mgammastart;
    Matrix<double> gamma_p = gamma;
    Matrix<double> Z(N,1);
    Matrix<double> Xbeta = MX * beta;
    

    // Gibbs loop
    int count = 0;
    for (int iter = 0; iter < tot_iter; ++iter)
      {



	// [gamma | Z, beta]
	for (int i=2; i<(ncat); ++i){
	  if (i==(ncat-1)){
	    gamma_p[i] = rtbnorm_combo(gamma[i], ::pow(tune[0], 2.0), 
				       gamma_p[i-1]);
	  }
	  else {
	    gamma_p[i] = rtnorm(gamma[i], ::pow(tune[0], 2.0), gamma_p[i-1], 
				gamma[i+1]);
	  }
	}
	double loglikerat = 0.0;
	double loggendenrat = 0.0;
	
	// loop over observations and construct the acceptance ratio
	for (int i=0; i<N; ++i){
	  if (MY[i] == ncat){
	    loglikerat = loglikerat 
	      + log(1.0  - 
		    pnorm1(gamma_p[MY[i]-1] - Xbeta[i]) ) 
	      - log(1.0 - 
		    pnorm1(gamma[MY[i]-1] - Xbeta[i]) );
	  }
	  else if (MY[i] == 1){
	    loglikerat = loglikerat 
	      + log(pnorm1(gamma_p[MY[i]] - Xbeta[i])  ) 
	      - log(pnorm1(gamma[MY[i]] - Xbeta[i]) );
	  }
	  else{
	    loglikerat = loglikerat 
	      + log(pnorm1(gamma_p[MY[i]] - Xbeta[i]) - 
		    pnorm1(gamma_p[MY[i]-1] - Xbeta[i]) ) 
	      - log(pnorm1(gamma[MY[i]] - Xbeta[i]) - 
		    pnorm1(gamma[MY[i]-1] - Xbeta[i]) );
	  }
	}
	for (int j=2; j<(ncat-1); ++j){	   
	  loggendenrat = loggendenrat 
	    + log(pnorm(gamma[j+1], gamma[j], tune[0]) - 
		  pnorm(gamma[j-1], gamma[j], tune[0]) )  
	    - log(pnorm(gamma_p[j+1], gamma_p[j], tune[0]) - 
		  pnorm(gamma_p[j-1], gamma_p[j], tune[0]) );
	}
	double logacceptrat = loglikerat + loggendenrat;
	if (runif() <= exp(logacceptrat)){
	  gamma = gamma_p;
	  ++accepts[0];
	}
	

	// [Z| gamma, beta, y] 
	Matrix<double> Z_mean = MX * beta;
	for (int i=0; i<N; ++i){
	  Z[i] = rtnorm(Z_mean[i], 1.0, gamma[MY[i]-1], gamma[MY[i]]);
	}
	
      
	// [beta|Z, gamma]
	Matrix<double> XpZ = t(MX) * Z;
	Matrix<double> sig_beta = invpd (MB0 + XpX );
	Matrix<double> C = cholesky (sig_beta);
	Matrix<double> betahat = sig_beta * (XpZ + MB0 * Mb0);
	beta = betahat + C * rnorm (k, 1);
	Xbeta = MX * beta;

	// store values in matrices
	if (iter >= burnin[0] && ((iter%thin[0])==0)){ 
	  for (int j = 0; j < k; j++)
	    betam (count,j) = beta[j];
	  for (int j = 0; j<(ncat+1); ++j)
	    gammam (count, j) = gamma[j];;
	  ++count;
	}
	  
	if (verbose[0] == 1 && iter % 500 == 0) {
	  cout << "MCMCoprobit iteration = " << iter + 1 << endl;
	  cout << " beta = " << endl << beta.toString();
	  cout << " acceptance rate = " << static_cast<double>(accepts[0])/
	    static_cast<double>(iter+1) << endl << endl;
	}
      }

    // return output
    Matrix<double> datam = cbind(betam, gammam);
    int loop = samrow[0] * samcol[0];
    for (int i=0; i<loop; ++i)
      sample[i] = datam[i];

  }
  
}
