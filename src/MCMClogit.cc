// sample from the posterior distribution of a logistic regression
// logistice regression model in R using linked C++ code in Scythe
//
// KQ 1/23/03

#include <iostream>
#include "Scythe_Matrix.h"
#include "Scythe_Simulate.h"
#include "Scythe_Stat.h"
#include "Scythe_LA.h"
#include "Scythe_IDE.h"
#include "Scythe_Math.h"

using namespace SCYTHE;
using namespace std;

double logit_logpost(Matrix<double>& Y, const Matrix<double>& X, 
		     const Matrix<double>& beta,
		     const Matrix<double>& beta_prior_mean, 
		     const Matrix<double>& beta_prior_var){

  Matrix<double> eta = X * beta;
  Matrix<double> p = exp(eta)/(1+exp(eta));
  double loglike = 0;
  for (int i=0; i<Y.rows(); ++i)
    loglike += Y[i]*::log(p[i]) + (1-Y[i])*::log(1-p[i]);
  double logprior = lndmvn(beta, beta_prior_mean, beta_prior_var);
  double logpost = loglike + logprior;
  return logpost;
}

extern "C"{

  using namespace SCYTHE;

  // simulate from posterior density and return a Gibbs by parameters matrix 
  // of the posterior density sample
  void
  logitpost (double* samdata, const int* samrow, const int* samcol, 
	     const double* Xdata, const int* Xrow, const int* Xcol,
	     const double* Ydata, const int* Yrow, const int* Ycol,
	     const int* burnin, const int* mcmc, const int* thin,
	     const int* seed, const int* verbose,
	     const int* bstartdata, const int* bstartrow, const int* bstartcol, 
	     const double* b0data, const int* b0row, const int* b0col,
	     const double* B0data, const int* B0row, const int* B0col,
	     const double* mdata, const int* mrow, const int* mcol,
	     const double* Vdata, const int* Vrow, const int* Vcol,
	     const double* tune, int* accepts){

    // initialize seed (mersenne twister / use default seed unless specified)
    if(seed==0) set_mersenne_seed(5489UL);
    else set_mersenne_seed(seed[0]);

    // form matrices
    Matrix<double> X(Xcol[0], Xrow[0], Xdata);
    X = t(X);
    Matrix<double> Y(Ycol[0], Yrow[0], Ydata);
    Y = t(Y);
    Matrix<double> b0(b0col[0], b0row[0], b0data);
    b0 = t(b0);
    Matrix<double> B0(B0col[0], B0row[0], B0data);
    B0 = t(B0);
    Matrix<double> m(mcol[0], mrow[0], mdata);
    m = t(m);
    Matrix<double> V(Vcol[0], Vrow[0], Vdata);
    V = t(V);
    Matrix<double> C = cholesky(V) * tune[0];

    Matrix<double> beta_prior_var = invpd(B0);

    // proposal parameters
    Matrix<double> propV = invpd(B0 + invpd(V));
    Matrix<double> propm = propV*(B0*b0 + invpd(V)*m);
    Matrix<double> propC = cholesky(propV) * tune[0];

    // define constants
    int k = X.cols();
    int tot_iter = burnin[0] + mcmc[0];

    // holding matrix
    Matrix<double> betam(k, mcmc[0]);

    // starting values
    Matrix<double> beta(bstartcol[0], bstartrow[0], bstartdata[0]);
    beta = t(beta);
  

    double logpost_cur = logit_logpost(Y,X,beta,b0,beta_prior_var);

    // MCMC loop
    int count = 0;
    for (int iter = 0; iter < tot_iter; ++iter){

      // sample beta
      Matrix<double> beta_can = beta + propC * rnorm(k,1);

      double logpost_can = logit_logpost(Y,X,beta_can, b0, beta_prior_var);
      double ratio = ::exp(logpost_can - logpost_cur); 

      if (runif() < ratio){
	beta = beta_can;
	logpost_cur = logpost_can;
	++accepts[0];
      }


      // store values in matrices
      if (iter >= burnin[0] && ((iter%thin[0])==0)){ 
	for (int j = 0; j < k; j++)
	  betam (j, count) = beta[j];
	++count;
      }


      if (verbose[0] == 1 && iter % 500 == 0) {
	cout << "MCMClogit iteration = " << iter + 1 << endl;
	cout << " beta = " << endl << beta.toString();
	cout << " acceptance rate = " << static_cast<double>(accepts[0])/
	    static_cast<double>(iter+1) << endl << endl;
      }
    

    }// end MCMC loop
  
    // return output
    Matrix<double> holder = t(betam);
    int number = samrow[0] * samcol[0];
    for (int i=0; i<number; ++i)
      samdata[i] = holder[i];
  }
  
}

