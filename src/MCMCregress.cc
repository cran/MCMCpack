// regress.cc is a program that simulates draws from the posterior
// density of a regression model with Gaussian errors.  It is written
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
// 10/10/2002 -- ported to Scythe0.3
// 

#include <iostream>
#include "Scythe_Matrix.h"
#include "Scythe_Simulate.h"
#include "Scythe_Stat.h"
#include "Scythe_LA.h"
#include "Scythe_IDE.h"


using namespace SCYTHE;
using namespace std;

// linear regression with Gaussian errors beta update (multivariate Normal 
// prior)
Matrix <double> 
regress_beta_update (const Matrix <double> &XpX, const Matrix <double> &XpY,
		     const Matrix <double> &b0, const Matrix <double> &B0,
		     const double sigma2)
{

  // this function gets the cross-product matrix X'X and the matrix X'Y
  // to minimize the amount of computation within the function

  int k = XpX.cols ();

  Matrix <double> sig_beta = invpd (B0 + XpX * (1.0 / sigma2));
  Matrix <double> C = cholesky (sig_beta);
  Matrix <double> betahat = sig_beta * (XpY * (1.0 / sigma2) + B0 * b0);
  Matrix <double> beta_update = betahat + C * rnorm (k, 1);

  return beta_update;
}

// linear regression with Gaussian errors sigma2 update (inverse-Gamma
// prior)
double
regress_sigma2_update (const Matrix <double> &X, const Matrix <double> &Y,
		       const Matrix <double> &beta, const double nu0,
		       const double delta0)
{

  Matrix <double> e = Y - X * beta;
  Matrix <double> SSE = crossprod (e);
  double nu2 = (nu0 + X.rows ()) * 0.5;
  double delta2 = (delta0 + SSE[0]) * 0.5;
  double sigma2_update = rigamma (nu2, delta2);

  return sigma2_update;
}



extern "C"
{

  using namespace SCYTHE;


  // simulate from posterior density and return a Gibbs by parameters matrix 
  // of the posterior density sample
  void
  regpost (double *sample, const int *samrow, const int *samcol,
	   const double *X, const int *Xrow, const int *Xcol,
	   const double *Y, const int *Yrow, const int *Ycol,
	   const int *burnin, const int *gibbs, const int *thin,
	   const int *seed, const int *verbose, const double *bstart,
	   const int *bstartrow, const int *bstartcol,
	   const double *sigma2start,
	   const double *b0, const int *b0row, const int *b0col,
	   const double *B0, const int *B0row, const int *B0col,
	   const double *nu, const double *delta)
  {
    // put together matrices
    Matrix <double> Msample (samcol[0], samrow[0], sample);
    Msample = t (Msample);
    Matrix <double> MX (Xcol[0], Xrow[0], X);
    MX = t (MX);
    Matrix <double> MY (Ycol[0], Yrow[0], Y);
    MY = t (MY);
    Matrix <double> Mbetastart (bstartcol[0], bstartrow[0], bstart);
    Mbetastart = t (Mbetastart);
    Matrix <double> Mb0 (b0col[0], b0row[0], b0);
    Mb0 = t (Mb0);
    Matrix <double> MB0 (B0col[0], B0row[0], B0);
    MB0 = t (MB0);

    // define constants and from cross-product matrices
    int k = MX.cols ();
    int tot_iter = burnin[0] + gibbs[0];
    Matrix <double> XpX = crossprod (MX);
    Matrix <double> XpY = t (MX) * MY;

    // holding matrices
    Matrix <double> betam (k, gibbs[0] / thin[0]);
    Matrix <double> sigmam (1, gibbs[0] / thin[0]);

    // initialize seed (mersenne twister / use default seed unless specified)
    if (seed == 0)
      set_mersenne_seed (5489UL);
    else
      set_mersenne_seed (seed[0]);

    // starting values
    Matrix <double> beta = Mbetastart;
    double sigma2 = sigma2start[0];

    // Gibbs loop
    int count = 0;
    for (int iter = 0; iter < tot_iter; ++iter)
      {

	beta = regress_beta_update (XpX, XpY, Mb0, MB0, sigma2);
	sigma2 = regress_sigma2_update (MX, MY, beta, nu[0], delta[0]);

	// store values in matrices
	if (iter >= burnin[0] && (iter % thin[0] == 0))
	  {
	    sigmam (0, count) = sigma2;
	    for (int j = 0; j < k; j++)
	      {
		betam (j, count) = beta[j];
	      }
	    ++count;
	  }

	if (verbose[0] == 1 && iter % 500 == 0)
	  {
	    cout << "MCMCregress iteration = " << iter + 1 << endl;
	    cout << " beta = " << endl << beta.toString ();
	    cout << " sigma2 = " << sigma2 << endl << endl;
	  }
      }

    // return output
    Matrix <double> storeagem = cbind (t (betam), t (sigmam));
    int loop = samrow[0] * samcol[0];
    for (int i = 0; i < loop; ++i)
      sample[i] = storeagem[i];

  }

}
