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
// 11/12/2002 -- ported to Scythe0.3
// 

#include <iostream>
#include "Scythe_Matrix.h"
#include "Scythe_Simulate.h"
#include "Scythe_Stat.h"
#include "Scythe_LA.h"
#include "Scythe_IDE.h"

using namespace SCYTHE;


extern "C"{

  using namespace SCYTHE;
	using namespace std;


  // simulate from posterior density and return a Gibbs by parameters matrix 
  // of the posterior density sample
  void
  probitpost (double* sample, const int* samrow, const int* samcol,
	      const double* X, const int* Xrow, const int* Xcol,
	      const double* Y, const int* Yrow, const int* Ycol,
	      const int* burnin, const int* gibbs,  const int* thin,
	      const int* seed, const int* verbose, const double* bstart,
	      const int* bstartrow, const int* bstartcol,
	      const double* b0, const int* b0row, const int* b0col,
	      const double* B0, const int* B0row, const int* B0col){
    // put together matrices
    Matrix<double> Msample(samcol[0], samrow[0], sample);
    Msample = t(Msample);
    Matrix<double> MX(Xcol[0], Xrow[0], X);
    MX = t(MX);
    Matrix<double> MY(Ycol[0], Yrow[0], Y);
    MY = t(MY);
    Matrix<double> Mbetastart(bstartcol[0],bstartrow[0], bstart);
    Mbetastart = t(Mbetastart);
    Matrix<double> Mb0(b0col[0], b0row[0], b0);
    Mb0 = t(Mb0);
    Matrix<double> MB0(B0col[0], B0row[0], B0);
    MB0 = t(MB0);

    // define constants and from cross-product matrices
    int k = MX.cols();
    int N = MX.rows();
    int tot_iter = burnin[0] + gibbs[0];
    Matrix<double> XpX = crossprod (MX);

    // holding matrices
    Matrix<double> betam(gibbs[0]/thin[0], k);

    // initialize seed (mersenne twister / use default seed unless specified)
    if(seed==0) set_mersenne_seed(5489UL);
    else set_mersenne_seed(seed[0]);
     
    // starting values
    Matrix<double> beta = Mbetastart;
    Matrix<double> Z(N,1);

    // Gibbs loop
    int count = 0;
    for (int iter = 0; iter < tot_iter; ++iter)
      {

	// [Z| beta, y]
	Matrix<double> Z_mean = MX * beta;
	for (int i=0; i<N; ++i){
	  if (MY[i] == 1.0)
	    Z[i] = rtbnorm_combo(Z_mean[i], 1.0, 0);
	  if (MY[i] == 0.0)
	    Z[i] = rtanorm_combo(Z_mean[i], 1.0, 0);
	}
      
	// [beta|Z]
	Matrix<double> XpZ = t(MX) * Z;
	Matrix<double> sig_beta = invpd (MB0 + XpX );
	Matrix<double> C = cholesky (sig_beta);
	Matrix<double> betahat = sig_beta * (XpZ + MB0 * Mb0);
	beta = betahat + C * rnorm (k, 1);
	

	// store values in matrices
	if (iter >= burnin[0] && ((iter%thin[0])==0)){ 
	  for (int j = 0; j < k; j++)
	    betam (count,j) = beta[j];
	  ++count;
	}
	  
	if (verbose[0] == 1 && iter % 500 == 0) {
	  cout << "MCMCprobit iteration = " << iter + 1 << endl;
	  cout << " beta = " << endl << beta.toString() << endl << endl;
	}
      }

    // return output
    int loop = samrow[0] * samcol[0];
    for (int i=0; i<loop; ++i)
      sample[i] = betam[i];

  }
  
}
