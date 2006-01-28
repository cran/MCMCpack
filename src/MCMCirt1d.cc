// MCMCirt1d.cc is C++ code to estimate a one-dimensional item response
// theory model. 
//
// ADM and KQ 1/15/2003
// ADM 7/28/2004 [updated to new Scythe version]
// completely rewritten and optimized for the 1-d case 8/2/2004 KQ
// storage changed to save memory KQ 1/27/2006

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

extern "C"{

// function called by R to fit model
  void
  MCMCirt1d(double* sampledata, 
	    const int* samplerow, 
	    const int* samplecol,
	    const int* Xdata, 
	    const int* Xrow, 
	    const int* Xcol,
	    const int* burnin, 
	    const int* mcmc,  
	    const int* thin,
	    const int *lecuyer, 
	    const int *seedarray,
	    const int *lecuyerstream, 
	    const int* verbose, 
	    const double* thetastartdata,
	    const int* thetastartrow, 
	    const int* thetastartcol, 
	    const double* astartdata, 
	    const int* astartrow, 
	    const int* astartcol,
	    const double* bstartdata, 
	    const int* bstartrow, 
	    const int* bstartcol,
	    const double* t0,
	    const double* T0,	       
	    const double* ab0data, 
	    const int* ab0row, 
	    const int* ab0col,
	    const double* AB0data, 
	    const int* AB0row, 
	    const int* AB0col,
	    const double* thetaeqdata, 
	    const int* thetaeqrow, 
	    const int* thetaeqcol,
	    const double* thetaineqdata, 
	    const int* thetaineqrow, 
	    const int* thetaineqcol,
	    const int* storei,
	    const int* storea
	    ) {

    using namespace SCYTHE; //Added by Matthew S. Fasman on 11/04/2004
    
    // put together matrices
    const Matrix<int> X = r2scythe(*Xrow, *Xcol, Xdata);
    Matrix<double> theta = r2scythe(*thetastartrow, *thetastartcol, 
				    thetastartdata);
    Matrix<double> alpha = r2scythe(*astartrow, *astartcol, 
				    astartdata);
    Matrix<double> beta = r2scythe(*bstartrow, *bstartcol, 
				    bstartdata);
    const Matrix<double> ab0 = r2scythe(*ab0row, *ab0col, ab0data);
    const Matrix<double> AB0 = r2scythe(*AB0row, *AB0col, AB0data);
    const Matrix<double> theta_eq = r2scythe(*thetaeqrow, 
					     *thetaeqcol, 
					     thetaeqdata);
    const Matrix<double> theta_ineq = r2scythe(*thetaineqrow, 
					       *thetaineqcol, 
					       thetaineqdata);

    // initialize rng stream
    rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);
    
    // constants
    const int J = X.rows();  // number of subjects (justices, legislators)
    const int K = X.cols();  // number of items (cases, roll calls)
    const int tot_iter = *burnin + *mcmc;  
    const int nsamp = *mcmc / *thin;
    

    // storage matrices (row major order)
    Matrix<double> theta_store;
    Matrix<double> eta_store;
    if (*storea == 1){
      theta_store = Matrix<double>(nsamp, J);
    }
    if (*storei == 1){
      eta_store = Matrix<double>(nsamp, K*2);
    }

    // starting values 
    Matrix<double> eta   = cbind(alpha, beta);
    Matrix<double> Z     = Matrix<double>(J,K);
    
    int count = 0;
    // MCMC sampling occurs in this for loop
    for (int iter=0; iter < tot_iter; ++iter){
      
      // sample latent utilities (Z)
      irt_Z_update1(Z, X, theta, eta, stream);

      // sample item (case, bill) parameters (eta)
      irt_eta_update1(eta, Z, theta, ab0, AB0, stream);

      // sample ability (ideal points) (theta)
      irt_theta_update1(theta, Z, eta, *t0, *T0, theta_eq,
			theta_ineq, stream);

      // print results to screen
      if (*verbose > 0 && iter % *verbose == 0){
	Rprintf("\n\nMCMCirt1d iteration %i of %i \n",
		(iter+1), tot_iter);
	//Rprintf("theta = \n");
	//for (int j=0; j<J; ++j)
	//  Rprintf("%10.5f\n", theta[j]);    
      }
      
      // store results
      if ((iter >= burnin[0]) && ((iter % thin[0]==0))) {
	
	// store ideal points
	if (*storea == 1){
	  for (int l=0; l<J; ++l)
	    theta_store(count, l) = theta[l];
	}
	
	// store bill parameters
	if (*storei == 1){
	  for (int l=0; l<K*2; ++l)
	    eta_store(count, l) = eta[l];
	}
	count++;	
      }
      
      R_CheckUserInterrupt(); // allow user interrupts
      
    } // end Gibbs loop
    
    delete stream; // clean up random number stream
    
    // return output
    Matrix<double> output;
    if(*storei == 0 && *storea == 1) {
      output = theta_store;
    }
    else if (*storei == 1 && *storea == 0){
      output = eta_store;
    }
    else {
      output = cbind(theta_store, eta_store);
    }
    
    const int size = *samplerow * *samplecol;
    for (int i=0; i<size; ++i)
      sampledata[i] = output[i];

  }
}


