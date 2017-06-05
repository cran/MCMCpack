//////////////////////////////////////////////////////////////////////////
// cMCMCirt1d.cc is C++ code to estimate a one-dimensional item response
// theory model.
//
// ADM and KQ 1/15/2003
// ADM 7/28/2004 [updated to new Scythe version]
// completely rewritten and optimized for the 1-d case 8/2/2004 KQ
// storage changed to save memory KQ 1/27/2006
// DBP 7/3/07 [ported to scythe 1.0.x]
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////

#ifndef CMCMCIRT1D_CC
#define CMCMCIRT1D_CC

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

/* cMCMCirt1d implementation. */
template <typename RNGTYPE>
void MCMCirt1d_impl (rng<RNGTYPE>& stream, const Matrix<int>& X,
    Matrix<>& theta, Matrix<>& eta, const Matrix<>& ab0,
    const Matrix<>& AB0, const Matrix<>& theta_eq,
    const Matrix<>& theta_ineq, double t0, double T0, unsigned int burnin,
    unsigned int mcmc, unsigned int thin, unsigned int verbose,
    bool storea, bool storei, double* sampledata, unsigned int samplesize)
{
  // constants
  const unsigned int J = X.rows();  // # subjects (justices, legislators)
  const unsigned int K = X.cols();  // number of items (cases, roll calls)
  const unsigned int tot_iter = burnin + mcmc;
  const unsigned int nsamp = mcmc / thin;

  // storage matrices (col major order)
  Matrix<> theta_store;
  Matrix<> eta_store;
  if (storea)
    theta_store = Matrix<>(nsamp, J);

  if (storei)
    eta_store = Matrix<>(nsamp, K*2);

  // starting values
  Matrix<> Z(J, K);

  // pre-compute
  const Matrix<> AB0ab0 = AB0 * ab0;

  unsigned int count = 0;
  // MCMC sampling occurs in this for loop
  for (unsigned int iter = 0; iter < tot_iter; ++iter){

    // sample latent utilities (Z)
    irt_Z_update1(Z, X, theta, eta, stream);

    // sample item (case, bill) parameters (eta)
    irt_eta_update1(eta, Z, theta,  AB0, AB0ab0, stream);

    // sample ability (ideal points) (theta)
    irt_theta_update1(theta, Z, eta, t0, T0, theta_eq, theta_ineq, stream);

    // print results to screen
    if (verbose > 0 && iter % verbose == 0) {
      Rprintf("\n\nMCMCirt1d iteration %i of %i \n", (iter+1), tot_iter);
      //Rprintf("theta = \n");
      //for (int j=0; j<J; ++j)
      //Rprintf("%10.5f\n", theta[j]);
    }

    // store results
    if ((iter >= burnin) && ((iter % thin == 0))) {

      // store ideal points
      if (storea)
        theta_store(count, _) = theta;

      // store bill parameters
      if (storei)
        eta_store(count, _) = t(eta);
      count++;
    }

    R_CheckUserInterrupt(); // allow user interrupts

  } // end Gibbs loop

  // return output
  Matrix<> output;
  if(! storei && storea) {
    output = theta_store;
  } else if (storei && ! storea){
    output = eta_store;
  } else {
    output = cbind(theta_store, eta_store);
  }

  for (unsigned int i = 0; i < samplesize; ++i)
    sampledata[i] = output[i];

}

extern "C" {

  void
  cMCMCirt1d(double* sampledata, const int* samplerow, const int* samplecol,
	    const int* Xdata, const int* Xrow, const int* Xcol,
	    const int* burnin, const int* mcmc,  const int* thin,
	    const int *uselecuyer, const int *seedarray,
      const int *lecuyerstream,
	    const int* verbose, const double* thetastartdata,
	    const int* thetastartrow, const int* thetastartcol,
	    const double* astartdata, const int* astartrow, const int* astartcol,
	    const double* bstartdata, const int* bstartrow, const int* bstartcol,
	    const double* t0, const double* T0,	const double* ab0data,
	    const int* ab0row, const int* ab0col, const double* AB0data,
	    const int* AB0row, const int* AB0col, const double* thetaeqdata,
	    const int* thetaeqrow, const int* thetaeqcol,
	    const double* thetaineqdata, const int* thetaineqrow,
	    const int* thetaineqcol, const int* storei, const int* storea)
  {
    // put together matrices
    const Matrix<int> X(*Xrow, *Xcol, Xdata);
    Matrix<> theta(*thetastartrow, *thetastartcol, thetastartdata);
    Matrix<> alpha(*astartrow, *astartcol, astartdata);
    Matrix<> beta(*bstartrow, *bstartcol, bstartdata);
    const Matrix<> ab0(*ab0row, *ab0col, ab0data);
    const Matrix<> AB0(*AB0row, *AB0col, AB0data);
    const Matrix<> theta_eq(*thetaeqrow, *thetaeqcol, thetaeqdata);
    const Matrix<> theta_ineq(*thetaineqrow, *thetaineqcol, thetaineqdata);
    Matrix<> eta  = cbind(alpha, beta);
    const int samplesize = (*samplerow) * (*samplecol);

    MCMCPACK_PASSRNG2MODEL(MCMCirt1d_impl, X, theta, eta, ab0,
        AB0, theta_eq, theta_ineq, *t0, *T0, *burnin, *mcmc, *thin,
        *verbose, *storea, *storei, sampledata, samplesize);
  }
}

#endif
