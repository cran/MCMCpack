//////////////////////////////////////////////////////////////////////////
// cMCMCfactanal.cc is C++ code to estimate a factor analysis model
//
// Andrew D. Martin
// Dept. of Political Science
// Washington University in St. Louis
// admartin@wustl.edu
//
// Kevin M. Quinn
// Dept. of Government
// Harvard University
// kevin_quinn@harvard.edu
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// revised version of older MCMCfactanal 5/11/2004 KQ
// updated to new verion of scythe 7/25/2004 ADM
// updated to Scythe 1.0.X 7/10/2007 ADM
// finished update to Scythe 1.0.X 7/30/2007 KQ
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////


#ifndef CMCMCFACTANAL_CC
#define CMCMCFACTANAL_CC

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

template <typename RNGTYPE>
void MCMCfactanal_impl (rng<RNGTYPE>& stream, const Matrix<>& X,
			Matrix<>& Lambda,
			Matrix<>& Psi, Matrix<>& Psi_inv,
			const Matrix<>& Lambda_eq,
			const Matrix<>&  Lambda_ineq,
			const Matrix<>& Lambda_prior_mean,
			const Matrix<>& Lambda_prior_prec,
			const Matrix<>& a0, const Matrix<>& b0,
			unsigned int burnin, unsigned int mcmc,
			unsigned int thin,
			unsigned int verbose,
			unsigned int storescores, Matrix<>& result) {

  // constants
  const unsigned int K = X.cols();  // number of manifest variables
  const unsigned int N = X.rows();  // number of observations
  const unsigned int D = Lambda.cols();  // number of factors
  const unsigned int tot_iter = burnin + mcmc;
  const unsigned int nsamp = mcmc / thin;
  const Matrix<> I = eye(D);
  const Matrix<> Lambda_free_indic = Matrix<>(K, D);
  for (unsigned int i=0; i<(K*D); ++i){
    if (Lambda_eq[i] == -999) Lambda_free_indic[i] = 1.0;
  }


  // starting value for phi
  Matrix<> phi = Matrix<>(N,D);

  // storage matrices (row major order)
  Matrix<> Lambda_store = Matrix<>(nsamp, K*D);
  Matrix<> Psi_store = Matrix<>(nsamp, K);
  Matrix<> phi_store;
  if (storescores==1){
    phi_store = Matrix<double>(nsamp, N*D);
  }


  unsigned int count    = 0;
  // sampling begins here
  for (unsigned int iter=0; iter < tot_iter; ++iter){

    // sample phi
    NormNormfactanal_phi_draw(phi, I, Lambda, Psi_inv,
			      X, N, D, stream);


    // sample Lambda
    NormNormfactanal_Lambda_draw(Lambda, Lambda_free_indic,
				 Lambda_prior_mean,
				 Lambda_prior_prec,
				 phi, X, Psi_inv,
				 Lambda_ineq,
				 D, K, stream);


    // sample Psi
    NormIGfactanal_Psi_draw(Psi, X, phi, Lambda,
			    a0, b0, K, N, stream);


    for (unsigned int i=0; i<K; ++i)
      Psi_inv(i,i) = 1.0 / Psi(i,i);



    // print results to screen
    if (verbose > 0 && iter % verbose == 0){
      Rprintf("\n\nMCMCfactanal iteration %i of %i \n", (iter+1), tot_iter);
      Rprintf("Lambda = \n");
      for (unsigned int i=0; i<K; ++i){
	for (unsigned int j=0; j<D; ++j){
	  Rprintf("%10.5f", Lambda(i,j));
	}
	Rprintf("\n");
      }
      Rprintf("diag(Psi) = \n");
      for (unsigned int i=0; i<K; ++i){
	Rprintf("%10.5f", Psi(i,i));
      }
      Rprintf("\n");
    }

    // store results
    if ((iter % thin)==0 && iter >= burnin ) {
      // store Lambda
      //Matrix<> Lambda_store_vec = Lambda.resize(1, K*D, true);
      //for (int l=0; l<K*D; ++l)
      //Lambda_store(count, l) = Lambda(l);
      rmview(Lambda_store(count, _)) = Lambda;

      // store Psi
      for (unsigned int i=0; i<K; ++i)
	Psi_store(count, i) = Psi(i,i);
      // stop phi
      if (storescores==1){
	//Matrix<> phi_store_vec = phi.resize(1, N*D, true);
	//for (int l=0; l<N*D; ++l)
	// phi_store(count, l) = phi(l);
	rmview(phi_store(count, _)) = phi;
      }
      count++;
    }

    // allow user interrupts
    R_CheckUserInterrupt();

  } // end Gibbs loop

  // return output
  //Matrix<double> result;
  result = cbind(Lambda_store, Psi_store);
  if(storescores == 1) {
    result = cbind(result, phi_store);
  }

}


extern "C" {

  void cMCMCfactanal(double *sampledata, const int *samplerow,
		    const int *samplecol, const double *Xdata,
		    const int *Xrow, const int *Xcol, const int *burnin,
		    const int *mcmc, const int *thin, const int *uselecuyer,
		    const int *seedarray, const int *lecuyerstream,
		    const int *verbose, const double *Lambdadata,
		    const int *Lambdarow, const int *Lambdacol,
		    const double *Psidata, const int *Psirow,
		    const int *Psicol, const double *Lameqdata,
		    const int *Lameqrow, const int *Lameqcol,
		    const double *Lamineqdata, const int *Lamineqrow,
		    const int *Lamineqcol, const double *Lampmeandata,
		    const int *Lampmeanrow, const int *Lampmeancol,
		    const double *Lampprecdata, const int *Lampprecrow,
		    const int *Lamppreccol, const double *a0data,
		    const int *a0row, const int *a0col,
		    const double *b0data, const int *b0row,
		    const int *b0col, const int *storescores) {

    // pull together Matrix objects
    const Matrix <> X(*Xrow, *Xcol, Xdata);
    Matrix <> Lambda(*Lambdarow, *Lambdacol, Lambdadata);
    Matrix <> Psi(*Psirow, *Psicol, Psidata);
    Matrix <> Psi_inv = invpd(Psi);
    const Matrix <> Lambda_eq(*Lameqrow, *Lameqcol,
			      Lameqdata);
    const Matrix <> Lambda_ineq(*Lamineqrow, *Lamineqcol,
				Lamineqdata);
    const Matrix <> Lambda_prior_mean(*Lampmeanrow,
				      *Lampmeancol,
				      Lampmeandata);
    const Matrix <> Lambda_prior_prec(*Lampprecrow,
				      *Lamppreccol,
				      Lampprecdata);
    const Matrix <> a0(*a0row, *a0col, a0data);
    const Matrix <> b0(*b0row, *b0col, b0data);

    Matrix<> storagematrix;

    MCMCPACK_PASSRNG2MODEL(MCMCfactanal_impl, X, Lambda, Psi, Psi_inv,
			   Lambda_eq, Lambda_ineq, Lambda_prior_mean,
			   Lambda_prior_prec,
			   a0, b0, *burnin, *mcmc, *thin,
			   *verbose, *storescores, storagematrix);

    const unsigned int size = *samplerow * *samplecol;
    for (unsigned int i = 0; i < size; ++i)
      sampledata[i] = storagematrix(i);
  }
}

#endif

