//////////////////////////////////////////////////////////////////////////
// HMMpanelFE.cc is C++ code
//
// Jong Hee Park
// Department of Political Science and International Relations
// Seoul National University
// jongheepark@snu.ac.kr

// Written 11/19/2008
// Modified 09/20/2009
//////////////////////////////////////////////////////////////////////////
#ifndef CHMMPANELFE_CC
#define CHMMPANELFE_CC

#include<vector>
#include<algorithm>

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
#include "lapack.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;


// used to access 1d arrays holding R matrices like a 2d array
#define M(ROW,COL,NROWS) (COL*NROWS+ROW)

template <typename RNGTYPE>
// For better identification, first two and last two states are constrained to be 1 and ns
// 10/31/2011 JHP
Matrix<int> hetero_state_sampler(rng<RNGTYPE>& stream,
				 const int m,
				 const int ntime_s,
				 const Matrix<>& Y,
				 const Matrix<>& delta,
				 const Matrix<>& Sigma,
				 const Matrix<>& P){
  const int ns = m + 1;
  const int ntime = ntime_s;

  Matrix<> F(ntime, ns);
  Matrix<> pr1(ns, 1);
  pr1[0] = 1;
  Matrix<> py(ns, 1);
  Matrix<> pstyt1(ns, 1);

  for (int tt=0; tt<ntime ; ++tt){
    for (int j = 0; j< ns; ++j){
      py[j]  =  dnorm(Y[tt], delta[j], sqrt(Sigma[j]));
    }
    if (tt==0) pstyt1 = pr1;
    else {
      pstyt1 =  ::t(F(tt-1,_)*P); // make it an ns by 1 matrix
    }

    Matrix<> unnorm_pstyt = pstyt1%py;
      const Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt); // pstyt = Pr(st|Yt)
      for (int j=0; j<ns ; ++j) {
	F(tt,j) = pstyt(j);
      }

  }// end of F matrix filtering
  Matrix<int> state(ntime, 1);
  Matrix<> ps = Matrix<>(ntime, ns);
  ps(ntime-1,_) = F(ntime-1,_);
  state(ntime-1) = ns;

  Matrix<> pstyn = Matrix<>(ns, 1);
  double pone = 0.0;
  int tt = ntime-2;

  while (tt >= 0){
    int st = state(tt+1);
    Matrix<> Pst_1 = ::t(P(_,st-1));
    Matrix<> unnorm_pstyn = F(tt,_)%Pst_1;
    pstyn = unnorm_pstyn/sum(unnorm_pstyn);

    if (st==1){
      state(tt) = 1;
    }
    else{
      pone = pstyn(st-2);
      if(stream.runif() < pone) state(tt) = st-1;
      else state(tt) = st;
    }
    ps(tt,_) = pstyn;
    --tt;
  }// end of while loop

  return state;

} // end of state sampler


template <typename RNGTYPE>
void HMMpanelFE_impl(rng<RNGTYPE>& stream,
		      int nsubj,  int ntime,
		      int mmax,  int mmin,
		     const Matrix<int>& mvector,
		      int totalstates,
		     const Matrix<>& Y, const Matrix<>& X,
		     const Matrix<int>& subjectid,
		      int  burnin,  int  mcmc,
		      int thin,  int  verbose,
		     Matrix<>& beta, double sigma2, Matrix<>& deltastart,
		     const Matrix<>& b0, const Matrix<>& B0,
		     const double delta0, const double Delta0,
		     const double c0, const double d0,
		     const Matrix<>& P0data,
		     const Matrix<>& Pstart,
		     const Matrix<>& subject_groupinfo,
		     Matrix<>& betastorage,
		     Matrix<>& statestorage,
		     Matrix<>& deltastorage,
		     Matrix<>& sigmastorage){

  // redefine constants
  const  int K = X.cols();; // ncol(X)
  const int NOBS = Y.rows();

  const int tot_iter = burnin + mcmc;
  vector< vector <double> > P0_holder;
  vector< vector <double> > P_holder;
  int count = 0;
  for ( int s=0; s< nsubj; ++s){
    const int nms = mvector[s] + 1;
    vector <double> P0mat;
    vector <double> Pmat;
    for(int ii=0; ii<(nms*nms); ++ii){
      P0mat.push_back(P0data(count + ii));
      Pmat.push_back(Pstart(count + ii));
    }
    count = count + nms*nms;

    P0_holder.push_back(P0mat);
    P_holder.push_back(Pmat);
  }

  vector< vector <double> > delta_holder;
  vector< vector <double> > sigma2_holder;
  vector< vector <int> > nstate;
  count = 0;
  for ( int s=0; s< nsubj; ++s){
    int nms = mvector[s] + 1;
    int ntime_s = subject_groupinfo(s, 2);
    vector <double> deltamat;
    vector <double> sigmamat;
    vector <int> nstatemat;
    for(int ii=0; ii<nms ;++ii){
      deltamat.push_back(0);
      nstatemat.push_back(ntime_s);
      sigmamat.push_back(sigma2);
    }
    count = count + nms;
    delta_holder.push_back(deltamat);
    nstate.push_back(nstatemat);
    sigma2_holder.push_back(sigmamat);
  }

  int *nobsk = new int[nsubj];
  for ( int k=0; k<nsubj; k++) {
    nobsk[k]=0;
    for ( int n=0; n<NOBS; n++) {
      if (subjectid[n]==k+1) {
	nobsk[k]+=1;
      }
    }
  }

  int **posk_arr = new int*[nsubj];
  for ( int k=0; k<nsubj; k++) {
    posk_arr[k] = new int[nobsk[k]];
    int repk=0;
    for ( int n=0; n<NOBS; n++) {
      if (subjectid[n]==k+1) {
	posk_arr[k][repk]=n;
	repk++;
      }
    }
  }

  Matrix<double> *Yk_arr = new Matrix<double>[nsubj];
  Matrix<double> *Xk_arr = new Matrix<double>[nsubj];
  for( int k=0; k<nsubj; k++) {
    Xk_arr[k] = Matrix<double>(nobsk[k], K);
    Yk_arr[k] = Matrix<double>(nobsk[k], 1);
    for (int l=0; l<nobsk[k]; l++) {
      for ( int p=0; p< K; p++) {
	Xk_arr[k](l, p) = X[p*NOBS + posk_arr[k][l]];
      }
      Yk_arr[k](l,0) = Y[posk_arr[k][l]];
    }
  }

  Matrix<double> *tXk_arr = new Matrix<double>[nsubj];
  for( int k=0; k<nsubj; k++) {
    tXk_arr[k] = ::t(Xk_arr[k]);
  }

  /////////////////////////////////////////////////
  // initialize newY for the first loop only
  /////////////////////////////////////////////////
  Matrix<>* newY = new Matrix<>[nsubj]; // newY = Y - delta
  Matrix<>* Yres = new Matrix<>[nsubj]; // Yres = Y - Xbeta
  for ( int s=0; s<nsubj; ++s){
    newY[s] = Yk_arr[s] - delta0;
    Yres[s] = Yk_arr[s];
  }

  // MCMC iterations start here
   int sampcount = 0;
   for ( int iter=0; iter < tot_iter; ++iter){
     double delta_sum = 0;
     for ( int s=0; s< nsubj; ++s) {
       int ntime_s = subject_groupinfo(s, 2);
       Matrix<> Zk_arr(ntime_s, 1);
       for (int tt=0; tt< ntime_s; ++tt) {
	 Zk_arr(tt) = 1;
       }
       Yres[s] = Yk_arr[s] - Xk_arr[s]*beta;

       if(mvector[s]==0){
	 double Dn = 1/(Delta0 + (double)ntime_s/sigma2_holder[s][0]);
	 double dn = Dn*(Delta0*delta0 + sum(Yres[s])/sigma2_holder[s][0]);
	 delta_holder[s][0] = stream.rnorm(dn, sqrt(Dn));
	 delta_sum = delta_sum + delta_holder[s][0];
	 newY[s] = Yk_arr[s] - Zk_arr*delta_holder[s][0];

	 // Sample sigma
	 double shape = (c0 + (double)ntime_s)/2;
	 const Matrix<> SSE = crossprod (newY[s]);
	 double scale =(d0 + SSE[0])/2;
	 sigma2_holder[s][0] = 1/stream.rgamma(shape, scale);
       }
       else {
	 const int nscur = mvector[s] + 1;
	 Matrix<> P(nscur, nscur);
	 Matrix<> P0(nscur, nscur);
	 Matrix<> delta(nscur, 1);
	 Matrix<> Sigma(nscur, 1);

	 for (int i=0;i<(nscur*nscur); ++i){
	   P0[i] = P0_holder[s][i];
	   P[i] = P_holder[s][i];
	 }

	 for (int i=0;i<nscur; ++i){
	   delta[i] = delta_holder[s][i];
	   Sigma[i] = sigma2_holder[s][i];
	 }
	 // Sample s
	 Matrix<int> state_s = hetero_state_sampler(stream, mvector[s], ntime_s, Yres[s], delta, Sigma, P);
	 // Sample delta and Sigma
	 int delta_count = 0;
	 for (int j = 0; j <nscur ; ++j){
	   nstate[s][j] = 0;
	   for (int i = 0; i<ntime_s; ++i){
	     if (state_s[i] == (j+1)) {
	       nstate[s][j] = nstate[s][j] + 1;
	     }// end of if
	   }// end of int i<n

	   delta_count = delta_count + nstate[s][j];
	   // sample delta
	   Matrix<> yj = Yres[s]((delta_count - nstate[s][j]), 0, (delta_count - 1), 0);
	   Matrix<> Yj = Yk_arr[s]((delta_count - nstate[s][j]), 0, (delta_count - 1), 0);
	   double Dn = 1/(Delta0 + (double)nstate[s][j]/sigma2_holder[s][j]);
	   double dn = Dn*(Delta0*delta0 + sum(yj)/sigma2_holder[s][j]);
	   delta_holder[s][j] = stream.rnorm(dn, sqrt(Dn));
	   delta_sum = delta_sum + delta_holder[s][j];
	   newY[s]((delta_count - nstate[s][j]), 0, (delta_count - 1), 0) = Yj - delta_holder[s][j];

	   // Sample sigma
	   double shape = (c0 + (double)nstate[s][j])/2;
	   const Matrix<> SSE = crossprod (newY[s]((delta_count - nstate[s][j]), 0, (delta_count - 1), 0));
	   double scale =(d0 + SSE[0])/2;
	   sigma2_holder[s][j] = 1/stream.rgamma(shape, scale);

	 }
	 // assure that there is no label switching problem
	 // the code needs to be added here

	 // Sample P
	 double shape1 = 0;
	 double shape2 = 0;

	 for (int j =0; j<(nscur-1); ++j){
	   shape1 =  std::abs(P0(j,j) + nstate[s][j] - 1);
	   shape2 =  P0(j,j+1) + 1; //
	   P(j,j) = stream.rbeta(shape1, shape2);
	   P(j,j+1) = 1 - P(j,j);
	 }
	 P(mvector[s], mvector[s]) = 1; //no jump at the last state

	 for(int ii=0; ii<(nscur*nscur) ;++ii) {
	   P_holder[s][ii] = P[ii];
	 }
       }//end of else (mvector!=0)

     }// end of subject specific looping
     // Sample beta
     Matrix<> XVX(K, K);
     Matrix<> XVY(K, 1);
     for( int s = 0; s<nsubj; ++s) {
       int ntime_s = subject_groupinfo(s, 2);
       int delta_count = 0;
       Matrix<> Vi = eye(ntime_s);
       for (int j = 0; j <(mvector[s] + 1); ++j){
	 delta_count = delta_count + nstate[s][j];
	 for(int i = (delta_count - nstate[s][j]); i<delta_count; ++i) {
	   Vi(i,i) = 1/sigma2_holder[s][j];
	 }
       }
       XVX = XVX + tXk_arr[s] * Vi * Xk_arr[s];
       XVY = XVY + tXk_arr[s] * Vi * newY[s];
     }
     Matrix<> beta_var = invpd(B0 + XVX);
     Matrix<> beta_mean = beta_var*(B0*b0 + XVY);
     beta = stream.rmvnorm(beta_mean, beta_var);
     // STORE
     if (iter >= burnin && ((iter % thin) == 0)) {
      for( int j=0;j<K; ++j) {
	betastorage(sampcount,j) = beta(j);
      }
      int count = 0;
      for (int s=0; s<nsubj; ++s){
	for (int j=0; j<(mvector[s] + 1); ++j){
	  sigmastorage(sampcount, count) = sigma2_holder[s][j];
	  deltastorage(sampcount, count) = delta_holder[s][j];
	  double nstat = nstate[s][j];
	  statestorage(sampcount, count) = nstat;
	  ++ count;
	}
      }
      ++sampcount;
     }

     // REPORT
     if(verbose > 0 && iter % verbose == 0){
       Rprintf("\n ----------------------------------------------------------------------- ");
       Rprintf("\n\n HMMpanelFE %i of %i \n", iter, tot_iter);
       Rprintf("\n beta = \n");
       for(int i=0;i<K; ++i) {
	 Rprintf("%10.5f\n", beta(i));
       }
     }
   }// END of MCMC loop

   delete [] Yres;
   delete [] tXk_arr;
   delete [] Xk_arr;
   delete [] Yk_arr;
   delete [] newY;
}


extern "C" {
  void cHMMpanelFE(double *deltadraws,  double* sigmadraws,
		  double *statedraws, //const int* statecol,
		  double* betadraws, const int* betarow, const int* betacol,
		  const int* totalstates,
		  const int* nsubj, const int* ntime, const int* nobs,
		  const int* subjectid,
		  const int* m,
		  const int* mmax, const int* mmin,
		  const double* Ydata, const int* Yrow, const int* Ycol,
		  const double* Xdata, const int* Xrow, const int* Xcol,
		  const int* burnin, const int* mcmc, const int* thin, const int* verbose,
		  const int *uselecuyer, const int *seedarray, const int *lecuyerstream,
		  const double* betastartdata, const double* sigma2start,
		  const double* deltastartdata, const int* deltastartrow,
		  const double* b0data, const double* B0data,
		  const double* delta0, const double* Delta0,
		  const double* c0, const double* d0,
		  const double* P0data, const int* P0row,
		  const double* Pstartdata,
		  const double* subject_groupinfodata){


    // pull together Matrix objects
    Matrix<> Y(*Yrow, *Ycol, Ydata);
    Matrix<> X(*Xrow, *Xcol, Xdata);
    Matrix<> betastart(*Xcol, 1, betastartdata);
    Matrix<> deltastart(*deltastartrow, 1, deltastartdata);
    Matrix<> b0(*Xcol, 1, b0data);
    Matrix<> B0(*Xcol, *Xcol, B0data);
    Matrix<int> subjectid_mat(*nobs, 1, subjectid);
    Matrix<> subject_groupinfo(*nsubj, 3, subject_groupinfodata);
    Matrix<> P0(*P0row, 1, P0data);
    Matrix<> Pstart(*P0row, 1, Pstartdata);
    Matrix<int> mvector(*nsubj, 1, m);

    Matrix<> betastorage(*betarow, *betacol);
    Matrix<> sigmastorage(*betarow,  *totalstates);
    Matrix<> deltastorage(*betarow,  *totalstates);
    Matrix<> statestorage(*betarow,  *totalstates);

    MCMCPACK_PASSRNG2MODEL(HMMpanelFE_impl,
			   *nsubj, *ntime, *mmax, *mmin, mvector,
			   *totalstates, Y, X,
			   subjectid_mat,
			   *burnin, *mcmc, *thin, *verbose,
			   betastart, *sigma2start, deltastart,
			   b0, B0, *delta0, *Delta0, *c0, *d0,
			   P0, Pstart, subject_groupinfo,
			   betastorage, statestorage, deltastorage, sigmastorage);
     int deltasize = *betarow * *totalstates;
    for ( int i=0; i < deltasize; ++i){
      deltadraws[i] = deltastorage(i);
      sigmadraws[i] = sigmastorage(i);
      statedraws[i] = statestorage(i);
    }
     int betasize = *betarow * *betacol;
    for ( int i=0; i < betasize; ++i){
      betadraws[i] = betastorage(i);
    }
  }// end of cHMMpanelFE

}// end of extern C

#endif /* CHMMPANELFE_CC  */
