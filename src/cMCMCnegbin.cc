//////////////////////////////////////////////////////////////////////////
// cMCMCnegbin.cc is C++ code to estimate a  Negative
// Binomial regression model
//
// Matthew Blackwell
// Department of Government
// Harvard University
// mblackwell@gov.harvard.edu

// Written 05/19/2017
//////////////////////////////////////////////////////////////////////////

#ifndef CMCMCNEGBIN_CC
#define CMCMCNEGBIN_CC

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

#include "MCMCnbutil.h"

#include <R.h>      
#include <R_ext/Utils.h> 

using namespace std;
using namespace scythe;


//////////////////////////////////////////// 
// MCMCnegbinReg implementation.  
//////////////////////////////////////////// 
template <typename RNGTYPE>
void MCMCnegbinReg_impl(rng<RNGTYPE>& stream, 
                        double *betaout, 
                        double *nuout,
                        double *rhoout,
                        double *tau1out,
                        double *tau2out,
                        int *comp1out,
                        int *comp2out,
                        double *sr1out,
                        double *sr2out,
                        double *mr1out,
                        double *mr2out,
                        double *rhosizes,
                        const double *Ydata, 
                        const int *Yrow, 
                        const int *Ycol, 
                        const double *Xdata, 
                        const int *Xrow, 
                        const int *Xcol, 
                        const int *burnin, 
                        const int *mcmc, 
                        const int *thin, 
                        const int *verbose, 
                        const double *betastart, 
                        const double *nustart,
                        const double *rhostart,
                        const double *tau1start, 
                        const double *tau2start, 
                        const double *component1start,
                        const double *e,
                        const double *f,
                        const double *g,
                        const double *rhostepdata,
                        const double *b0data, 
                        const double *B0data, 
                        double *logmarglikeholder, 
                        double *loglikeholder, 
                        const int *chib)
{
  const Matrix <> Y(*Yrow, *Ycol, Ydata);  
  const Matrix <> X(*Xrow, *Xcol, Xdata);
  const int tot_iter = *burnin + *mcmc;  
  const int nstore = *mcmc / *thin;     
  const int n = Y.rows();
  const int k = X.cols();
  const int max_comp = 10;
  const Matrix <> b0(k, 1, b0data);
  const Matrix <> B0(k, k, B0data);   
  const Matrix <> B0inv = invpd(B0);   
  Matrix<> wr1(max_comp, 1);
  Matrix<> mr1(max_comp, 1);
  Matrix<> sr1(max_comp, 1);
  Matrix<> wr2(n, max_comp);
  Matrix<> mr2(n, max_comp);
  Matrix<> sr2(n, max_comp);
  Matrix<> nr2(n, 1);

 
  Matrix <> nu(n, 1, nustart);
  double rho = *rhostart;
  double step_out = *rhostepdata;
  Matrix <> beta(k, 1, betastart);
  Matrix <> tau1(n, 1, tau1start);
  Matrix <> tau2(n, 1, tau2start);
  Matrix <> component1(n, 1, component1start);
  Matrix <> component2(n, 1);
  
  Matrix<> beta_store(nstore, k);
  Matrix<int> component1_store(nstore, n);
  Matrix<int> component2_store(nstore, n);
  Matrix<> tau1_store(nstore, n);
  Matrix<> tau2_store(nstore, n);
  Matrix<> sr1_store(nstore, n);
  Matrix<> sr2_store(nstore, n);
  Matrix<> mr1_store(nstore, n);
  Matrix<> mr2_store(nstore, n);
  Matrix<> nu_store(nstore, n);
  Matrix<> rho_store(nstore, 1);
  Matrix<> rhostep;

  init_aux(stream, Y, wr1, mr1, sr1, wr2, mr2, sr2, nr2, component2);
  int nplus = 0;
  for (int t = 0; t < n; t++) {
    int yt = (int) Y[t];
    if (yt > 0) {
      nplus = nplus + 1;
    }
  }

  Matrix<> Xplus(nplus, k);
  int xp_count = 0;
  for (int t = 0; t < nplus; t++) {
    int yt = (int) Y[t];
    Matrix<> xt = X(t, _);
    if (yt > 0) {
      Xplus(xp_count, _) = xt;
      xp_count++;
    }
  }
  //MCMC loop
  int count = 0;
  int debug = 0;
  for (int iter = 0; iter < tot_iter; ++iter){
    if (debug > 0) {
      Rprintf("\niter: %i\n-----\n", iter);
    }

    if (debug > 0) {
      Rprintf("Sampling rho...\n");
    }

    
    //////////////////////
    // 1. Sample rho
    //////////////////////
    Matrix<> lambda = exp(X * beta);    
    rhostep = rho_slice_sampler(stream, Y, lambda, rho, step_out, *g, *e, *f);
    rho = rhostep[0];
    

    if (debug > 0) {
      Rprintf("Sampling nu...\n");
    }
    //////////////////////
    // 2. Sample nu
    //////////////////////
    for (int t = 0; t < n; t++) {
      int yt = (int) Y[t];
      nu[t] = stream.rgamma(rho + yt, rho + lambda[t]);
    }

    if (debug > 0) {
      Rprintf("Sampling tau...\n");
    }
    //////////////////////
    // 4. Sample tau 
    //////////////////////
    Matrix<> TAUout;
    for (int t = 0; t < n; t++) {
      int yt = (int) Y[t];
      double mut = exp(log(nu[t]) + log(lambda[t]));
      
      TAUout = tau_comp_sampler(stream, yt, mut, wr1, mr1, sr1, 
                                wr2(t,_), mr2(t,_), sr2(t,_), nr2[t]);
      tau1[t] = TAUout[0];
      tau2[t] = TAUout[1];
      component1[t] = TAUout[2];
      component2[t] = TAUout[3];
    }


    if (debug > 0) {
      Rprintf("Sampling beta...\n");
    }
    //////////////////////
    // 2. Sample beta 
    //////////////////////
    Matrix<> y_tilde(n,1);
    Matrix<> Sigma_inv_sum(n, 1);
    Matrix<> yp_tilde(n,1);
    Matrix<> Sigma_plus_inv_sum(n, 1);
    Matrix<> sr1_hold(n,1);
    Matrix<> sr2_hold(n,1);
    Matrix<> mr1_hold(n,1);
    Matrix<> mr2_hold(n,1);
    for (int t = 0; t<n ; ++t) {
      int yt = (int) Y[t]; 
      int comp1 = (int) component1[t];
        if (yt > 0) {
          int comp2 = (int) component2[t];
          sr2_hold[t] = sr2(t, comp2 - 1);
          mr2_hold[t] = mr2(t, comp2 - 1);
          Sigma_plus_inv_sum[t] = 1/sqrt(sr2(t, comp2 - 1));
          yp_tilde[t] = (-log(tau2[t]) - log(nu[t]) - mr2(t, comp2 - 1))/sqrt(sr2(t, comp2-1));
          
        }
        sr1_hold[t] = sr1[comp1 - 1];
        mr1_hold[t] = mr1[comp1 - 1];          
        Sigma_inv_sum[t] = 1/sqrt(sr1[comp1 - 1]);
        y_tilde[t] = (-log(tau1[t]) - log(nu[t])  - mr1[comp1 - 1])/sqrt(sr1[comp1 - 1]);

    }
    Matrix<> yjp = rbind(y_tilde, selif(yp_tilde, Y > 0));
    Matrix<> Xjp = rbind(X, selif(X, Y > 0));
    Matrix<> wip = rbind(Sigma_inv_sum, selif(Sigma_plus_inv_sum, Y > 0));
    Matrix<> Xwj(Xjp.rows(), k);
    for (unsigned int h = 0; h<Xjp.rows(); ++h){
      Xwj(h, _) = Xjp(h,_)*wip[h];
    }
    Matrix<> Bn = invpd(B0 + ::t(Xwj)*Xwj);
    Matrix<> bn = Bn*gaxpy(B0, b0, ::t(Xwj)*yjp);
    if (k == 1) {
      beta[0] = stream.rnorm(bn[0], sqrt(Bn[0]));
    } else {
      beta = stream.rmvnorm(bn, Bn);
    }

    if (debug > 0) {
      Rprintf("Sampling P...\n");
    }


    if (iter >= *burnin && ((iter % *thin)==0)){
      beta_store(count,_) = ::t(beta);
      tau1_store(count,_) = tau1;
      tau2_store(count,_) = tau2;
      component1_store(count,_) = component1;
      component2_store(count,_) = component2;
      sr1_store(count,_) = sr1_hold;
      sr2_store(count,_) = sr2_hold;
      mr1_store(count,_) = mr1_hold;
      mr2_store(count,_) = mr2_hold;
      nu_store(count,_) = nu;
      rho_store(count, _) = rho;
      ++count; 
    }   

    if(*verbose > 0 && iter % *verbose == 0){
      Rprintf("\n\n MCMCnegbin iteration %i of %i \n", (iter+1), tot_iter);
      Rprintf("rho is %10.5f\n", rho);
      for (int j = 0; j<k; ++j){
        Rprintf("beta(%i) is %10.5f\n", j+1, beta[j]);
      }
    }
    
  }// end MCMC loop 
  
  
  
  //////////////////////////////////////
  if (*chib == 1) {
  //////////////////////////////////////
    if(*verbose > 0){
      Rprintf("\nCalculating marginal likelihood...\n");
    }
    Matrix<> beta_st = ::t(meanc(beta_store));
    Matrix<> lambda = exp(X * beta_st);        
    double rho_st = exp(mean(log(rho_store)));

    //////////////////////
    // beta
    //////////////////////      
    Matrix<> density_beta(nstore, 1);  
    for (int iter = 0; iter<nstore; ++iter){
      Matrix<> y_tilde(n,1);
      Matrix<> Sigma_inv_sum(n, 1);
      Matrix<> yp_tilde(n,1);
      Matrix<> Sigma_plus_inv_sum(n, 1);
      for (int t = 0; t<n ; ++t) {
        int yt = (int) Y[t]; 
        int comp1 = (int) component1_store(iter, t);
        if (yt > 0) {
          int comp2 = (int) component2_store(iter, t);
          Sigma_plus_inv_sum[t] = 1/sqrt(sr2(t, comp2 - 1));
          yp_tilde[t] = (-log(tau2_store(iter,t)) - log(nu_store(iter,t)) - mr2(t, comp2 - 1))/sqrt(sr2(t, comp2-1));
        }
        Sigma_inv_sum[t] = 1/sqrt(sr1[comp1 - 1]);
        y_tilde[t] = (-log(tau1_store(iter,t)) - log(nu_store(iter, t))  - mr1[comp1 - 1])/sqrt(sr1[comp1 - 1]);
      }
      Matrix<> yjp = rbind(y_tilde, selif(yp_tilde, Y > 0));
      Matrix<> Xjp = rbind(X, selif(X, Y > 0));
      Matrix<> wip = rbind(Sigma_inv_sum, selif(Sigma_plus_inv_sum, Y > 0));
      Matrix<> Xwj(Xjp.rows(), k);
      for (unsigned int h = 0; h<Xjp.rows(); ++h){
        Xwj(h, _) = Xjp(h,_)*wip[h];
      }
      Matrix<> Bn = invpd(B0 + ::t(Xwj)*Xwj);
      Matrix<> bn = Bn*gaxpy(B0, b0, ::t(Xwj)*yjp);
      if (k == 1) {
        density_beta[iter] = exp(lndnorm(beta_st[0], bn[0], sqrt(Bn[0])));
      } else {
        density_beta[iter] = exp(lndmvn(beta_st, bn, Bn));
      }      
    }
    double pdf_beta = log(mean(density_beta));

    //////////////////////
    // rho/nu
    //////////////////////
    Matrix<> density_rho(nstore, 1);
    for (int iter = 0; iter < nstore; iter++) {

      rhostep = rho_slice_sampler(stream, Y, lambda, rho,
                                  step_out, *g, *e, *f);
      rho = rhostep[0];
      double L = rhostep[3];
      double R = rhostep[4];

      // We're going to assume that the distribution is unimodal and we can use the slice width.
      if (rho_st > L && rho_st < R) {
        density_rho[iter] = 1/(R-L);
      } else {
        density_rho[iter] = 0.0;
      }
    }
    
    double pdf_rho = log(mean(density_rho));
    
    
    //////////////////////
    // likelihood
    //////////////////////
    double loglike = 0;
    for (int t = 0; t < n ; ++t){
      int yt = (int) Y[t];
      loglike += lngammafn(rho_st + yt) - lngammafn(rho_st) - lngammafn(yt + 1);
      loglike += rho_st * log(rho_st) + yt*log(lambda[t]) - (rho_st + yt) * log(rho_st + lambda[t]);
    }
        
    //////////////////////
    // prior
    //////////////////////
    double density_beta_prior;
    double density_rho_prior;
    
    if (k == 1) {
      density_beta_prior = lndnorm(beta_st[0], b0[0], B0inv[0]); 
    } else {
      density_beta_prior = lndmvn(beta_st, b0, B0inv);
    }
    
    density_rho_prior = (*e - 1) * log(rho_st) - (*e + *f) * log(rho_st + *g);
    
    // compute marginal likelihood
    const double logprior = density_beta_prior + density_rho_prior; 
    const double logmarglike = (loglike + logprior) - (pdf_beta + pdf_rho);
    if (*verbose >0 ){
      Rprintf("\nlogmarglike = %10.5f\n", logmarglike);
      Rprintf("loglike = %10.5f\n", loglike);
      Rprintf("log_prior = %10.5f\n", logprior);
      Rprintf("log_beta = %10.5f\n", pdf_beta);
      Rprintf("log_rho = %10.5f\n", pdf_rho);
    }
    logmarglikeholder[0] = logmarglike;
    loglikeholder[0] = loglike;
    
  }
        
  R_CheckUserInterrupt();      
  
  for (int i = 0; i<(nstore*k); ++i){
    betaout[i] = beta_store[i]; 
  }
  for (int i = 0; i<(nstore*n); ++i){
    nuout[i] = nu_store[i];
    tau1out[i] = tau1_store[i];
    tau2out[i] = tau2_store[i];
    comp1out[i] = component1_store[i];
    comp2out[i] = component2_store[i];
    sr1out[i] = sr1_store[i];
    sr2out[i] = sr2_store[i];
    mr1out[i] = mr1_store[i];
    mr2out[i] = mr2_store[i];
  }
  for (int i = 0; i<nstore; ++i){
    rhoout[i] = rho_store[i];
  }
  *rhosizes = step_out;
}
 
extern "C" {
  void cMCMCnegbin(double *betaout, 
                   double *nuout,
                   double *rhoout,
                   double *tau1out,
                   double *tau2out,
                   int *comp1out,
                   int *comp2out,
                   double *sr1out,
                   double *sr2out,
                   double *mr1out,
                   double *mr2out,
                   double *rhosizes,
                   const double *Ydata, 
                   const int *Yrow, 
                   const int *Ycol, 
                   const double *Xdata,
                   const int *Xrow, 
                   const int *Xcol, 
                   const int *burnin, 
                   const int *mcmc, 
                   const int *thin, 
                   const int *verbose, 
                   const double *betastart, 
                   const double *nustart,
                   const double *rhostart,
                   const double *tau1start, 
                   const double *tau2start, 
                   const double *component1start,
                   const double *e,
                   const double *f,
                   const double *g,
                   const double *rhostepdata,
                   const int* uselecuyer, 
                   const int* seedarray, 
                   const int* lecuyerstream,
                   const double *b0data, 
                   const double *B0data, 
                   double *logmarglikeholder, 
                   double *loglikeholder, 
                   const int *chib){

    MCMCPACK_PASSRNG2MODEL(MCMCnegbinReg_impl, 
                           betaout, nuout, rhoout, 
                           tau1out, tau2out, comp1out, comp2out,
                           sr1out, sr2out, mr1out, mr2out,
                           rhosizes,
                           Ydata, Yrow, Ycol, 
                           Xdata, Xrow, Xcol, 
                           burnin, mcmc, thin, verbose, 
                           betastart, nustart, rhostart,
                           tau1start, tau2start, component1start,
                           e, f, g, rhostepdata, 
                           b0data, B0data, 
                           logmarglikeholder, loglikeholder, chib);
  }//end MCMC
} // end extern "C"


#endif
