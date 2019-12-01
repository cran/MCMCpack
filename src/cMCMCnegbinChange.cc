//////////////////////////////////////////////////////////////////////////
// cMCMCnegbinChange.cc is C++ code to estimate a Chib-style changepoint 
// with Negative Binomial emission distribution
//
// Matthew Blackwell
// Department of Government
// Harvard University
// mblackwell@gov.harvard.edu

// Written 05/19/2017
//////////////////////////////////////////////////////////////////////////

#ifndef CMCMCNEGBINCHANGE_CC
#define CMCMCNEGBINCHANGE_CC

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



template <typename RNGTYPE>
hmm_state negbin_reg_state_sampler(rng<RNGTYPE>& stream, 
                                   const int m, 
                                   const Matrix<>& Y,
                                   const Matrix<>& X,
                                   const Matrix<>& beta,
                                   const Matrix<>& P,
                                   const Matrix<>& rho,
                                   double &loglike,
                                   const int& fixed_m){
  
  const int ns = m + 1;
  const int n = Y.rows();

  Matrix<> trans(ns, ns);
  Matrix<> nstate(ns, 1);
  Matrix<> F(n, ns);
  Matrix<> pr1(ns, 1);
  if (fixed_m == 1) {
    pr1[0] = 1;
  } else {
    pr1 = 1;
  }
  Matrix<> py(ns, 1);
  Matrix<> pstyt1(ns, 1);
  Matrix<> logP = log(P);

  for (int t=0; t<n ; ++t){
    int yt = (int) Y[t];
    Matrix<> lambda = exp(X(t,_)*::t(beta)); 
    for (int j = 0; j< ns; ++j){
      py[j] = 0;
      py[j] += lngammafn(rho[j] + yt) - lngammafn(rho[j]) - lngammafn(yt +1);
      py[j] += rho[j] * log(rho[j]) + yt*log(lambda[j])- (rho[j] + yt) * log(rho[j] + lambda[j]);
    }
    if (t==0) {
      pstyt1 = log(pr1);
    } else {
      for (int j = 0; j < ns; j++) {
        Matrix<> Phold = ::t(F(t-1,_)) + logP(_,j);
        if (R_finite(max(Phold))) {
          pstyt1[j] =  log(sum(exp(Phold - max(Phold)))) + max(Phold);
        } else {
          pstyt1[j] = log(0.0);
        }
      }
    }
    Matrix<> unnorm_pstyt = pstyt1 + py;
    double log_norm = log(sum(exp(unnorm_pstyt - max(unnorm_pstyt)))) + max(unnorm_pstyt);
    Matrix<> pstyt = unnorm_pstyt - log_norm;
    for (int j=0; j<ns ; ++j) F(t,j) = pstyt(j);
    loglike += log_norm;
  }
  
  Matrix<int> s(n, 1);                        
  Matrix<> ps = Matrix<>(n, ns);  
  ps(n-1,_) = exp(F(n-1,_));

  Matrix<> unnorm_pstyn(ns, 1);
  Matrix<> pstyn(ns, 1);

  if (fixed_m == 1) {
    s(n-1) = ns;
    nstate[ns - 1] = 1;
    trans(ns - 1, ns - 1) = 1;
  } else {
    unnorm_pstyn = ::t(exp(F(n-1,_)));
    pstyn = unnorm_pstyn/sum(unnorm_pstyn);
    s[n-1] = sample_discrete(stream, pstyn);
    nstate[s[n-1]-1] += 1;
  }

  double pone = 0.0;
  int t = n-2;
  while (t >= 0){
    int st = s(t+1);
    Matrix<> Pst_1 = ::t(logP(_,st-1)); 
    unnorm_pstyn = F(t,_)+Pst_1;
    pstyn = unnorm_pstyn - log(sum(exp(unnorm_pstyn - max(unnorm_pstyn)))) - max(unnorm_pstyn);
    pstyn = exp(pstyn);
    if (st==1) {
      s(t) = 1;
      trans(0, 0) += 1;
      nstate[0] += 1;
    } else {
      pone = pstyn(st-2);
      if(stream.runif() < pone) {
        s(t) = st-1;
        trans(st - 2, st - 1) += 1;        
      }  else {
        s(t) = st;
        trans(st - 1, st - 1) += 1;
      }
      nstate[s(t) - 1] += 1;
    }
    ps(t,_) = pstyn;
    --t;
  }

  hmm_state result;
  result.s = s;
  result.ps = ps;
  result.trans = trans;
  result.nstate = nstate;
  return result;
} 



//////////////////////////////////////////// 
// MCMCnegbinRegChange implementation.  
//////////////////////////////////////////// 
template <typename RNGTYPE>
void MCMCnegbinRegChangepoint_impl(rng<RNGTYPE>& stream, 
                                   double *betaout, 
                                   double *Pout, 
                                   double *psout, 
                                   double *sout,   
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
                                   const int *m, 
                                   const int *burnin, 
                                   const int *mcmc, 
                                   const int *thin, 
                                   const int *verbose, 
                                   const double *betastart, 
                                   const double *Pstart, 
                                   const double *nustart,
                                   const double *rhostart,
                                   const double *tau1start, 
                                   const double *tau2start, 
                                   const double *component1start,
                                   const double *a,
                                   const double *b,
                                   const double *e,
                                   const double *f,
                                   const double *g,
                                   const double *rhostepdata,
                                   const double *b0data, 
                                   const double *B0data, 
                                   const double *A0data,
                                   const int *fixed_m,
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
  const int ns = *m + 1;        
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

  const Matrix <> A0(ns, ns, A0data);
 
  Matrix <> nu(n, 1, nustart);
  Matrix <> rho(ns, 1, rhostart);
  Matrix <> rho_slice(ns, 1, true, 0.0);
  Matrix <> step_out(ns, 1, rhostepdata);
  Matrix <> beta(ns, k, betastart);
  Matrix <> tau1(n, 1, tau1start);
  Matrix <> tau2(n, 1, tau2start);
  Matrix <> component1(n, 1, component1start);
  Matrix <> component2(n, 1);
  Matrix <> P(ns, ns, Pstart);

  Matrix<> beta_st(ns, k);
  Matrix<> P_st(ns, ns);
  Matrix<> rho_st(ns, 1);
  Matrix<> beta_store(nstore, ns*k);
  Matrix<> P_store(nstore, ns*ns);
  Matrix<> ps_store(n, ns);
  Matrix<> s_store(nstore, n);
  Matrix<> nstate(ns, 1);   
  Matrix<int> component1_store(nstore, n);
  Matrix<int> component2_store(nstore, n);
  Matrix<> tau1_store(nstore, n);
  Matrix<> tau2_store(nstore, n);
  Matrix<> sr1_store(nstore, n);
  Matrix<> sr2_store(nstore, n);
  Matrix<> mr1_store(nstore, n);
  Matrix<> mr2_store(nstore, n);
  Matrix<> nu_store(nstore, n);
  Matrix<> rho_store(nstore, ns);
  
  hmm_state Sout;
  Matrix<> s(n,1);
  Matrix<> ps(n,ns);
  Matrix<> V(ns,ns);
  Matrix<> rhostep;
  double curr_lp_max = 0.0;
  double logpost = 0.0;
    
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
      Rprintf("Sampling s...\n");
    }

    //////////////////////
    // 1. Sample s
    //////////////////////
    if (ns > 1) {
      double loglike = 0.0;
      Sout = negbin_reg_state_sampler(stream, *m, Y, X, beta, P, rho, loglike, *fixed_m);
      s = Sout.s;
      ps = Sout.ps;
      nstate = Sout.nstate;
      logpost += loglike;
      if ((iter == *burnin) || ((iter > *burnin) && (logpost > curr_lp_max))) {
        beta_st = beta;
        P_st = P;
        rho_st = rho;
        curr_lp_max = logpost;
      }
      
    } else {
      s = 1;
      ps = 1;
    }
    if (debug > 0) {
      Rprintf("Sampling rho...\n");
    }

    //////////////////////
    // 5. Sample rho
    //////////////////////  
    // We need to evaluate the density in the first iteration for the slice sampler
    for (int j = 0; j < ns; j++) {
      if (nstate[j] > 0) {
        Matrix<> yj = selif(Y, s == (j + 1));
        Matrix<> xj = selif(X, s == (j + 1));
        Matrix<> lam = exp(xj * ::t(beta(j,_)));
        rhostep = rho_slice_sampler(stream, yj, lam, rho[j], step_out[j], *g, *e, *f);
      } else {
        rhostep = rho_prior_sampler(stream, rho[j], step_out[j], *g, *e, *f);
      }
      rho[j] = rhostep[0];
      if (iter > 10) {
        step_out[j] = (1/((double) iter)) * rhostep[2] + ((((double) iter)-1)/((double) iter)) * step_out[j];
      }
    }
    

    if (debug > 0) {
      Rprintf("Sampling nu...\n");
    }

    //////////////////////
    // 6. Sample nu
    //////////////////////
    Matrix<> lambda(n, 1);
    for (int t = 0; t < n; t++) {
      int st = (int) s[t];
      int yt = (int) Y[t];
      Matrix<> xt = X(t, _);
      Matrix<> mu_t = exp(xt* ::t(beta(st-1,_)));
      lambda[t] = mu_t[0];
      nu[t] = stream.rgamma(rho[st-1] + yt, rho[st-1] + lambda[t]);
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

    Matrix<int> pnstate(ns, 1);
    
    for (int j = 0; j <ns ; ++j){
      for (int i = 0; i<n; ++i){
	if (s[i] == (j+1)) { 
          int yt = (int) Y[i];
          if (yt > 0) {
            pnstate[j] = pnstate[j] + 1;
          }
	}
      }
      if (nstate[j] == 0) {
        if (k == 1) {
          beta(j, _) = stream.rnorm(b0[0], sqrt(1/B0[0]));
        } else {
          beta(j,_) = stream.rmvnorm(b0, invpd(B0));
        }
      } else {
        int tot_rows = nstate[j] + pnstate[j];
        Matrix<> yjp(tot_rows, 1);
        Matrix<> Xjp(tot_rows, k);
        Matrix<> wip(tot_rows, 1);

        yjp(0, 0, nstate[j] - 1, 0) = selif(y_tilde, s == (j + 1));
        Xjp(0, 0, nstate[j] - 1, k-1) = selif(X, s == (j + 1));
        wip(0, 0, nstate[j] - 1, 0) = selif(Sigma_inv_sum, s == (j + 1));
        if (pnstate[j] > 0) {
          yjp(nstate[j], 0, tot_rows - 1, 0) = selif(yp_tilde, (s == (j+1)) & (Y > 0));
          Xjp(nstate[j], 0, tot_rows - 1, k - 1) = selif(X, (s == (j + 1)) & (Y > 0));
          wip(nstate[j], 0, tot_rows - 1, 0) = selif(Sigma_plus_inv_sum, (s == (j + 1)) & (Y > 0));
        }
        Matrix<> Xwj(Xjp.rows(), k);
        for (unsigned int h = 0; h<Xjp.rows(); ++h){
          Xwj(h, _) = Xjp(h,_)*wip[h];
        }

        Matrix<> Bn = invpd(B0 + ::t(Xwj)*Xwj); 
        Matrix<> bn = Bn*gaxpy(B0, b0, ::t(Xwj)*yjp);

        if (k == 1) {
          beta(j, _) = stream.rnorm(bn[0], sqrt(Bn[0]));
        } else {
          beta(j,_) = stream.rmvnorm(bn, Bn);
        }
      }
    }

    if (debug > 0) {
      Rprintf("Sampling P...\n");
    }

    //////////////////////
    // 3. Sample P
    //////////////////////        
    double shape1 = 0;
    double shape2 = 0;    
    P(ns-1, ns-1) = 1; 
    
    for (int j =0; j< (ns-1); ++j){
      shape1 =  A0(j,j) + nstate[j] - 1;  
      shape2 =  A0(j,j+1) + 1;
      P(j,j) = stream.rbeta(shape1, shape2);
      P(j,j+1) = 1 - P(j,j);
    }

    if (iter >= *burnin && ((iter % *thin)==0)){
      Matrix<> tbeta = ::t(beta); 
      for (int i=0; i<(ns*k); ++i){
	beta_store(count,i) = tbeta[i];
      }
      for (int j=0; j<ns*ns; ++j){    
	P_store(count,j)= P[j];
      }
      for (int l=0; l<n ; ++l){           
	ps_store(l,_) = ps_store(l,_) + ps(l,_);           
      }
      s_store(count,_) = s(_, 0);
      tau1_store(count,_) = tau1;
      tau2_store(count,_) = tau2;
      component1_store(count,_) = component1;
      component2_store(count,_) = component2;
      sr1_store(count,_) = sr1_hold;
      sr2_store(count,_) = sr2_hold;
      mr1_store(count,_) = mr1_hold;
      mr2_store(count,_) = mr2_hold;
      nu_store(count,_) = nu(_, 0);
      rho_store(count, _) = rho(_, 0);
      ++count; 
    }   

    // calculate prior
    logpost = 0.0;
    for (int j=0; j<ns ; ++j){
      if (k == 1) {
        logpost += lndnorm(beta[j], b0[0], B0inv[0]); 
      }
      else {
        logpost += lndmvn(::t(beta(j,_)), b0, B0inv);
      }
      logpost += (*e - 1) * log(rho[j]) - (*e + *f) * log(rho[j] + *g);
    }   
    if (ns > 1) {
      for (int j =0; j< (ns-1); ++j){
        logpost += log(dbeta(P(j,j), A0(j,j), A0(j,j+1))); 
      }
    }
    
    if(*verbose > 0 && iter % *verbose == 0){
      Rprintf("\n\n MCMCnegbinChange iteration %i of %i \n", (iter+1), tot_iter);
      for (int j = 0;j<ns; ++j){
	Rprintf("The number of observations in state %i is %10.5f\n", j+1, static_cast<double>(nstate[j]));     
      }
      for (int i = 0; i<ns; ++i){
        if (nstate[i] > 0) {
          Rprintf("rho in state %i is %10.5f\n", i+1, rho[i]); 
          for (int j = 0; j<k; ++j){
            Rprintf("beta(%i) in state %i is %10.5f\n", j+1, i+1, beta(i, j));
          }
        }
      }
    }
    
  }// end MCMC loop 
  
  
  
  //////////////////////////////////////
  if (*chib == 1) {
  //////////////////////////////////////
    if(*verbose > 0) {
      Rprintf("\nCalculating marginal likelihood...\n");
    }

    //////////////////////
    // beta
    //////////////////////      
    Matrix<> density_beta(nstore, ns);  
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
      
      Matrix<int> pnstate(ns, 1);
      for (int j = 0; j <ns ; ++j){
        nstate[j] = 0;
        for (int i = 0; i<n; ++i){
          if (s_store(iter,i) == (j+1)) { 
            int yt = (int) Y[i];
            nstate[j] += 1;
            if (yt > 0) {
              pnstate[j] = pnstate[j] + 1;
            }
          }
        }
        int tot_rows = nstate[j] + pnstate[j];
        Matrix<> yjp(tot_rows, 1);
        Matrix<> Xjp(tot_rows, k);
        Matrix<> wip(tot_rows, 1);

        s = ::t(s_store(iter,_));
        yjp(0, 0, nstate[j] - 1, 0) = selif(y_tilde, s == (j + 1));
        Xjp(0, 0, nstate[j] - 1, k-1) = selif(X, s == (j + 1));
        wip(0, 0, nstate[j] - 1, 0) = selif(Sigma_inv_sum, s == (j + 1));
        if (pnstate[j] > 0) {
          yjp(nstate[j], 0, tot_rows - 1, 0) = selif(yp_tilde, (s == (j+1)) & (Y > 0));
          Xjp(nstate[j], 0, tot_rows - 1, k - 1) = selif(X, (s == (j + 1)) & (Y > 0));
          wip(nstate[j], 0, tot_rows - 1, 0) = selif(Sigma_plus_inv_sum, (s == (j + 1)) & (Y > 0));
        }
        Matrix<> Xwj(Xjp.rows(), k);
        for (unsigned int h = 0; h<Xjp.rows(); ++h) {
          Xwj(h, _) = Xjp(h,_)*wip[h];
        }

        Matrix<> Bn = invpd(B0 + ::t(Xwj)*Xwj); 
        Matrix<> bn = Bn*gaxpy(B0, b0, ::t(Xwj)*yjp);
        
        if (k == 1) {
          density_beta(iter, j) = exp(lndnorm(beta_st[j], bn[0], sqrt(Bn[0])));
        } else {
          density_beta(iter, j) = exp(lndmvn(::t(beta_st(j,_)), bn, Bn));
        }
      }
      
    }

    double pdf_beta = sum(log(meanc(density_beta)));
    //////////////////////
    // P
    //////////////////////
    Matrix<> density_P(nstore, ns);
    double pdf_P;
    double ll = 0.0;

    if (ns > 1) {
    for (int iter = 0; iter < nstore; ++iter){
      Sout = negbin_reg_state_sampler(stream, *m, Y, X, beta_st, P, rho, ll, *fixed_m);
      s = Sout.s;
      ps = Sout.ps; 
      nstate = Sout.nstate;

      for (int j = 0; j < ns; j++) {
        if (nstate[j] > 0) {
          Matrix<> yj = selif(Y, s == (j + 1));
          Matrix<> xj = selif(X, s == (j + 1));
          Matrix<> lam = exp(xj * ::t(beta(j,_)));
          rhostep = rho_slice_sampler(stream, yj, lam, rho[j], step_out[j], *g, *e, *f);
        } else {
          rhostep = rho_prior_sampler(stream, rho[j], step_out[j], *g, *e, *f);
        }
        rho[j] = rhostep[0];
      }

      double shape1 = 0;
      double shape2 = 0;    
      P(ns-1, ns-1) = 1; 
      
      for (int j =0; j< (ns-1); ++j){
        shape1 =  A0(j,j) + nstate[j] - 1;  
        shape2 =  A0(j,j+1) + 1;         
        P(j,j) = stream.rbeta(shape1, shape2);
        P(j,j+1) = 1 - P(j,j);
        density_P(iter, j) = dbeta(P_st(j,j), shape1, shape2); 
      } 
      density_P(iter, ns-1) = 1; 
      
    }           
    pdf_P = log(prod(meanc(density_P)));
    } else {
      pdf_P = 0;
    }

    //////////////////////
    // rho/nu
    //////////////////////
    Matrix<> density_rho(nstore, ns);  
    for (int iter = 0; iter < nstore; iter++) {
      Sout = negbin_reg_state_sampler(stream, *m, Y, X, beta_st, P_st, rho, ll, *fixed_m);
      s = Sout.s;
      ps = Sout.ps; 
      nstate = Sout.nstate;

      for (int j = 0; j < ns; j++) {
        if (nstate[j] > 0) {
          Matrix<> yj = selif(Y, s == (j + 1));
          Matrix<> xj = selif(X, s == (j + 1));
          Matrix<> lam = exp(xj * ::t(beta(j,_)));
          rhostep = rho_slice_sampler(stream, yj, lam, rho[j], step_out[j], *g, *e, *f);
        } else {
          rhostep = rho_prior_sampler(stream, rho[j], step_out[j], *g, *e, *f);
        }
        rho[j] = rhostep[0];
        double L = rhostep[3];
        double R = rhostep[4];
        if (rho_st[j] > L && rho_st[j] < R) {
          if (((R-L) + 1) <= 1) {
            density_rho(iter, j) = 1/std::numeric_limits<double>::epsilon();
          } else {
            density_rho(iter, j) = 1/(R-L);
          }
        } else {
          density_rho(iter, j) = 0.0;
        }
      }
      
    }
    
    double pdf_rho = sum(log(meanc(density_rho)));
    
  
    //////////////////////
    // likelihood
    //////////////////////       
    Matrix<> F = Matrix<>(n, ns);
    Matrix<> like(n, 1);
    Matrix<> pr1 = Matrix<>(ns, 1);
    pr1[0] = 1;
    Matrix<> py(ns, 1);
    Matrix<> pstyt1(ns, 1);
    
    for (int t=0; t<n ; ++t){
      int yt = (int) Y[t];
      Matrix<> lambda = exp(X(t,_)*::t(beta_st)); 
      for (int j = 0; j< ns; ++j){
        py[j] = 0;
        py[j] += lngammafn(rho_st[j] + yt) - lngammafn(rho_st[j]) - lngammafn(yt +1);
        py[j] += rho_st[j] * log(rho_st[j]) + yt*log(lambda[j])- (rho_st[j] + yt) * log(rho_st[j] + lambda[j]);
        py[j] = exp(py[j]);
      }
      if (t==0) pstyt1 = pr1;
      else {
        pstyt1 =  ::t(F(t-1,_)*P_st); 
      }
      Matrix<> unnorm_pstyt = pstyt1%py;
      Matrix<> pstyt = unnorm_pstyt/sum(unnorm_pstyt); 
      for (int j=0; j<ns ; ++j){
        F(t,j) = pstyt(j);
      }
      like[t] = sum(unnorm_pstyt);
    }
    
    const double loglike = sum(log(like));

    //////////////////////
    // prior
    //////////////////////
    Matrix<> density_beta_prior(ns, 1);
    Matrix<> density_P_prior(ns, 1);
    Matrix<> density_rho_prior(ns, 1);
    density_P[ns-1] = 1; //
  
    for (int j=0; j<ns ; ++j){
      if (k == 1) {
        density_beta_prior[j] = lndnorm(beta_st[j], b0[0], B0inv[0]); 
      }
      else {
        density_beta_prior[j] = lndmvn(::t(beta_st(j,_)), b0, B0inv);
      }
      density_rho_prior[j] = (*e - 1) * log(rho_st[j]) - (*e + *f) * log(rho_st[j] + *g);
    }   
    if (ns > 1) {
      for (int j =0; j< (ns-1); ++j){
        density_P_prior[j] = log(dbeta(P_st(j,j), A0(j,j), A0(j,j+1))); 
      }
    }
    
    // compute marginal likelihood
    const double logprior = sum(density_beta_prior) + sum(density_P_prior) + sum(density_rho_prior); 
    const double logmarglike = (loglike + logprior) - (pdf_beta + pdf_P + pdf_rho);
    if (*verbose >0 ){
      Rprintf("\nlogmarglike = %10.5f\n", logmarglike);
      Rprintf("loglike = %10.5f\n", loglike);
      Rprintf("log_prior = %10.5f\n", logprior);
      Rprintf("log_beta = %10.5f\n", pdf_beta);
      Rprintf("log_P = %10.5f\n", pdf_P);
      Rprintf("log_rho = %10.5f\n", pdf_rho);
    }
    logmarglikeholder[0] = logmarglike;
    loglikeholder[0] = loglike;
    
  }
        
  R_CheckUserInterrupt();      
  
  for (int i = 0; i<(nstore*ns*k); ++i){
    betaout[i] = beta_store[i]; 
  }
  for (int i = 0; i<(nstore*ns*ns); ++i){
    Pout[i] = P_store[i]; 
  }
  for (int i = 0; i<(n*ns); ++i){
    psout[i] = ps_store[i]; 
  }
  for (int i = 0; i<(nstore*n); ++i){
    sout[i] = s_store[i];
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
  for (int i = 0; i<(nstore*ns); ++i){
    rhoout[i] = rho_store[i];
  }
  for (int i = 0; i<ns; ++i){
    rhosizes[i] = step_out[i];
  }
}
 
extern "C" {
  void cMCMCnegbinChange(double *betaout, 
                         double *Pout, 
                         double *psout, 
                         double *sout,   
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
                         const int *m, 
                         const int *burnin, 
                         const int *mcmc, 
                         const int *thin, 
                         const int *verbose, 
                         const double *betastart, 
                         const double *Pstart, 
                         const double *nustart,
                         const double *rhostart,
                         const double *tau1start, 
                         const double *tau2start, 
                         const double *component1start,
                         const double *a,
                         const double *b,
                         const double *e,
                         const double *f,
                         const double *g,
                         const double *rhostepdata,
                         const int* uselecuyer, 
                         const int* seedarray, 
                         const int* lecuyerstream,
                         const double *b0data, 
                         const double *B0data, 
                         const double *A0data,
                         const int *fixed_m,
                         double *logmarglikeholder, 
                         double *loglikeholder, 
                         const int *chib){

    MCMCPACK_PASSRNG2MODEL(MCMCnegbinRegChangepoint_impl, 
                           betaout, Pout, psout, sout, nuout, rhoout, 
                           tau1out, tau2out, comp1out, comp2out,
                           sr1out, sr2out, mr1out, mr2out,
                           rhosizes,
                           Ydata, Yrow, Ycol, 
                           Xdata, Xrow, Xcol, 
                           m, burnin, mcmc, thin, verbose, 
                           betastart, Pstart, nustart, rhostart,
                           tau1start, tau2start, component1start,
                           a, b, e, f, g, rhostepdata, 
                           b0data, B0data, A0data, fixed_m,
                           logmarglikeholder, loglikeholder, chib);
  }//end MCMC
} // end extern "C"


#endif
