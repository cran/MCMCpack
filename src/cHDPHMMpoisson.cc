//////////////////////////////////////////////////////////////////////////
// cHDPHMMpoisson.cc is C++ code to estimate a (sticky) HDP-HMM with
// Poisson emission distribution
// (Fox, Sudderth, Jordan, & Willsky, 2013)
//
// Matthew Blackwell
// Department of Government
// Harvard University
// mblackwell@gov.harvard.edu

// Written 05/19/2017
//////////////////////////////////////////////////////////////////////////

#ifndef CHDPHMMPOISSON_CC
#define CHDPHMMPOISSON_CC

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
hmm_state poisson_hdp_reg_state_sampler(rng<RNGTYPE>& stream,
                                          const int& ns,
                                          const Matrix<>& Y,
                                          const Matrix<>& X,
                                          const Matrix<>& beta,
                                          const Matrix<>& P){

  const int n = Y.rows();
  Matrix<> trans(ns, ns);
  Matrix<> nstate(ns, 1);
  Matrix<> M(ns, n);
  Matrix<> pr1 = ones(ns, 1);
  Matrix<> py(ns, n);
  Matrix<> pstyt1(ns, 1);
  Matrix<> logP = log(P);
  for (int t = (n-1); t >= 0; --t){
    int yt = (int) Y[t];
    Matrix<> lambda = exp(X(t,_) * ::t(beta));
    for (int j = 0; j< ns; ++j){
      py(j, t) = log(dpois(yt, lambda[j]));
    }
    Matrix<> unnorm_pstyt(ns, 1);
    Matrix<> pstyt(ns, 1);

    if (t == (n - 1)) {
      pstyt1 = py(_,t);
    } else {
      pstyt1 = M(_,t+1) + py(_,t);
    }
    for (int j = 0; j < ns; j++) {
      Matrix<> Phold = ::t(logP(j,_)) + pstyt1;
      unnorm_pstyt[j] = log(sum(exp(Phold-max(Phold)))) + max(Phold);
    }
    M(_,t) = unnorm_pstyt;
  }

  int st = 1;
  Matrix<int> s(n, 1);
  Matrix<> ps(n, ns);
  Matrix<> pstyn(ns, 1);
  Matrix<> Pst_1(ns,1);
  Matrix<> unnorm_pstyn(ns,1);
  Matrix<> cump(ns,1);
  Matrix<> durs(n, 1);
  Matrix<> s_norep(n, 1);
  for (int t = 0; t < n; ++t) {
    if (t == 0) {
      unnorm_pstyn = M(_,t+1) + py(_,t);

    } else {
      st = s(t-1);
      Pst_1 = ::t(logP(st-1,_));

      if (t == (n - 1)) {
        unnorm_pstyn = py(_,t) + Pst_1;
      } else {
        unnorm_pstyn = M(_,t+1) + py(_,t) + Pst_1;
      }
    }
    pstyn = unnorm_pstyn - log(sum(exp(unnorm_pstyn - max(unnorm_pstyn)))) - max(unnorm_pstyn);
    pstyn = exp(pstyn);
    s[t] = sample_discrete(stream, pstyn);
    if (t > 0) {
      trans(st-1, s(t)- 1) += 1;
    }
    nstate[s(t)-1] += 1;
    ps(t,_) = pstyn;
  }

  hmm_state result;
  result.s = s;
  result.ps = ps;
  result.trans = trans;
  result.nstate = nstate;

  return result;
}





////////////////////////////////////////////
// HDPHMMpoissonReg implementation.
////////////////////////////////////////////
template <typename RNGTYPE>
void HDPHMMpoissonReg_impl(rng<RNGTYPE>& stream,
                           double *betaout,
                           double *Pout,
                           double *psout,
                           double *sout,
                           double *tau1out,
                           double *tau2out,
                           int *comp1out,
                           int *comp2out,
                           double *sr1out,
                           double *sr2out,
                           double *mr1out,
                           double *mr2out,
                           double *gammaout,
                           double *akout,
                           double *thetaout,
                           const double *Ydata,
                           const int *Yrow,
                           const int *Ycol,
                           const double *Xdata,
                           const int *Xrow,
                           const int *Xcol,
                           const int *K,
                           const int *burnin,
                           const int *mcmc,
                           const int *thin,
                           const int *verbose,
                           const double *betastart,
                           const double *Pstart,
                           const double *tau1start,
                           const double *tau2start,
                           const double *component1start,
                           const double *gammastart,
                           const double *akstart,
                           const double *thetastart,
                           const double *a_alpha,
                           const double *b_alpha,
                           const double *a_gamma,
                           const double *b_gamma,
                           const double *a_theta,
                           const double *b_theta,
                           const double *b0data,
                           const double *B0data)
{
  const Matrix <> Y(*Yrow, *Ycol, Ydata);
  const Matrix <> X(*Xrow, *Xcol, Xdata);
  const int tot_iter = *burnin + *mcmc;
  const int nstore = *mcmc / *thin;
  const int n = Y.rows();
  const int k = X.cols();
  const int ns = *K;
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

  Matrix <> nu(n, 1, true, 1.0);
  Matrix <> beta(ns, k, betastart);
  Matrix <> tau1(n, 1, tau1start);
  Matrix <> tau2(n, 1, tau2start);
  Matrix <> component1(n, 1, component1start);
  Matrix <> component2(n, 1);
  Matrix <> P(ns, ns, Pstart);
  Matrix <> gamma_prime(ns, 1);
  double gamma = *gammastart;
  double alpha_p_kappa = *akstart;
  double theta = *thetastart;
  double alpha = (1 - theta) * alpha_p_kappa;
  double kappa = theta * alpha_p_kappa;

  Matrix<> beta_store(nstore, ns*k);
  Matrix<> P_store(nstore, ns*ns);
  Matrix<> ps_store(n, ns);
  Matrix<> s_store(nstore, n);
  Matrix<int> nstate(ns, 1);
  Matrix<int> component1_store(nstore, n);
  Matrix<int> component2_store(nstore, n);
  Matrix<> tau1_store(nstore, n);
  Matrix<> tau2_store(nstore, n);
  Matrix<> sr1_store(nstore, n);
  Matrix<> sr2_store(nstore, n);
  Matrix<> mr1_store(nstore, n);
  Matrix<> mr2_store(nstore, n);
  Matrix<> gamma_store(nstore, 1);
  Matrix<> theta_store(nstore, 1);
  Matrix<> ak_store(nstore, 1);

  hmm_state Sout;
  Matrix<> s(n,1);
  Matrix<> ps(n,ns);
  Matrix<> V(ns,ns);

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

    if (debug > 0) Rprintf("iter: %i\n-----\n", iter);
    if (debug > 0) Rprintf("Sampling s...\n");
    //////////////////////
    // 1. Sample s
    //////////////////////
    if (ns > 1) {
      Sout = poisson_hdp_reg_state_sampler(stream, ns, Y, X, beta, P);
      s = Sout.s;
      nstate = Sout.nstate;
    } else {
      s = 1;
      ps = 1;
    }


    if (debug > 0) Rprintf("Sampling tau...\n");
    //////////////////////
    // 4. Sample tau
    //////////////////////
    Matrix<> TAUout;
    for (int t = 0; t < n; t++) {
      int yt = (int) Y[t];
      int st = (int) s[t];
      Matrix<> xt = X(t, _);
      Matrix<> mu_t = exp(xt* ::t(beta(st-1,_)));
      double mut = mu_t[0];

      TAUout = tau_comp_sampler(stream, yt, mut, wr1, mr1, sr1,
                                wr2(t,_), mr2(t,_), sr2(t,_), nr2[t]);
      tau1[t] = TAUout[0];
      tau2[t] = TAUout[1];
      component1[t] = TAUout[2];
      component2[t] = TAUout[3];
    }

    if (debug > 0) Rprintf("Sampling beta...\n");
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

    int beta_count = 0;
    int pcount = 0;
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
      beta_count = beta_count + nstate[j];
      pcount = pcount + pnstate[j];
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

    if (debug > 0) Rprintf("Sampling dish counts...\n");
    //////////////////////
    // 3. Sample dish counts
    //////////////////////
    Matrix<> rest_dishes(ns, ns);
    Matrix<> rest_dishes_over(ns, ns);
    Matrix<> sum_w(ns, 1);
    int n_jk;
    int m_jk;
    double m_num;
    double pp;
    double gam_thet;
    for (int j = 0; j < ns; ++j) {
      for (int k = 0; k < ns; ++k) {
        n_jk = Sout.trans(j,k);
        if (n_jk == 0) {
          rest_dishes(j,k) = 0;
        } else {
          m_jk = 1;
          for (int h = 1; h < n_jk; ++h) {
            m_num = alpha * gamma_prime(k);
            if (j == k) m_num += kappa;
            m_jk += stream.runif() < (m_num/(m_num + h));
          }
          rest_dishes(j,k) = m_jk;
          rest_dishes_over(j,k) = m_jk;
        }
      }
      if (rest_dishes(j,j) > 0) {
        gam_thet = gamma_prime(j) * (1 - theta);
        pp =  theta/(gam_thet + theta);
        sum_w(j) = stream.rbinom(rest_dishes(j,j), pp);
        rest_dishes_over(j,j) = rest_dishes(j,j) - sum_w(j);
      }
    }

    if (debug > 0) Rprintf("Sampling gamma_prime...\n");
    //////////////////////
    // 3. Sample P
    //////////////////////
    Matrix<> gamma_prime_dir = gamma/ns + ::t(sumc(rest_dishes_over));
    gamma_prime = stream.rdirich(gamma_prime_dir);
    // sometimes with single regime stuff, you get nans from this draw.
    if (std::isnan(gamma_prime(1))) {
      gamma_prime = gamma/ns;
    }
    if (debug > 0) Rprintf("Sampling P...\n");
    for (int j = 0; j < ns; ++j) {
      Matrix<> p_dirich_params(ns,1);
      for (int i = 0; i < ns; ++i) {
        p_dirich_params(i) = alpha * gamma_prime(i) + Sout.trans(j,i);
        if (i == j) p_dirich_params(i) += kappa;
      }
      if (min(p_dirich_params + 1) <= 1) p_dirich_params += std::numeric_limits<double>::epsilon();

      P(j,_) = ::t(stream.rdirich(p_dirich_params));
    }


    if (debug > 0) Rprintf("Sampling concentration params...\n");
    // //////////////////////
    // // 3. Sample hyperparams
    // //////////////////////
    Matrix<> Nkdot = sumc(::t(Sout.trans));
    Matrix<> Mkdot = sumc(::t(rest_dishes));
    Matrix<> Nkdot_valid = selif(::t(Nkdot), ::t(Nkdot) > 0);
    Matrix<> Mkdot_valid = selif(::t(Mkdot), ::t(Mkdot) > 0);
    double ak0 = alpha_p_kappa;
    double gamma0 = gamma;
    Matrix<> Mdotk = (sumc(rest_dishes_over) > 0);
    int Kbar = sum(Mdotk);
    int Mbar_tot = sum(rest_dishes_over);
    // only sample these if the parameters make sense
    if (debug > 0) Rprintf("Sampling alpha+kappa...\n");
    if (*a_alpha > 0 && *b_alpha > 0) {
      alpha_p_kappa = sample_conparam(stream, ak0, Nkdot_valid, sum(Mkdot_valid), *a_alpha, *b_alpha, 50);
    }
    if (debug > 0) Rprintf("Sampling gamma...\n");
    if (*a_gamma > 0 && *b_gamma > 0) {
      gamma = sample_conparam(stream, gamma0, Mbar_tot, Kbar, *a_gamma, *b_gamma, 50);
    }
    if (debug > 0) Rprintf("Sampling theta...\n");
    if (*a_theta > 0 && *b_theta > 0) {
      theta = stream.rbeta(*a_theta + sum(sum_w), *b_theta + (sum(rest_dishes) - sum(sum_w)));
    }
    kappa = theta * alpha_p_kappa;
    alpha = (1 - theta) * alpha_p_kappa;


    if (iter >= *burnin && ((iter % *thin)==0)){
      Matrix<> tbeta = ::t(beta);
      for (int i=0; i<(ns*k); ++i){
	beta_store(count,i) = tbeta[i];
      }
      for (int j=0; j<ns*ns; ++j){
	P_store(count,j)= P[j];
      }
      s_store(count,_) = s(_, 0);
      tau1_store(count,_) = tau1(_, 0);
      tau2_store(count,_) = tau2(_, 0);
      component1_store(count,_) = component1(_, 0);
      component2_store(count,_) = TAUout(_, 3);
      sr1_store(count,_) = sr1_hold;
      sr2_store(count,_) = sr2_hold;
      mr1_store(count,_) = mr1_hold;
      mr2_store(count,_) = mr2_hold;
      theta_store(count) = theta;
      gamma_store(count) = gamma;
      ak_store(count) = alpha_p_kappa;
      ++count;
    }

    if(*verbose > 0 && iter % *verbose == 0){
      Rprintf("\n\n HDPHMMpoisson iteration %i of %i \n\n", (iter+1), tot_iter);
      for (int j = 0;j<ns; ++j){
	Rprintf("The number of observations in state %i is %10.5f\n", j+1, static_cast<double>(nstate[j]));
      }
      for (int i = 0; i<ns; ++i){
        if (nstate[i] > 0) {
          for (int j = 0; j<k; ++j){
            Rprintf("beta(%i) in state %i is %10.5f\n", j+1, i+1, beta(i, j));
          }
        }
      }
    }

  }// end MCMC loop

  R_CheckUserInterrupt();

  for (int i = 0; i<(nstore*ns*k); ++i){
    betaout[i] = beta_store[i];
  }
  for (int i = 0; i<(nstore*ns*ns); ++i){
    Pout[i] = P_store[i];
  }
  for (int i = 0; i<(nstore*n); ++i){
    sout[i] = s_store[i];
    tau1out[i] = tau1_store[i];
    tau2out[i] = tau2_store[i];
    comp1out[i] = component1_store[i];
    comp2out[i] = component2_store[i];
    sr1out[i] = sr1_store[i];
    sr2out[i] = sr2_store[i];
    mr1out[i] = mr1_store[i];
    mr2out[i] = mr2_store[i];

  }
  for (int i = 0; i < nstore; ++i) {
    gammaout[i] = gamma_store(i);
    thetaout[i] = theta_store(i);
    akout[i] = ak_store(i);
  }
}

extern "C" {
  void cHDPHMMpoisson(double *betaout,
                      double *Pout,
                      double *psout,
                      double *sout,
                      double *tau1out,
                      double *tau2out,
                      int *comp1out,
                      int *comp2out,
                      double *sr1out,
                      double *sr2out,
                      double *mr1out,
                      double *mr2out,
                      double *gammaout,
                      double *akout,
                      double *thetaout,
                      const double *Ydata,
                      const int *Yrow,
                      const int *Ycol,
                      const double *Xdata,
                      const int *Xrow,
                      const int *Xcol,
                      const int *K,
                      const int *burnin,
                      const int *mcmc,
                      const int *thin,
                      const int *verbose,
                      const double *betastart,
                      const double *Pstart,
                      const double *tau1start,
                      const double *tau2start,
                      const double *component1start,
                      const double *gammastart,
                      const double *akstart,
                      const double *thetastart,
                      const double *a_alpha,
                      const double *b_alpha,
                      const double *a_gamma,
                      const double *b_gamma,
                      const double *a_theta,
                      const double *b_theta,
                      const int* uselecuyer,
                      const int* seedarray,
                      const int* lecuyerstream,
                      const double *b0data,
                      const double *B0data) {

    MCMCPACK_PASSRNG2MODEL(HDPHMMpoissonReg_impl,
                           betaout, Pout, psout, sout,
                           tau1out, tau2out, comp1out, comp2out,
                           sr1out, sr2out, mr1out, mr2out,
                           gammaout, akout, thetaout,
                           Ydata, Yrow, Ycol,
                           Xdata, Xrow, Xcol,
                           K, burnin, mcmc, thin, verbose,
                           betastart, Pstart,
                           tau1start, tau2start, component1start,
                           gammastart, akstart, thetastart,
                           a_alpha, b_alpha, a_gamma, b_gamma,
                           a_theta, b_theta,
                           b0data, B0data);
  }//end MCMC
} // end extern "C"


#endif
