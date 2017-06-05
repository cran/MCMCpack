//////////////////////////////////////////////////////////////////////////
// cHDPHSMMnegbin.cc is C++ code to estimate a HDP-HSMM with Negative
// Binomial emission distribution (Johnson & Willsky, 2013)
//
// Matthew Blackwell
// Department of Government
// Harvard University
// mblackwell@gov.harvard.edu

// Written 05/19/2017
//////////////////////////////////////////////////////////////////////////

#ifndef CHDPHSMMNEGBIN_CC
#define CHDPHSMMNEGBIN_CC

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



struct hsmm_state {
  Matrix<> s;
  Matrix<> s_norep;
  Matrix<> ps;
  Matrix<> durs;
  Matrix<> trans;
  Matrix<> nstate;
};

// define cumsumc
template <matrix_order RO, matrix_style RS, typename T,
          matrix_order PO, matrix_style PS> 
Matrix<T,RO,RS>
cumsumc (const Matrix<T,PO,PS>& A)
{
  Matrix<T,RO,RS> res (A.rows(), A.cols(), false);
  
  for (unsigned int j = 0; j < A.cols(); ++j) {
    res(0,j) = A(0,j);
    for (unsigned int i = 1; i < A.rows(); ++i) {
      res(i,j) = A(i,j) + res(i-1,j);
    }
  }
    
  
  return res;
}

template <typename T, matrix_order O, matrix_style S>
Matrix<T,O,Concrete>
cumsumc (const Matrix<T,O,S>& A)
{
  return cumsumc<O,Concrete>(A);
}         



template <typename RNGTYPE>
hsmm_state negbin_hdphsmm_reg_state_sampler(rng<RNGTYPE>& stream, 
                                            const int ns,
                                            const Matrix<>& Y,
                                            const Matrix<>& X,
                                            const Matrix<>& beta,
                                            const Matrix<>& P,
                                            const Matrix<>& nu,
                                            const Matrix<>& rho,
                                            const Matrix<>& omega,
                                            const double r){
  
  const int n = Y.rows();
  Matrix<int> trans(ns, ns);
  Matrix<int> nstate(ns, 1);
  Matrix<> B(n, ns);
  Matrix<> Bstar(n, ns);
  Matrix<> pD(n, ns);
  Matrix<> pDs(n, ns);
  Matrix<> pr1 = ones(ns, 1);
  Matrix<> py(n, ns);
  Matrix<> cy(n, ns);
  Matrix<> pstyt1(ns, 1);

  // no self-transitions
  Matrix<> Pbar = P - P % eye(ns);
  Pbar = log(Pbar) - log(Pbar * ones(ns, ns));
  for (int t = 0; t < n; t++) {
    int yt = (int) Y[t];
    Matrix<> lambda = exp(X(t,_) * ::t(beta)); 
    for (int j = 0; j< ns; ++j){
      pD(t, j) = log(dnbinom(t+1, r, omega[j]));
      pDs(t, j) = log(1 - pnbinom(t+1, r, omega[j]));
      py(t, j) += lngammafn(rho[j] + yt) - lngammafn(rho[j]) - lngammafn(yt +1);
      py(t, j) += rho[j] * log(rho[j]) + yt*log(lambda[j])- (rho[j] + yt) * log(rho[j] + lambda[j]);
    }
  }
  B(n-1,_) = 0.0;
  for (int t = (n-1); t >= 0; --t){
    Matrix<> unnorm_pstyt(ns, 1);
    Matrix<> pstyt(ns, 1);
    Matrix<> censored(1,ns);
    Matrix<> result((n-t),ns);
    Matrix<> Bhold(1,ns);
    Matrix<> cy = cumsumc(py(t, 0, n-1, ns-1));
    result = B(t,0,n-1,ns-1) + cy + pD(0, 0, (n-t-1), ns-1);
    Matrix<> maxes = maxc(result);
    for (int j = 0; j < ns; j++) {
      result(_,j) -= maxes[j];
    }
    Bstar(t,_) = log(sumc(exp(result))) + maxes;
    censored = pDs(n-t-1,_) + cy(cy.rows()-1, _);
    for (int j = 0; j < ns; j++) {
      maxes[j] = std::max(Bstar(t,j), censored[j]);
      Bstar(t,j) -= maxes[j];
      censored[j] -= maxes[j];
    }
    Bstar(t,_) = log(exp(Bstar(t,_)) + exp(censored)) + maxes;
    if (t > 0) {
      for (int j = 0; j < ns; j++) {
        Bhold = Pbar(j,_) + Bstar(t,_);        
        B(t-1,j) = log(sum(exp(Bhold - max(Bhold)))) + max(Bhold); 
      }
    }
  }

  int st;  
  Matrix<int> s(n, 1);
  Matrix<bool> ch(n,1,true,false);
  Matrix<int> durs(n,1);
  Matrix<> pstyn(ns, 1);
  Matrix<> Pst_1(ns,1);
  Matrix<> unnorm_pstyn(ns,1);
  Matrix<> cump(ns,1);
  int t = 0;
  int dur = 0;
  double durprob, this_pd, prior_pd;
  while (t < n) {
    ch[t] = true;
    if (t == 0) {
      unnorm_pstyn = Bstar(t,_);
    } else {
      st = s(t-1);
      unnorm_pstyn = Bstar(t,_) + Pbar(st-1,_);
    }
    pstyn = unnorm_pstyn - log(sum(exp(unnorm_pstyn - max(unnorm_pstyn)))) - max(unnorm_pstyn);
    pstyn = exp(pstyn);
    cump(0) = pstyn(0); 
    for (int j = 1; j < ns; j++) {
      cump(j) = cump(j-1) + pstyn(j);
    }
    double UU = stream.runif();
    s(t) = 1;
    for (int j = 1; j < ns; j++) {
      if (UU >= cump(j-1) && UU < cump(j)) {
        s(t) = j+1;
      }
    }
    if (t > 0) trans(st-1, s(t) - 1) += 1;    

    durprob = stream.runif();
    for (dur = 0; durprob > 0.0 && t+dur < n; dur++) {
      prior_pd = pD(dur, s(t) - 1);
      this_pd = prior_pd + sum(py(t, s(t) - 1, t+dur, s(t)-1));
      this_pd += B(t+dur, s(t) - 1) - Bstar(t, s(t) - 1);
      durprob -= exp(this_pd);
      s(t+dur) = s(t);            
    }
    durs[t] = dur;
    nstate[s(t) - 1] += dur;

    t += dur;        
  }
  Matrix<> s_norep = selif(s, durs > 0);
  Matrix<> durs_norep = selif(durs, durs > 0);
  hsmm_state result;
  result.s = s;
  result.s_norep = s_norep;
  result.durs = durs_norep;
  result.trans = trans;
  result.nstate = nstate;  
  return result;
} 




//////////////////////////////////////////// 
// HDPHSMMnegbinRegChange implementation.  
//////////////////////////////////////////// 
template <typename RNGTYPE>
void HDPHSMMnegbinReg_impl(rng<RNGTYPE>& stream, 
                           double *betaout, 
                           double *Pout, 
                           double *omegaout, 
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
                           double *gammaout,
                           double *alphaout,
                           double *rhosizes,
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
                           const double *nustart,
                           const double *rhostart,
                           const double *tau1start, 
                           const double *tau2start, 
                           const double *component1start,
                           const double *alphastart,
                           const double *gammastart,
                           const double *omegastart,
                           const double *a_alpha, 
                           const double *b_alpha,
                           const double *a_gamma, 
                           const double *b_gamma,
                           const double *a_omega, 
                           const double *b_omega,
                           const double *e,
                           const double *f,
                           const double *g,
                           const double *r,
                           const double *rhostepdata,
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
  Matrix <> omega(ns, 1, omegastart);
  double gamma = *gammastart;
  double alpha = *alphastart;
  Matrix <> gamma_prime(ns, 1, true, gamma/ns);
  
  Matrix<> beta_store(nstore, ns*k);
  Matrix<> P_store(nstore, ns*ns);
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
  Matrix<> nu_store(nstore, n);
  Matrix<> rho_store(nstore, ns);
  Matrix<> omega_store(nstore, ns);
  Matrix<> gamma_store(nstore, 1);
  Matrix<> alpha_store(nstore, 1);
  
  hsmm_state Sout;
  Matrix<> s(n,1);
  Matrix<> (n,ns);
  Matrix<int> trans_counts(ns,ns);
  Matrix<> V(ns,ns);
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
    if (debug > 0) Rprintf("\niter: %i\n-----\n", iter);
    if (debug > 0) Rprintf("before s \n");
    //////////////////////
    // 1. Sample s
    //////////////////////
    Sout = negbin_hdphsmm_reg_state_sampler(stream, ns, Y, X, beta, P, nu, rho, omega, *r);
    s = Sout.s;
    Matrix<> s_norep = Sout.s_norep;
    Matrix<> durs = Sout.durs;
    nstate = Sout.nstate;
    trans_counts = Sout.trans;

    if (debug > 0) Rprintf("before rho \n");
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
    
    if (debug > 0) Rprintf("before nu \n");
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

    if (debug > 0) Rprintf("before tau \n");
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


    if (debug > 0) Rprintf("before beta \n");
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




    if (debug > 0) Rprintf("before DA for P \n");
    // Augmented variables to make sampling P easier
    Matrix<int> froms = sumc(::t(Sout.trans));
    for (int j = 0; j < ns; j++) {
      if (debug > 1) Rprintf("j = %i, froms = %i P = %10.5f tc = %i\n", j, froms[j], P(j,j), trans_counts(j,j));
      int draw;
      for (int nj = 0; nj < froms[j]; nj++) {
        draw = std::floor(stream.rexp(-log(P(j,j))));
        // sometimes we get wonky draws here
        if (draw > 0) trans_counts(j,j) += draw;
      }
      if (debug > 1) Rprintf("DA tc[j,j] = %i\n", trans_counts(j,j));
    }

    //////////////////////
    // 3. Sample dish counts
    //////////////////////
    if (debug > 0) Rprintf("before dish counts \n");
    Matrix<int> rest_dishes(ns, ns);
    int n_jk;
    int m_jk;
    double m_num;
    for (int j = 0; j < ns; ++j) {
      for (int k = 0; k < ns; ++k) {
        n_jk = trans_counts(j,k);
        if (n_jk == 0) {
          rest_dishes(j,k) = 0;
        } else {
          m_jk = 0;
          for (int h = 0; h < n_jk; ++h) {
            m_num = alpha * gamma_prime(k);
            if (stream.runif() < (m_num/(m_num + h))) {
              m_jk++;
            }
          }
          rest_dishes(j,k) = m_jk;
        }
      }
    }

    //////////////////////
    // 3. Sample P
    //////////////////////
    if (debug > 0) Rprintf("before concentration params \n");
    // Gamma prime is beta in all HDP-H(S)MM notation
    Matrix<int> Nkdot = ::t(sumc(::t(trans_counts)));
    Matrix<int> Mkdot = ::t(sumc(::t(rest_dishes)));
    Matrix<int> Mdotk = ::t(sumc(rest_dishes));
    double alpha0 = alpha;
    double gamma0 = gamma;
    int Mbar_tot = sum(rest_dishes);

    if (debug > 1) {
      for (int j = 0; j < ns; j++) {
        Rprintf("Nkdot[%i] = %i  Mkdot[%i] = %i  Mdotk[%i] = %i  g_prime[%i] = %10.5f\n", j, Nkdot[j], j, Mkdot[j], j, Mdotk[j], j, gamma_prime[j]);
      }
    }
    // only sample these if the parameters make sense
    if (debug > 0) Rprintf("before alpha \n");
    if (*a_alpha > 0 && *b_alpha > 0) {
      alpha = sample_conparam(stream, alpha0, Nkdot, sum(Mkdot), *a_alpha, *b_alpha, 50);
    }
    if (debug > 0) Rprintf("alpha = %10.5f\nbefore gamma\n", alpha);
    if (*a_gamma > 0 && *b_gamma > 0) {
      gamma = sample_conparam(stream, gamma0, Mbar_tot, ns, *a_gamma, *b_gamma, 50);
    }
    if (debug > 0) Rprintf("gamma = %10.5f\nbefore gamma prime\n", gamma);
    Matrix<double> gamma_prime_dir(ns,1);
    for (int j = 0; j < ns; j++) {
      gamma_prime_dir[j] = gamma/ns + Mdotk[j];
    }
    if (min(gamma_prime_dir + 1) <= 1) gamma_prime_dir += std::numeric_limits<double>::epsilon();
    
    gamma_prime = stream.rdirich(gamma_prime_dir);
    while (std::isnan(gamma_prime[0])) {
      gamma_prime = stream.rdirich(gamma_prime_dir);
    }

    if (debug > 0) Rprintf("before P\n", alpha);
    for (int j = 0; j < ns; ++j) {
      Matrix<double> p_dirich_params(ns,1);
      for (int i = 0; i < ns; ++i) {
        p_dirich_params(i) = alpha * gamma_prime(i) + trans_counts(j,i);
      }
      if (min(p_dirich_params + 1) <= 1) p_dirich_params += std::numeric_limits<double>::epsilon();
      P(j,_) = ::t(stream.rdirich(p_dirich_params));
      
      // the rdirichlet function has numerical stability problems so sometimes we have to redraw
      while (std::isnan(P(j,0))) {
        if (debug > 0) Rprintf("\nRedrawing P due to numerical problems (iter %i)...\n", iter);
        P(j,_) = ::t(stream.rdirich(p_dirich_params));
      }
      // again, due to numerical issues we can get numerically close to 1, just use the mean instead
      if (P(j,j) == 1.0) {
        P(j,_) = p_dirich_params/sum(p_dirich_params);
      }
      
    }

    if (debug > 0) Rprintf("before durs \n");
    // //////////////////////
    // // 3. Sample duration parameters
    // //////////////////////
    int last_state = s(n-1) - 1;
    int cen_dur = durs[durs.rows()-1];

    for (int j = 0; j < ns; j++) {
      int times_visited = sum(s_norep == (j + 1));
      int uncen_total_dur = nstate[j];
    
      // augment the last censored duration with acceptance sampling
      if (j == last_state) {
        int uncen_dur = 0;
        double tail_prob;
        tail_prob = log(1.0 - pnbinom(cen_dur, *r, omega[j]));
        if (exp(tail_prob) > 0.1) {
           while (uncen_dur < cen_dur) {
            uncen_dur = stream.rnbinom(*r, omega[j]);
            R_CheckUserInterrupt();      
          }
        } else {
          if (debug > 0) Rprintf("direct sampling \n");
          double u = stream.runif();
          uncen_dur = cen_dur;
          while (u > 0) {
            u -= exp(log(dnbinom(uncen_dur, *r, omega[j])) - tail_prob);
            uncen_dur++;            
          }
        }
        uncen_total_dur += uncen_dur - cen_dur;
      }
      omega[j] = stream.rbeta(*a_omega + *r * times_visited, *b_omega + uncen_total_dur);
    }
    


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
      nu_store(count,_) = nu(_, 0);
      rho_store(count, _) = rho(_, 0);
      omega_store(count, _) = omega(_, 0);
      gamma_store(count) = gamma;
      alpha_store(count) = alpha;
      ++count; 
    }   

    if(*verbose > 0 && iter % *verbose == 0){
      Rprintf("\n\n HDPHSMMnegbinChange iteration %i of %i \n\n", (iter+1), tot_iter);
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
  
  
        
  R_CheckUserInterrupt();      
  
  for (int i = 0; i<(nstore*ns*k); ++i){
    betaout[i] = beta_store[i]; 
  }
  for (int i = 0; i<(nstore*ns*ns); ++i){
    Pout[i] = P_store[i]; 
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
    omegaout[i] = omega_store[i];
  }
  for (int i = 0; i<ns; ++i){
    rhosizes[i] = step_out[i];
  }
  for (int i = 0; i < nstore; ++i) {
    gammaout[i] = gamma_store(i);
    alphaout[i] = alpha_store(i);
  }
}
 
extern "C" {
  void cHDPHSMMnegbin(double *betaout, 
                      double *Pout,
                      double *omegaout,
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
                      double *gammaout,
                      double *alphaout,
                      double *rhosizes,
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
                      const double *nustart,
                      const double *rhostart,
                      const double *tau1start, 
                      const double *tau2start, 
                      const double *component1start,
                      const double *alphastart,
                      const double *gammastart,
                      const double *omegastart,
                      const double *a_alpha, 
                      const double *b_alpha,
                      const double *a_gamma, 
                      const double *b_gamma,
                      const double *a_omega, 
                      const double *b_omega,
                      const double *e,
                      const double *f,
                      const double *g,
                      const double *r,
                      const double *rhostepdata,
                      const int* uselecuyer, 
                      const int* seedarray, 
                      const int* lecuyerstream,
                      const double *b0data, 
                      const double *B0data) {

    MCMCPACK_PASSRNG2MODEL(HDPHSMMnegbinReg_impl, 
                           betaout, Pout, omegaout, sout, nuout, rhoout, 
                           tau1out, tau2out, comp1out, comp2out,
                           sr1out, sr2out, mr1out, mr2out,
                           gammaout, alphaout, rhosizes,
                           Ydata, Yrow, Ycol, 
                           Xdata, Xrow, Xcol, 
                           K, burnin, mcmc, thin, verbose, 
                           betastart, Pstart, nustart, rhostart,
                           tau1start, tau2start, component1start,
                           alphastart, gammastart, omegastart,
                           a_alpha, b_alpha, a_gamma, b_gamma,
                           a_omega, b_omega, e, f, g, r, rhostepdata, 
                           b0data, B0data);
  }//end MCMC
} // end extern "C"


#endif
