#ifndef MCMCNBUTIL_H
#define MCMCNBUTIL_H

#include "matrix.h"
#include "algorithm.h"
#include "distributions.h"
#include "la.h"
#include "stat.h"
#include "smath.h"
#include "rng.h"

using namespace std;
using namespace scythe;

struct hmm_state {
  Matrix<> s;
  Matrix<> ps;
  Matrix<> trans;
  Matrix<> nstate;
};

Matrix<> component_selector(const int yt);

template <typename RNGTYPE>
inline void init_aux(rng<RNGTYPE>& stream, 
                     const Matrix<>& Y,
                     Matrix<>& wr1, 
                     Matrix<>& mr1,
                     Matrix<>& sr1,
                     Matrix<>& wr2,
                     Matrix<>& mr2,
                     Matrix<>& sr2,
                     Matrix<>& nr2,
                     Matrix<>& component2) {
 
  Matrix<> rmat1 = component_selector(1);
  wr1 = rmat1(_, 0);
  mr1 = rmat1(_, 1);
  sr1 = rmat1(_, 2);
  const int n = Y.rows();

  for (int t = 0; t < n; t++) {
    int yt = (int) Y[t];
    if (yt > 0) {
      Matrix<> rmat2 = component_selector(yt);
      nr2[t] = rmat2.rows();
      for (unsigned int j = 0; j < rmat2.rows(); j++) {
        wr2(t, j) = rmat2(j, 0);
        mr2(t, j) = rmat2(j, 1);
        sr2(t, j) = rmat2(j, 2);
      }
      int comp = ceil((rmat2.rows()) * stream.runif());
      component2[t] = comp;
    }
  }

  return;
}

template <typename RNGTYPE>
inline int sample_discrete(rng<RNGTYPE>& stream,
                    const Matrix<>& probs) {
  Matrix <> cprobs(probs.rows(), 1);
  cprobs[0] = probs[0];
  for (unsigned int i = 1; i < probs.rows(); i++) {
    cprobs[i] = cprobs[i-1] + probs[i];
  }
  double U = stream.runif(); 
  int draw = 1;
  for (unsigned int h = 0; h < probs.rows(); h++) {
    if (cprobs[h] <= U && U < cprobs[h+1]) {
      draw = h + 2;
    }
  }
  return draw;
}


double rho_conditional_log_density(const double rho, 
                                   const Matrix<>& y,
                                   const Matrix<>& lambda,
                                   const double g,
                                   const double e,
                                   const double f);


template <typename RNGTYPE>
Matrix<> rho_slice_sampler(rng<RNGTYPE>& stream,
                     const Matrix<>& y,
                     const Matrix<>& lambda,
                     const double rho, 
                     const double step_out,
                     const double g, 
                     const double e, 
                     const double f) {

  double  U, L, R, V, end_val, this_slice, new_rho, new_slice, old_rho;
  const int M = 100;
  int J, K;

  old_rho = rho;
  this_slice = rho_conditional_log_density(old_rho, y, lambda, g, e, f);
 
  // Step (a): draw the slice
  U = stream.rexp(1);
  this_slice = this_slice - U;
  
  // Step (b): step out the from the old rho to create intervals
  U = stream.runif();
  L = rho - step_out * U;
  R = L + step_out;
  L = max(0.0, L);
  V = stream.runif();
  J = floor(M * V);
  K = (M - 1) - J;
      
  // Step out left
  end_val = rho_conditional_log_density(L, y, lambda, g, e, f);
  while (this_slice < end_val && J > 0) {
    L = max(0.0, L - step_out);
    end_val = rho_conditional_log_density(L, y, lambda, g, e, f);
    J--;
    R_CheckUserInterrupt();
  }

  // Step out right
  end_val = rho_conditional_log_density(R, y, lambda, g, e, f);
  while (this_slice < end_val && K > 0) {
    R += step_out;
    end_val = rho_conditional_log_density(R, y, lambda, g, e, f);
    K--;
    R_CheckUserInterrupt();
  }
  
  // Step (c): Draw from the slice with shrinkage
  U = stream.runif();
  new_rho = L + U * (R - L);
  new_slice = rho_conditional_log_density(new_rho, y, lambda, g, e, f);
  while (new_slice <= this_slice) {
    if (new_rho > rho) { 
      R = new_rho;
    } else {
      L = new_rho;
    }
    U = stream.runif();
    new_rho = L + U * (R - L);
    new_slice = rho_conditional_log_density(new_rho, y, lambda, g, e, f);
    R_CheckUserInterrupt();
  }

  Matrix<> rhostep(5, 1);
  rhostep[0] = new_rho;
  rhostep[1] = new_slice;
  rhostep[2] = abs(new_rho - rho);
  rhostep[3] = L;
  rhostep[4] = R;

  return rhostep;
  
}

//rho slice sampler -- prior only
template <typename RNGTYPE>
Matrix<> rho_prior_sampler(rng<RNGTYPE>& stream,
                           const double rho, 
                           const double step_out,
                           const double g, 
                           const double e, 
                           const double f) {

  double  U, L, R, V, end_val, this_slice, new_rho, new_slice, old_rho;
  const int M = 100;
  int J, K;

  old_rho = rho;
  this_slice = (e - 1) * log(old_rho) - (e+f) * log(old_rho + g);
 
  // Step (a): draw the slice
  U = stream.rexp(1);
  this_slice = this_slice - U;
  
  // Step (b): step out the from the old rho to create intervals
  U = stream.runif();
  L = rho - step_out * U;
  R = L + step_out;
  L = max(0.0, L);
  V = stream.runif();
  J = floor(M * V);
  K = (M - 1) - J;
      
  // Step out left
  end_val = (e - 1) * log(L) - (e+f) * log(L + g);
  while (this_slice < end_val && J > 0) {
    L = max(0.0, L - step_out);
    end_val = (e - 1) * log(L) - (e+f) * log(L + g);
    J--;
    R_CheckUserInterrupt();
  }

  // Step out right
  end_val = (e - 1) * log(R) - (e+f) * log(R + g);
  while (this_slice < end_val && K > 0) {
    R += step_out;
    end_val = (e - 1) * log(R) - (e+f) * log(R + g);
    K--;
    R_CheckUserInterrupt();
  }
  
  // Step (c): Draw from the slice with shrinkage
  U = stream.runif();
  new_rho = L + U * (R - L);
  new_slice = (e - 1) * log(new_rho) - (e+f) * log(new_rho + g);
  while (new_slice <= this_slice) {
    if (new_rho > rho) { 
      R = new_rho;
    } else {
      L = new_rho;
    }
    U = stream.runif();
    new_rho = L + U * (R - L);
    new_slice = (e - 1) * log(new_rho) - (e+f) * log(new_rho + g);
    R_CheckUserInterrupt();
  }

  Matrix<> rhostep(5, 1);
  rhostep[0] = new_rho;
  rhostep[1] = new_slice;
  rhostep[2] = abs(new_rho - rho);
  rhostep[3] = L;
  rhostep[4] = R;

  return rhostep;
  
}

//tau and component sampler
template <typename RNGTYPE>
Matrix<> tau_comp_sampler(rng<RNGTYPE>& stream, 
                          const int y,
                          const double mu,
                          const Matrix<>& wr1,
                          const Matrix<>& mr1,
                          const Matrix<>& sr1,
                          const Matrix<>& wr2,
                          const Matrix<>& mr2,
                          const Matrix<>& sr2,
                          const int nr2){

  const int nr1 = wr1.rows();
  Matrix<> post_tau1_mat(nr1, 1);
  double tau_t1;
  double tau_t2;
  int component1, component2;

  double xi = stream.rexp(mu);
  if (y == 0) {
    tau_t1 = 1 + xi;
    tau_t2 = 0.0;
    component2 = 0;
  } else {
    Matrix<> post_tau2_mat(nr2, 1);
    tau_t2 = stream.rbeta(y, 1);
    tau_t1 = 1-tau_t2+xi;
    for (int h=0; h < nr2; ++h){   
      double psi = dnorm(-log(tau_t2) - log(mu), mr2[h], sqrt(sr2[h]));
      post_tau2_mat[h] = wr2[h] * psi;
    }
    Matrix <> norm_post_tau2 = post_tau2_mat/sum(post_tau2_mat);
    component2 = sample_discrete(stream, norm_post_tau2);
  }
  for (int h = 0; h < nr1; ++h){   
    double psi = dnorm(-log(tau_t1) - log(mu), mr1[h], sqrt(sr1[h]));
    post_tau1_mat[h] = wr1[h] * psi;
  }
  Matrix <> norm_post_tau1 = post_tau1_mat/sum(post_tau1_mat);
  component1 = sample_discrete(stream, norm_post_tau1);

  Matrix<> TAUout(4, 1);
  TAUout[0] = tau_t1;
  TAUout[1] = tau_t2;
  TAUout[2] = component1;
  TAUout[3] = component2;

  return TAUout;
  
}

// derived from
template <typename RNGTYPE>
inline double sample_conparam(rng<RNGTYPE>& stream,
                       double curr_param,
                       const Matrix<>& numdata,
                       const int tot_class,
                       const double aa,
                       const double bb,
                       const int numiter) {

  const int n_group = numdata.rows();
  double gammaa, gammab;
  double alpha = curr_param;

  for (int i = 0; i < numiter; ++i) {
    double wvec_sum = 0;
    double svec_sum = 0;
    double u;
    
    for (int j = 0; j < n_group; ++j) {
      if (numdata[j] > 0) {
        wvec_sum += log(stream.rbeta(alpha + 1, numdata[j]));
        u = stream.runif();
        if (u < (numdata[j]/(alpha + numdata[j]))) {
          svec_sum++;
        }
      }
    }
    gammaa = aa + tot_class - svec_sum;
    gammab = bb - wvec_sum; 
    alpha = stream.rgamma(gammaa, gammab);
  }
  return alpha;
}

#endif
