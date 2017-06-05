// cMCMCirtHier1d.cc is C++ code to estimate a one-dimensional item response
// theory model with subject-level predictors beta and common variance
// sigma2.
//
// ADM and KQ 1/15/2003
// ADM 7/28/2004 [updated to new Scythe version]
// completely rewritten and optimized for the 1-d case 8/2/2004 KQ
// storage changed to save memory KQ 1/27/2006
// DBP 7/3/07 [ported to scythe 1.0.x]
//
// MJM added second level and marginal likelihood thereof
// MJM implemented parameter expansion (alpha) for latent variance 2008-11-18

#ifndef CMCMCIRTHIER1D_CC
#define CMCMCIRTHIER1D_CC

#include "MCMCrng.h"
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

#include "MCMCfcds.h"

using namespace std;
using namespace scythe;

/* This fn basically defined in MCMCregress.cc
and should be moved to Scythe proper */
static double lndigamma(double theta, double a, double b) {
  double logf =  a * log(b) - lngammafn(a) + -(a+1) * log(theta) +
                 -b/theta;
  return logf;
  //pow(b, a) / gammafn(a) * pow(theta, -(a+1)) * exp(-b/theta);
}


// Parameter-Expanded Latent Data Update for
// 1d IRT. mjm, 2008-11-18
template <typename RNGTYPE>
double irt_W_update(Matrix<>& Z, const Matrix<>& X, const Matrix<>& theta,
		    const Matrix<>& eta, const double& alpha,
		    const double& px_a0, const double& px_b0,
		    const Matrix<>& etahat, const Matrix<>& thetahat,
		    rng<RNGTYPE>& stream) {

 // define constants
  const unsigned int J = theta.rows();
  const unsigned int K = eta.rows();

  double RSS=0.0;
  int df=0;
  // perform update from truncated Normal / standard Normals
  for (unsigned int i = 0; i < J; ++i) {
    for (unsigned int j = 0; j < K; ++j){
      const double Z_mean = alpha*( -eta(j,0) + theta(i) * eta(j,1) );
      const double Zhat = ( -etahat(j,0) + thetahat(i) * etahat(j,1) );
      if (X(i,j) == 1) {
        Z(i,j) = stream.rtbnorm_combo(Z_mean, alpha, 0);
	++df;
      } else if (X(i,j) == 0) {
        Z(i,j) = stream.rtanorm_combo(Z_mean, alpha, 0);
	++df;
      } else {
        Z(i,j) = stream.rnorm(Z_mean, std::pow(alpha, 2.0));
      }
      Z(i,j) /= alpha;
      const double e = Z(i,j) - Zhat;
      RSS += std::pow(e , 2.0);
    }
  }
  // Liu and Wu 1999, p1272:
  // draw a0 ~ IG(a0,b0). Then draw a1 ~ IG(nu0+a0RSS/2, a0+n/2)
  const double c_post = (px_a0 + df) * 0.5;
  const double d_post = (px_b0 + RSS) * 0.5;
  double alpha1 = stream.rigamma(c_post,d_post);

  //Rprintf("\nRSS: %5f alpha0: %5f alpha1: %5f\n",RSS,alpha,alpha1);
  return(std::sqrt(alpha1/alpha));
}
/*
template <typename RNGTYPE>
void hirt_level0_metrop (Matrix<>& theta, Matrix<>& thetahat,
			 const Matrix<>& Z,
			 const Matrix<>& eta,
			 const Matrix<>& beta, const Matrix<>& Xj,
			 const double& sigma2,
			 const double& alpha,
			 rng<RNGTYPE>& stream){
  //construct proposal sum(Zjk ~N( eta1k*(-eta0k + thetaj), alpha? ) )
  //proposal logdensity
  //current logdensity
  //exp( prop - cur)
  //runif & accept if > ratio
  //note the acceptance
  }
*/


/* MCMCirt1d implementation. */
template <typename RNGTYPE>
void MCMCirtHier1d_impl (rng<RNGTYPE>& stream,
			 const Matrix<int>& X, Matrix<>& theta,
			 Matrix<>& eta,
			 Matrix<>& thetahat, Matrix<>& etahat,
			 const Matrix<>& ab0,
			 const Matrix<>& AB0, const Matrix<>& Xj,
			 Matrix<>& beta, const Matrix<>& b0,
			 const Matrix<>& B0,
			 const double c0, const double d0,
			 unsigned int burnin,
			 unsigned int mcmc, unsigned int thin, unsigned int verbose,
			 bool storea, bool storei,
			 double* sampledata, unsigned int samplesize,
			 bool chib, double* logmarglike,
			 bool px, const double px_a0, const double px_b0,
			 bool metromix
			 )
{
  // constants
  const unsigned int J = X.rows();  // # subjects (justices, legislators)
  const unsigned int K = X.cols();  // number of items (cases, roll calls)
  const unsigned int L = Xj.cols(); // covariates Xj for mutheta
  const unsigned int tot_iter = burnin + mcmc;
  const unsigned int nsamp = mcmc / thin;

  // storage matrices (col major order)
  Matrix<> theta_store;
  Matrix<> eta_store;
  Matrix<> beta_store(nsamp,L);
  Matrix<> sigma2_store(nsamp,1);
  if (storea)
    theta_store = Matrix<>(nsamp, J);

  if (storei)
    eta_store = Matrix<>(nsamp, K*2);

  // starting values
  Matrix<> Z(J, K);

  // pre-compute
  const Matrix<> AB0ab0 = AB0 * ab0;
  const Matrix<> XpX = crossprod(Xj);
  Matrix<> XpY;
  double alpha;
  if(!px) {
     alpha = 1.0;
  } else {
     alpha = stream.rigamma(px_a0,px_b0);
  }
  double sigma2 = NormIGregress_sigma2_draw ( Xj, theta, beta, c0, d0, stream);

  unsigned int count = 0;
  // MCMC sampling occurs in this for loop
  for (unsigned int iter = 0; iter < tot_iter; ++iter){
    // sample latent data (Z) OR parameter-expanded data (W)
    if(!px) {
      irt_Z_update1(Z, X, theta,eta,stream);
    } else {
      // alpha here is the post variance of W
      alpha = irt_W_update(Z, X, theta,eta,alpha,px_a0,px_b0,
			   etahat,thetahat,stream);
    }

    // Following Tierney(1994) it would be nice to 'restart' the chain
    // using a metropolis jump on the entire level0, something like
    // once every thin*10 iterations or so.
    if( metromix && (iter % (thin*10) == 0)) {
      //      hirt_level0_metrop(theta, thetahat, Z, eta, etahat,
      //		beta, Xj, sigma2, alpha, stream);
    } else {
    // sample ability (ideal points) (theta)
    hirt_theta_update1(theta, thetahat, Z, eta, beta, Xj, sigma2,
    				  alpha, stream);

    // sample item (case, bill) parameters (eta)
    hirt_eta_update1(eta, etahat, Z, theta,  AB0, AB0ab0,
    		     alpha, stream);
    }

    XpY = t(Xj) * theta;
    beta = NormNormregress_beta_draw(XpX, XpY, b0, B0, sigma2, stream);

    // update level2 sigma2
    sigma2 = NormIGregress_sigma2_draw (Xj, theta, beta, c0, d0, stream);

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

      beta_store(count, _) = t(beta);
      sigma2_store(count,0) = sigma2;
      // store beta
      count++;
    }

    R_CheckUserInterrupt(); // allow user interrupts

  } // end Gibbs loop

  // return output
  Matrix<> output;
  if(! storei && storea) {        // only theta
    output = theta_store;
  } else if (storei && ! storea){ // only eta
    output = eta_store;
  } else {                        // everything
    output = cbind(theta_store, eta_store);
  }
  output = cbind(output, beta_store);  // always return beta,
  output = cbind(output, sigma2_store);// and sigma2.

  for (unsigned int i = 0; i < samplesize; ++i)
    sampledata[i] = output[i];

  // BEGIN MARGINAL LIKELIHOOD COMPUTATION
    if (chib == 1) {
     // marginal likelihood calculation stuff starts here
     const double sigma2star = meanc(sigma2_store)(0);
     const Matrix<> thetastar = t(meanc(theta_store));
     Matrix<> betastar = t(meanc(beta_store));
     double sigma2fcdsum = 0.0;
     XpY = t(Xj) * thetastar;
     // second set of Gibbs scans
     for (unsigned int iter = 0; iter < tot_iter; ++iter) {
       double sigma2 = NormIGregress_sigma2_draw (Xj, thetastar, beta, c0, d0,
                                                  stream);

       beta = NormNormregress_beta_draw (XpX, XpY, b0, B0, sigma2,
                                         stream);
       const Matrix<> e = gaxpy(Xj, (-1*beta), thetastar);
       const Matrix<> SSE = crossprod (e);
       const double c_post = (c0 + X.rows ()) * 0.5;
       const double d_post = (d0 + SSE(0)) * 0.5;
       sigma2fcdsum += lndigamma(sigma2star, c_post, d_post);

       // print output to stdout
       if(verbose > 0 && iter % verbose == 0) {
         Rprintf("\n\nMCMCregress (reduced) iteration %i of %i \n",
                 (iter+1), tot_iter);
       }
       R_CheckUserInterrupt(); // allow user interrupts
     } // end MCMC loop

     double sigma2fcdmean = sigma2fcdsum / static_cast<double>(tot_iter);
     const double sig2_inv = 1.0 / sigma2star;
     const Matrix<> sig_beta = invpd (B0 + XpX * sig2_inv);
     Matrix<> betahat = sig_beta * gaxpy(B0, b0, XpY*sig2_inv);
     double logbetafcd = 0.0;
     for (unsigned int i=0; i<L; ++i) {
       logbetafcd += lndnorm(betastar[i], betahat[i], sig_beta[i]);
     }
     Rprintf("logPost: beta %10.5f sigma2 %10.5f\n",logbetafcd,sigma2fcdmean);
     // calculate loglikelihood at (betastar, sigma2star)
     double sigmastar = sqrt(sigma2star);
     Matrix<> xi = Xj * betastar;
     double loglike = 0.0;
     // is there a reason to do this and not lndmvn?
     for (unsigned int i = 0; i < L; ++i) {
      loglike += lndnorm(thetastar(i), xi(i), sigmastar);
     }
     Rprintf("logLike: %10.5f\n",loglike);

     // calculate log prior ordinate
     const double logpriorsig = (lndigamma(sigma2star, c0/2.0, d0/2.0));
     double logpriorbeta = lndmvn(betastar, b0, invpd(B0));
     if (L==1) { logpriorbeta = lndnorm(betastar[0], b0[0], 1.0/B0(0)); }
     const double logprior = logpriorsig+logpriorbeta;
     Rprintf("logPrior: %10.5f\n",logprior);

     // put pieces together and print the marginal likelihood
     logmarglike[0] = loglike + logprior - logbetafcd - (sigma2fcdmean);
     Rprintf("logM: %10.5f\n",logmarglike[0]);

   }
}

extern "C" {

  void
  cMCMCirtHier1d(double* sampledata, const int* samplerow, const int* samplecol,
		const int* Xdata, const int* Xrow, const int* Xcol,
		const int* burnin, const int* mcmc,  const int* thin,
		const int *uselecuyer, const int *seedarray,
		const int *lecuyerstream,
		const int* verbose, const double* thetastartdata,
		const int* thetastartrow, const int* thetastartcol,
		const double* astartdata, const int* astartrow, const int* astartcol,
		const double* bstartdata, const int* bstartrow, const int* bstartcol,
		const double* ab0data, const int* ab0row, const int* ab0col,
		const double* AB0data, const int* AB0row, const int* AB0col,
		const double* Xjdata, const int* Xjrow, const int* Xjcol,
		const double* betastartdata, const int* betastartrow, const int* betastartcol,
		const double* b0data, const int* b0row, const int* b0col,
		const double* B0data, const int* B0row, const int* B0col,
		const double* c0, const double* d0,
		const int* storei, const int* storea,
		double* logmarglikeholder, const int* chib, const int* px,
		const double* px_a0, const double* px_b0)
  {
    // put together matrices
    const Matrix<int> X(*Xrow, *Xcol, Xdata);
    Matrix<> theta(*thetastartrow, *thetastartcol, thetastartdata);
    Matrix<>thetahat(*thetastartrow, *thetastartcol, thetastartdata);
    Matrix<> a(*astartrow, *astartcol, astartdata);
    Matrix<> b(*bstartrow, *bstartcol, bstartdata);
    const Matrix<> ab0(*ab0row, *ab0col, ab0data);
    const Matrix<> AB0(*AB0row, *AB0col, AB0data);
    Matrix<> eta  = cbind(a,b);
    Matrix<> etahat = cbind(a,b);
    const Matrix<> Xj(*Xjrow, *Xjcol, Xjdata);
    Matrix<> beta(*betastartrow, *betastartcol, betastartdata);
    const Matrix<> b0(*b0row,*b0col,b0data);
    const Matrix<> B0(*B0row,*B0col,B0data);
    const int samplesize = (*samplerow) * (*samplecol);
    const bool metromix = 0;

    MCMCPACK_PASSRNG2MODEL(MCMCirtHier1d_impl, X, theta,
			   eta,
			   thetahat, etahat, ab0, AB0,
			   Xj, beta, b0, B0, *c0, *d0,
			   *burnin, *mcmc, *thin,
			   *verbose, *storea, *storei,
			   sampledata, samplesize,
			   *chib, logmarglikeholder, *px, *px_a0, *px_b0,
			   metromix);
  }
}

#endif
