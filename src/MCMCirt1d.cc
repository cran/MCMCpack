// MCMCirt1d.cc is C++ code to estimate a one-dimensional item response
// theory model. 

#include <iostream>
#include "Scythe_Matrix.h"
#include "Scythe_Simulate.h"
#include "Scythe_Stat.h"
#include "Scythe_LA.h"
#include "Scythe_IDE.h"
#include "Scythe_Math.h"

using namespace SCYTHE;
using namespace std;

// update latent data for standard item response models
Matrix<double>
irt_Z_update (Matrix<double> &Y, Matrix<double> &theta,
                            Matrix<double> &eta) {
   // define constants
   int J = theta.rows();
   int K = eta.rows();
   Matrix<double> Z = Matrix<double>(J,K);
    
   // perform update from truncated Normal / standard Normals
   Matrix<double> theta_star = cbind(ones<double>(J,1), theta);
   for (int i=0; i<J; ++i) {
      for (int j=0; j<K; ++j){
         Matrix<double> Z_mean = theta_star(i, 0, i,
                 theta_star.cols() - 1) *
                 t(eta(j, 0, j, eta.cols() - 1));
         if (Y(i,j) == 1.0)
            Z(i,j) = rtbnorm_combo(Z_mean(0,0), 1.0, 0);
         if (Y(i,j) == 0.0)
            Z(i,j) = rtanorm_combo(Z_mean(0,0), 1.0, 0);
         if (Y(i,j) != 1.0 && Y(i,j) != 0.0)
            Z(i,j) = rnorm(Z_mean(0,0), 1.0);
       }
   }
   return Z;
}


// update bill parameters for item response model
// this function takes two K x 1 vectors of prior means (Mamean and Mbmean)
// two K x 1 vectors of prior variances (Mavar and Mbvar) and
// a K x 1 vector of prior covariances (Mabcov).
// note: works only for one-dimensional case
Matrix<double>
   irt_eta_update (Matrix<double> & Z, Matrix<double> & theta,
   const Matrix<double> & Mamean, const Matrix<double> & Mbmean,
   const Matrix<double> & Mavar, const Matrix<double> & Mbvar,
   const Matrix<double> & Mabcov ) {

    // define constants
    int J = theta.rows();
    int D = 1;
    int K = Z.cols();
    Matrix<double> eta = Matrix<double>(K, D+1);
   
    // perform update 
    Matrix<double> theta_star = cbind(ones<double>(J,1),theta);
    for (int j=0; j<K; ++j){

       // form bill specific priors
       Matrix<double> b0 = Matrix<double>(D+1,1);
       Matrix<double> B0 = eye<double>(D+1);
       b0(0,0) = Mamean(j,0);
       b0(1,0) = Mbmean(j,0);
       B0(0,0) = Mavar(j,0);
       B0(1,1) = Mbvar(j,0);
       B0(0,1) = Mabcov(j,0);
       B0(1,0) = Mabcov(j,0);
       
       Matrix<double> Emat = invpd(t(theta_star) * theta_star + invpd(B0));    
       Matrix<double> emat = Emat * (t(theta_star) * 
          Z(0, j, Z.rows() - 1, j) + invpd(B0) * b0);
       Matrix<double> new_eta = emat + cholesky(Emat) * rnorm(D+1,1);
       for(int l = 0; l<D+1; l++)
          eta(j,l) = new_eta(l,0);
       }
   return eta; 
}


// update ideal points for one dimensional item response model 
// this function takes a J x 1 vector of prior means (t0)
// and a J x 1 vector of prior variances (T0)
// note: works only for one-dimensional case
Matrix<double>
irt_theta_update (Matrix<double> & Z, Matrix<double> & eta, const 
   Matrix<double> & t0, const Matrix<double> & T0) {
   
   int J = Z.rows();
   Matrix<double> theta = Matrix<double>(J,1);
   
   // perform update from multivariate Normal
   Matrix<double> alpha = eta(0, 0, eta.rows() - 1, 0);
   Matrix<double> beta = eta(0, 1, eta.rows() - 1, 1);

   Matrix<double> Tmat = pow(t(beta) * beta + T0, -1);   
   for (int i=0; i<J; ++i) {
      Matrix<double> tmat = Tmat[i]*(T0[i] * t0[i] + t(beta)*(alpha -  
						t(Z(i,0,i,Z.cols() - 1))));
      Matrix<double> new_theta = t(tmat + sqrt(Tmat[i]) * rnorm(1,1));
      for(int l = 0; l<1; l++)
	     theta(i,l) = new_theta(0,l);
      }
   return theta;
}

extern "C"{

// function called by R to fit model
void
irt1dpost (double* sample, const int* samrow, const int* samcol,
           const double* votedata, const int* voterow, const int* votecol,
           const int* burnin, const int* gibbs,  const int* thin,
           const int* seed, const int* verbose, const double* tstart,
	       const int* tstartrow, const int* tstartcol, const int* thetafixed,
	       const double* bstart, const int* bstartrow, const int* bstartcol,
	       const double* astart, const int* astartrow, const int* astartcol,
	       const double* t0data, const int* t0row, const int* t0col,
	       const double* T0data, const int* T0row, const int* T0col,	       
	       const double* ameandata, const int* ameanrow, const int* ameancol,
	       const double* bmeandata, const int* bmeanrow, const int* bmeancol,
	       const double* avardata, const int* avarrow, const int* avarcol,
	       const double* bvardata, const int* bvarrow, const int* bvarcol,
	       const double* abcovdata, const int* abcovrow, const int* abcovcol,
	       const int* store
	       ) {

  // put together matrices
  Matrix<double> Msample(samcol[0], samrow[0], sample);
  Msample = t(Msample);
  Matrix<double> Mvote(votecol[0], voterow[0], votedata);
  Mvote = t(Mvote);  
  Matrix<double> Y = t(Mvote);
  Matrix<double> Mthetastart(tstartcol[0],tstartrow[0], tstart);
  Mthetastart = t(Mthetastart);
  Matrix<double> Malphastart(astartcol[0],astartrow[0], astart);
  Malphastart = t(Malphastart);  
  Matrix<double> Mbetastart(bstartcol[0],bstartrow[0], bstart);
  Mbetastart = t(Mbetastart);
  Matrix<double> Mt0(t0col[0], t0row[0], t0data);
  Mt0 = t(Mt0);
  Matrix<double> MT0(T0col[0], T0row[0], T0data);
  MT0 = t(MT0);
  Matrix<double> Mamean(ameancol[0], ameanrow[0], ameandata);
  Mamean = t(Mamean);
  Matrix<double> Mbmean(bmeancol[0], bmeanrow[0], bmeandata);
  Mbmean = t(Mbmean);
  Matrix<double> Mavar(avarcol[0], avarrow[0], avardata);
  Mavar = t(Mavar);  
  Matrix<double> Mbvar(bvarcol[0], bvarrow[0], bvardata);
  Mbvar = t(Mbvar);
  Matrix<double> Mabcov(abcovcol[0], abcovrow[0], abcovdata);
  Mabcov = t(Mabcov);
    
  // initialize seed (mersenne twister / use default seed unless specified)
  if(seed==0) set_mersenne_seed(5489UL);
  else set_mersenne_seed(seed[0]);
      
  // parameters
  int J = Mvote.cols();  // number of actors (justices)
  int K = Mvote.rows();  // number of bills (cases)
  int D = 1;             // dimensionality of the issue space
  int tot_iter = burnin[0] + gibbs[0];  
   
  // storage matrices (row major order)
  Matrix<double> theta_store = Matrix<double>(gibbs[0]/thin[0],J*D);
  Matrix<double> eta_store = Matrix<double>(gibbs[0]/thin[0],K*(D+1));
  Matrix<double> Z     = Matrix<double>(J,K);
  
  // starting values 
  Matrix<double> theta = Mthetastart;
  Matrix<double> eta   = cbind(Malphastart, Mbetastart);
  int count    = 0;
 
  ///////////////////
  // Gibbs Sampler //
  ///////////////////
  
  for (int iter=0; iter < tot_iter; ++iter){
    
    // sample latent utilities (Z)
    Z = irt_Z_update(Y, theta, eta);
        
    // sample bill parameters (eta)
    eta = irt_eta_update(Z, theta, Mamean, Mbmean, Mavar, Mbvar, Mabcov);
    
    // sample ideal points (theta)
    theta = irt_theta_update(Z, eta, Mt0, MT0);
      
    // fix rotational invariance by forcing a member to
    // have a negative ideal point
    if(theta(thetafixed[0]-1,0) > 0) {
        for(int i=0; i<J; ++i)
            theta(i,0) = -1 * theta(i,0);
        }
       
    // print results to screen
    if (verbose[0] == 1 && iter % 50 == 0){
      cout << " MCMCirt1d Iteration = " << iter << endl;
      cout << " theta = " << t(theta).toString() << endl;     
    }

    // store results
    if ((iter >= burnin[0]) && ((iter % thin[0]==0))) {
    
       // store ideal points
       Matrix<double> theta_store_vec = reshape(theta,1,J*D);
       for (int l=0; l<J*D; ++l)
          theta_store(count, l) = theta_store_vec(0,l);
        
       // store bill parameters
       Matrix<double> eta_store_vec = reshape(eta,1,K*(D+1));
       for (int l=0; l<K*(D+1); ++l)
              eta_store(count, l) = eta_store_vec(0,l);
       count++;
        }
     } // end Gibbs loop

    // return output
    Matrix<double> output;
    if(store[0] == 0) {
       output = theta_store;
    }
    else {
       output = cbind(theta_store, eta_store);
    }
    
    int loop = samrow[0] * samcol[0];
    for (int i=0; i<loop; ++i)
      sample[i] = output[i];

  }
}


