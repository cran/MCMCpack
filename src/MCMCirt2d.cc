// MCMCirt2d.cc is C++ code to estimate a 2-dimensional item response
// theory model.  It is currently used with the MCMCirt2d.R front-end,
// but can also be used for higher dimensional models with other
// front-ends.

#include <iostream>
#include "Scythe_Matrix.h"
#include "Scythe_LA.h"
#include "Scythe_Simulate.h"
#include "Scythe_Stat.h"
#include "Scythe_IDE.h"
#include "Scythe_Math.h"

using namespace SCYTHE;
using namespace std;

// update latent data for standard item response models
Matrix<double>
irt_Z_update2 (Matrix<double> &Y, Matrix<double> &theta,
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
Matrix<double>
   irt_eta_update2 (Matrix<double> & Z, Matrix<double> & theta,
   const Matrix<double> & Mb0, const Matrix<double> & MB0) {

    // define constants
    int J = theta.rows();
    int D = 2;
    int K = Z.cols();
    Matrix<double> eta = Matrix<double>(K, D+1);
   
    // perform update 
    Matrix<double> theta_star = cbind(ones<double>(J,1),theta);
    for (int j=0; j<K; ++j){

       // form bill specific priors
       Matrix<double> b0 = Matrix<double>(D+1,1);
       b0(0,0) = Mb0(j,0);
       b0(1,0) = Mb0(j,1);
       b0(2,0) = Mb0(j,2);
           
       Matrix<double> toxpnd = Matrix<double>(6,1);
       toxpnd(0,0) = MB0(j,0);
       toxpnd(1,0) = MB0(j,1);
       toxpnd(2,0) = MB0(j,2);
       toxpnd(3,0) = MB0(j,3);
       toxpnd(4,0) = MB0(j,4);
       toxpnd(5,0) = MB0(j,5);
 
       Matrix<double> B0 = Matrix<double>(D+1);
       B0 = xpnd(toxpnd);
       
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
irt_theta_update2 (Matrix<double> & Z, Matrix<double> & eta, const 
   Matrix<double> & Mt0, const Matrix<double> & MT0) {
   
   int D = 2;
   int J = Z.rows();
   Matrix<double> theta = Matrix<double>(J,2);
   
   // perform update from multivariate Normal
   Matrix<double> alpha = eta(0, 0, eta.rows() - 1, 0);
   Matrix<double> beta = eta(0, 1, eta.rows() - 1, 1);
    
   for (int i=0; i<J; ++i) {
   
     // form justice specific priors
     Matrix<double> t0 = Matrix<double>(D,1);
     t0(0,0) = Mt0(i,0);
     t0(1,0) = Mt0(i,1);
   
     Matrix<double> toxpnd = Matrix<double>(3,1);
     toxpnd(0,0) = MT0(i,0);
     toxpnd(1,0) = MT0(i,1);
     toxpnd(2,0) = MT0(i,2);
     Matrix<double> T0 = Matrix<double>(D,D);
     T0 = xpnd(toxpnd);
     
     Matrix<double> Tmat = invpd(t(beta)*beta + T0);
     Matrix<double> tmat = Tmat*(T0 * t0 + t(beta)*(alpha -  
						t(Z(i,0,i,Z.cols() - 1))));
     Matrix<double> new_theta = t(tmat + cholesky(Tmat) * rnorm(D,1));  
     for(int l = 0; l<D; l++)
	     theta(i,l) = new_theta(0,l);
      }
   return theta;
}

extern "C"{

// function called by R to fit model
void
irt2dpost (double* sample, const int* samrow, const int* samcol,
           const double* votedata, const int* voterow, const int* votecol,
           const int* burnin, const int* gibbs,  const int* thin,
           const int* seed, const int* verbose, const double* tstart,
	       const int* tstartrow, const int* tstartcol, const int* thetafixed1,
	       const double* abstart,
	       const int* abstartrow, const int* abstartcol,
	       const double* t0data, const int* t0row, const int* t0col,
	       const double* T0data, const int* T0row, const int* T0col,	       
	       const double* b0data, const int* b0row, const int* b0col,
	       const double* B0data, const int* B0row, const int* B0col,
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
  Matrix<double> Malphabetastart(abstartcol[0],abstartrow[0], abstart);
  Malphabetastart = t(Malphabetastart);  
  Matrix<double> Mt0(t0col[0], t0row[0], t0data);
  Mt0 = t(Mt0);
  Matrix<double> MT0(T0col[0], T0row[0], T0data);
  MT0 = t(MT0);
  Matrix<double> Mb0(b0col[0], b0row[0], b0data);
  Mb0 = t(Mb0);
  Matrix<double> MB0(B0col[0], B0row[0], B0data);
  MB0 = t(MB0);
    
  // initialize seed (mersenne twister / use default seed unless specified)
  if(seed==0) set_mersenne_seed(5489UL);
  else set_mersenne_seed(seed[0]);
      
  // parameters
  int J = Mvote.cols();  // number of actors (justices)
  int K = Mvote.rows();  // number of bills (cases)
  int D = 2;             // dimensionality of the issue space
  int tot_iter = burnin[0] + gibbs[0];  
   
  // storage matrices (row major order)
  Matrix<double> theta_store = Matrix<double>(gibbs[0]/thin[0],J*D);
  Matrix<double> eta_store = Matrix<double>(gibbs[0]/thin[0],K*(D+1));
  Matrix<double> Z     = Matrix<double>(J,K);
  
  // starting values 
  Matrix<double> theta = Mthetastart;
  Matrix<double> eta   = Malphabetastart;
  int count    = 0;
 
  ///////////////////
  // Gibbs Sampler //
  ///////////////////
  
  for (int iter=0; iter < tot_iter; ++iter){
    
    // sample latent utilities (Z)
    Z = irt_Z_update2(Y, theta, eta);
        
    // sample bill parameters (eta)
    eta = irt_eta_update2(Z, theta, Mb0, MB0);
    
    // sample ideal points (theta)
    theta = irt_theta_update2(Z, eta, Mt0, MT0);
      
    // fix rotational invariance by forcing member
    // thetafixed1 to be negative on both dimensions
    if(theta(thetafixed1[0]-1,0) > 0) {
       for(int i=0; i<J; ++i)
           theta(i,0) = -1 * theta(i,0);
        }
    if(theta(thetafixed1[0]-1,1) > 0) {
       for(int i=0; i<J; ++i)
           theta(i,1) = -1 * theta(i,1);
        }

     // fix labeling invariance by forcing member
     // thetafixed2 to be lower on dim1 than dim2
     // if(theta(thetafixed2[0]-1,0) > theta(thetafixed2[0]-1,1)) {
     //   Matrix<double> holder = theta;
     //   for(int i=0; i<J; ++i) {
     //      theta(i,0) = holder(i,1);
     //      theta(i,1) = holder(i,0);
     //   } 
     // }
     

    // print results to screen
    if (verbose[0] == 1 && iter % 50 == 0){
      cout << " MCMCirt2d Iteration = " << iter << endl;
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

