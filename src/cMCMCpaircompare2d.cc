//////////////////////////////////////////////////////////////////////////
// MCMCpaircompare.cc is C++ code to estimate a pairwise comparison model. 
//
// KQ 3/18/2015
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////

#ifndef MCMCPAIRCOMPARE2D_CC
#define MCMCPAIRCOMPARE2D_CC

#include<vector>

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

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;

/* MCMCpaircompare2d implementation. */

// update latent Ystar values in 2-d paired comparisons model
template <typename RNGTYPE>
void paircompare2d_Ystar_update (Matrix<>& Ystar,
				 const Matrix<unsigned int>& MD, 
				 const Matrix<>& theta, const Matrix<>& gamma,
				 rng<RNGTYPE>& stream){

  const unsigned int N = MD.rows();
  for (unsigned int i = 0; i < N; ++i){
    const double gamma_i = gamma(MD(i,0));
    const double mean = std::cos(gamma_i)*theta(MD(i, 1), 0) + std::sin(gamma_i)*theta(MD(i, 1), 1) -
      std::cos(gamma_i)*theta(MD(i, 2), 0) - std::sin(gamma_i)*theta(MD(i, 2), 1);

    const bool above = MD(i,1) == MD(i,3); // cand 1 chosen
    const bool below = MD(i,2) == MD(i,3); // cand 2 chosen
    if (above){
      Ystar(i) = stream.rtbnorm_combo(mean, 1.0, 0.0);
    }
    else if (below){
      Ystar(i) = stream.rtanorm_combo(mean, 1.0, 0.0);
    }
    else{
      Ystar(i) = stream.rnorm(mean, 1.0);
    }
  }
}





// update theta values in 2-d paired comparisons model
template <typename RNGTYPE>
void paircompare2d_theta_update (Matrix<>& theta, const Matrix<>& Ystar, 
				 const Matrix<unsigned int>& MD,
				 const Matrix<>& gamma,
				 const Matrix<unsigned int>& theta_n,
				 const Matrix<>& theta_eq,
				 const Matrix<>& theta_ineq,
				 const vector< vector < double* > >& theta_Ystar_ptr,
				 const vector< vector < double* > >& theta_gamma_ptr,
				 const vector< vector < vector < double* > > >& theta_theta_ptr,
				 const vector< vector < double > >& theta_sign,
				 rng<RNGTYPE>& stream){
  
  const unsigned int J = theta.rows(); 
  const Matrix<> I = eye<double>(2);  
  for (unsigned int j = 0; j < J; ++j){
    Matrix<> X(theta_n[j],2);
    Matrix<> z(theta_n[j],1);

    // case 1: no equality constraints at all
    if (theta_eq(j,0) == -999 && theta_eq(j,1) == -999){
      for (unsigned int i = 0; i < theta_n[j]; ++i){	  
	X(i,0)= std::cos(*theta_gamma_ptr[j][i]) * theta_sign[j][i];
	X(i,1)= std::sin(*theta_gamma_ptr[j][i]) * theta_sign[j][i];
	z(i,0) = *theta_Ystar_ptr[j][i] + theta_sign[j][i] * 
	  (std::cos(*theta_gamma_ptr[j][i]) * *theta_theta_ptr[j][i][0]+ 
	   std::sin(*theta_gamma_ptr[j][i]) * *theta_theta_ptr[j][i][1]);
      }
      Matrix<> v = invpd(crossprod(X) + I); 
      Matrix<> m = v * t(X)*z;
      Matrix<> v_C = cholesky(v);      
      // case 1a: no inequality constraints at all
      if (theta_ineq(j,0) == 0 && theta_ineq(j,1) == 0){
	Matrix<> theta_j_free = gaxpy(v_C, stream.rnorm(2, 1, 0, 1),m);
	theta(j, _) =  t(theta_j_free);
      }
      else{ // case 1b: some inequality constraint holds
	Matrix<> theta_j_free = gaxpy(v_C, stream.rnorm(2, 1, 0, 1),m);
	double ineq_holds = 0.0;
        for (unsigned int jj = 0; jj < 2; ++jj){			
	  ineq_holds = std::min(ineq_holds,
				theta_ineq(j,jj) * theta_j_free(jj)); 
	}		
	while (ineq_holds < 0){
	  theta_j_free = gaxpy(v_C, stream.rnorm(2, 1, 0, 1), m);
	  double test = 0;
	  for (unsigned int jj = 0; jj < 2; ++jj) {
	    test = std::min(test, theta_ineq(j,jj) * theta_j_free(jj)); 
	  } 
	  ineq_holds = test; 
        }
	theta(j,_) = t(theta_j_free);  
      }      
    }
    else if (theta_eq(j,0) != -999 && theta_eq(j,1) == -999){ 
    // case 2: equality constraint on theta_j1 but not on theta_j2
      double xx = 0.0;
      double xz = 0.0;
      for (unsigned int i = 0; i < theta_n[j]; ++i){
	const double x_ji = std::sin(*theta_gamma_ptr[j][i]) * theta_sign[j][i];
	const double z_ji = *theta_Ystar_ptr[j][i] + theta_sign[j][i] *
	  (std::cos(*theta_gamma_ptr[j][i]) * *theta_theta_ptr[j][i][0]+ 
	   std::sin(*theta_gamma_ptr[j][i]) * *theta_theta_ptr[j][i][1]) + -1*theta_sign[j][i] * theta_eq(j,0) * std::cos(*theta_gamma_ptr[j][i]); 
	xx += x_ji * x_ji;
	xz += x_ji * z_ji;
      }
      const double v = 1.0/(xx + 1.0); 
      const double m = v * (xz);	
      // case 2a: no inequality constraint on theta_j2
      if (theta_ineq(j,1) == 0){
	theta(j,1) =  stream.rnorm(m, std::sqrt(v)); 
      }
      else{ // case 2b: inequality constraint on theta_j2
	if (theta_ineq(j,1) > 0) { // theta_j2 > 0
	  theta(j,1) = stream.rtbnorm_combo(m, v, 0);  
	} else { // theta_j2 < 0
	  theta(j,1) = stream.rtanorm_combo(m, v, 0);  	  
	}
      }
    }
    else if (theta_eq(j,0) == -999 && theta_eq(j,1) != -999){
    // case 3: equality constraint on theta_j2 but not on theta_j1
      double xx = 0.0;
      double xz = 0.0;
      for (unsigned int i = 0; i < theta_n[j]; ++i){
	const double x_ji = std::cos(*theta_gamma_ptr[j][i]) * theta_sign[j][i];
	const double z_ji = *theta_Ystar_ptr[j][i] + theta_sign[j][i] *
	  (std::cos(*theta_gamma_ptr[j][i]) * *theta_theta_ptr[j][i][0]+ 
	   std::sin(*theta_gamma_ptr[j][i]) * *theta_theta_ptr[j][i][1]) + -1*theta_sign[j][i] * theta_eq(j,1) * std::sin(*theta_gamma_ptr[j][i]); 
	xx += x_ji * x_ji;
	xz += x_ji * z_ji;
      }
      const double v = 1.0/(xx + 1.0); 
      const double m = v * (xz);	
      // case 3a: no inequality constraint on theta_j1
      if (theta_ineq(j,0) == 0){
	theta(j,0) =  stream.rnorm(m, std::sqrt(v)); 
      }
      else{ // case 3b: inequality constraint on theta_j1
	if (theta_ineq(j,0) > 0) { // theta_j1 > 0
	  theta(j,0) = stream.rtbnorm_combo(m, v, 0);  
	} else { // theta_j1 < 0
	  theta(j,0) = stream.rtanorm_combo(m, v, 0);  	  
	}
      }
    }
    else if (theta_eq(j,0) != -999 && theta_eq(j,1) != -999){
    // case 4: equality constraints on both theta_j1 and theta_j2
      theta(j,0) = theta_eq(j,0);
      theta(j,1) = theta_eq(j,1);
    }
   
  } // end j loop
  
}


// update gamma values in 2-d paired comparisons model
template <typename RNGTYPE>
void paircompare2d_gamma_update (Matrix<>& gamma,
				 const Matrix<unsigned int>& gamma_n,
				 const vector< vector < double* > >& gamma_Ystar_ptr,
				 const vector< vector < vector< double* > > >& gamma_theta1_ptr,
				 const vector< vector < vector< double* > > >& gamma_theta2_ptr,
				 // const vector< vector < int > >& gamma_choice,
				 const double& tune, // A random walk step in the MH sampler : step \sim runif(-tune, tune). 
				 vector<double>& gamma_trial,
				 vector<double>& gamma_accept,
				 rng<RNGTYPE>& stream){
  
  const unsigned int I = gamma.rows();  
  for (unsigned int i = 0; i < I; ++i){
    double gamma_i=gamma(i);
    double gamma_i_new=gamma(i)+(-2*stream.runif()+1)*tune;
    //I use 1.5707959999999 to denote PI/2
    while(gamma_i_new<0 || gamma_i_new>1.5707959999999){
      gamma_i_new=gamma(i)+(-2*stream.runif()+1)*tune;
    }
    //I assume a uniform prior for gamma, so I only consider likelihood for MH rejection rate.
    double log_lik_old=0;
    double log_lik_new=0;
    double eta_i_j=0;//eta_i_j is the linear function value within the Probit link for judge i's j'th comparison. 
    double eta_i_j_new=0; 
	
    for (unsigned int j = 0; j < gamma_n[i]; ++j){
      eta_i_j=std::cos(gamma_i)* *gamma_theta1_ptr[i][j][0] + std::sin(gamma_i)**gamma_theta1_ptr[i][j][1] -
	std::cos(gamma_i)* *gamma_theta2_ptr[i][j][0]  - std::sin(gamma_i)* *gamma_theta2_ptr[i][j][1] ;
      log_lik_old += lndnorm(*gamma_Ystar_ptr[i][j], eta_i_j, 1);//lndnorm is the log scale normal density
		
      eta_i_j_new=std::cos(gamma_i_new)* *gamma_theta1_ptr[i][j][0] + std::sin(gamma_i_new)**gamma_theta1_ptr[i][j][1] -
	std::cos(gamma_i_new)* *gamma_theta2_ptr[i][j][0]  - std::sin(gamma_i_new)* *gamma_theta2_ptr[i][j][1] ;
      log_lik_new += lndnorm(*gamma_Ystar_ptr[i][j], eta_i_j_new, 1);
    }
    gamma_trial[i]+=1;//A new value has been proposed for gamma i
    if(stream.runif()<std::exp(log_lik_new-log_lik_old)){
      gamma(i) = gamma_i_new;
      gamma_accept[i]+=1;//A new value has been accepted for gamma i
    }   
  }
}



      


template <typename RNGTYPE>
void MCMCpaircompare2d_impl (rng<RNGTYPE>& stream,
			     const Matrix<unsigned int>& MD, 
			     Matrix<>& theta, Matrix<>& gamma,
			     const Matrix<>& theta_eq, 
			     const Matrix<>& theta_ineq,
			     const double tune,
			     const unsigned int burnin,
			     const unsigned int mcmc, 
			     const unsigned int thin,
			     const unsigned int verbose,
			     const bool storegamma, const bool storetheta, 
			     double* sampledata,
			     const unsigned int samplesize,
			     double* gammaacceptrate){


  // constants
  const unsigned int N = MD.rows();//data matrix
  const unsigned int J = theta.rows();
  const unsigned int I = gamma.rows();
  const unsigned int tot_iter = burnin + mcmc;  
  const unsigned int nsamp = mcmc / thin;
  
  
  // starting values for Ystar
  Matrix<> Ystar(N,1);
  

  // pre-compute what can be pre-computed
  Matrix<unsigned int> theta_n(J,1,true,0); // vector of 0s. Item j's total instances of being compared.
  Matrix<unsigned int> gamma_n(I,1,true,0); // vector of 0s. Judge i's total instances of making comparisons.
  for (unsigned int i = 0; i < N; ++i){
    gamma_n(MD(i,0)) += 1;
    theta_n(MD(i,1)) += 1;
    theta_n(MD(i,2)) += 1;
  }

  vector< vector< double* > > theta_Ystar_ptr;
  vector< vector< double* > > theta_gamma_ptr;
  vector< vector< vector< double* > > > theta_theta_ptr;//the inner vector has two pointers for theta1 and theta2
  vector< vector< double > > theta_sign;

  vector< vector< double* > > gamma_Ystar_ptr;
  vector< vector< vector< double* >  > > gamma_theta1_ptr; //the inner vector has two pointers for theta1 and theta2
  vector< vector< vector< double* >  > > gamma_theta2_ptr; //the inner vector has two pointers for theta1 and theta2
  vector< vector< int > > gamma_choice; 
  // One inner vector of the above vector of vectors stores judge i's choice between all the pairs he works on.
  // 1 indicates picking c1, 0 indicates picking c2
  vector<double> gamma_trial;
  vector<double> gamma_accept;
  
  theta_Ystar_ptr.reserve(J);
  theta_gamma_ptr.reserve(J);
  theta_theta_ptr.reserve(J);
  theta_sign.reserve(J);

  gamma_Ystar_ptr.reserve(I);
  gamma_theta1_ptr.reserve(I);
  gamma_theta2_ptr.reserve(I);
  // gamma_choice.reserve(I);
  gamma_trial.reserve(I);
  gamma_accept.reserve(I);
  for(unsigned int i=0; i<I; ++i){
    gamma_trial.push_back(0);
    gamma_accept.push_back(0);
  }
  

  for (unsigned int j = 0; j < J; ++j){
    vector< double* > Ystar_j_ptr;
    vector< double* > gamma_j_ptr;
    vector< vector <double*> > theta_j_ptr;
    vector< double > sign_j;
	 
    Ystar_j_ptr.reserve(theta_n(j));
    gamma_j_ptr.reserve(theta_n(j));
    theta_j_ptr.reserve(theta_n(j));
    sign_j.reserve(theta_n(j));
    /*
      Allocate proper length for the empty vectors (Ystar_j_ptr, gamma_j_ptr, theta_j_ptr, sign_j)
      Then we save these pre-allocated vectors inside the outer layer vectors.
    */
    theta_Ystar_ptr.push_back(Ystar_j_ptr);
    theta_gamma_ptr.push_back(gamma_j_ptr);
    theta_theta_ptr.push_back(theta_j_ptr);
    theta_sign.push_back(sign_j);
  }

  for (unsigned int i = 0; i < I; ++i){
    vector< double* > Ystar_i_ptr;
    vector< vector< double* > > theta1_i_ptr;
    vector< vector< double* > > theta2_i_ptr;
    // vector< int > choice_i;

    Ystar_i_ptr.reserve(gamma_n(i));
    theta1_i_ptr.reserve(gamma_n(i));
    theta2_i_ptr.reserve(gamma_n(i));
    // choice_i.reserve(gamma_n(i));
    /*
      Allocate proper length for the empty vectors (Ystar_j_ptr, theta1_i_ptr, theta2_i_ptr)
      Then we save these pre-allocated vectors inside the outer layer vectors.
    */
    gamma_Ystar_ptr.push_back(Ystar_i_ptr);
    gamma_theta1_ptr.push_back(theta1_i_ptr);
    gamma_theta2_ptr.push_back(theta2_i_ptr);
    // gamma_choice.push_back(choice_i);
  }

  
  for (unsigned int i = 0; i < N; ++i){
    unsigned int resp = MD(i,0); 
    unsigned int c1 = MD(i,1);
    unsigned int c2 = MD(i,2);
    /*
      Assign the starting values of Ystar, gamma, theta and sign to the above vectors.
      This is convenient for accessing wanted objects in later updation.
    */
    theta_Ystar_ptr[c1].push_back(&Ystar(i,0));
    theta_gamma_ptr[c1].push_back(&gamma(MD(i,0)));
    vector<double*> theta_temp1;
    theta_temp1.push_back(&theta(MD(i,2),0));
    theta_temp1.push_back(&theta(MD(i,2),1));
    theta_theta_ptr[c1].push_back(theta_temp1);
    theta_sign[c1].push_back(1.0);
    
    theta_Ystar_ptr[c2].push_back(&Ystar(i,0));
    theta_gamma_ptr[c2].push_back(&gamma(MD(i,0)));
    vector<double*> theta_temp2;
    theta_temp2.push_back(&theta(MD(i,1),0));
    theta_temp2.push_back(&theta(MD(i,1),1));
    theta_theta_ptr[c2].push_back(theta_temp2);
    theta_sign[c2].push_back(-1.0);

    gamma_Ystar_ptr[resp].push_back(&Ystar(i,0));
    gamma_theta1_ptr[resp].push_back(theta_temp2);//theta_temp2 represents c1
    gamma_theta2_ptr[resp].push_back(theta_temp1);//theta_temp1 represents c2
	
	
    // const bool above_i = MD(i,1) == MD(i,3); // cand 1 chosen
    // const bool below_i = MD(i,2) == MD(i,3); // cand 2 chosen
    // if (above_i){
      // gamma_choice[resp].push_back(1);
    // }
    // else if (below_i){
      // gamma_choice[resp].push_back(0);
    // }
	
	
    theta_temp1.clear();
    theta_temp2.clear();
  }


  
  
  
  // storage matrices (col major order)
  Matrix<> theta_store;
  Matrix<> gamma_store;
  if (storetheta)
    theta_store = Matrix<>(nsamp, J*2);
  
  if (storegamma)
    gamma_store = Matrix<>(nsamp, I);

  

  unsigned int count = 0;
  // MCMC sampling occurs in this for loop
  for (unsigned int iter = 0; iter < tot_iter; ++iter){

    // The following three functions are defined in MCMCfcds.h
    // sample Ystar
    paircompare2d_Ystar_update(Ystar, MD, theta, gamma, stream);

    // sample gamma
    paircompare2d_gamma_update(gamma,
			       gamma_n,
			       gamma_Ystar_ptr,
			       gamma_theta1_ptr,
			       gamma_theta2_ptr,
			       //gamma_choice,
			       tune, 
			       gamma_trial,
			       gamma_accept,
			       stream);
    

    // sample theta
    paircompare2d_theta_update(theta, Ystar, MD, gamma, theta_n, theta_eq,
			       theta_ineq,
			       theta_Ystar_ptr,
			       theta_gamma_ptr,
			       theta_theta_ptr,
			       theta_sign,
			       stream);
	
    
    // print results to screen
    if (verbose > 0 && iter % verbose == 0) {
      Rprintf("\n\nMCMCpaircompare2d iteration %i of %i \n", (iter+1), tot_iter);
      Rprintf("\nSummary of acceptance rates for gamma:\n");
      double gamma_accept_min = 1.0;
      double gamma_accept_max = 0.0;
      double gamma_accept_num = 0.0;
      double gamma_accept_denom = 0.0;
      for (unsigned int i=0; i<I; ++i){
	double local_accrate = gamma_accept[i]/gamma_trial[i];
	gamma_accept_min = (local_accrate < gamma_accept_min) ? local_accrate : gamma_accept_min;
	gamma_accept_max = (local_accrate > gamma_accept_max) ? local_accrate : gamma_accept_max;
	gamma_accept_num += local_accrate;
	gamma_accept_denom += 1.0;
      }
      Rprintf("Minimum gamma acceptance rate: %7.3f", gamma_accept_min, " ");
      Rprintf("\nMean gamma acceptance rate:    %7.3f", gamma_accept_num/gamma_accept_denom, " ");
      Rprintf("\nMaximum gamma acceptance rate: %7.3f", gamma_accept_max, " ");
      Rprintf("\n\n");
     
    }
    
    // store results
    if ((iter >= burnin) && ((iter % thin == 0))) {
      
      // store theta
      if (storetheta) {
	for (unsigned int j = 0; j < J; ++j) {
	  theta_store(count, j) = theta(j,0);
	  theta_store(count, J+j) = theta(j, 1);
	}
      }
        

      // store gamma
      if (storegamma)
        gamma_store(count, _) = gamma;

      count++;	
    }
    
    R_CheckUserInterrupt(); // allow user interrupts
       
  } // end MCMC sampling loop
 


  
  // put output back into sampledata
  Matrix<> output;
  if(storetheta && ! storegamma) {
    output = theta_store;
  } else if (storegamma && ! storetheta){
    output = gamma_store;
  } else {
    output = cbind(theta_store, gamma_store);
  }

  for (unsigned int i = 0; i < samplesize; ++i)
    sampledata[i] = output[i];

  for (unsigned int i = 0; i < I; ++i)
    gammaacceptrate[i] = gamma_accept[i]/gamma_trial[i];

} // end MCMCpaircompare implementation









extern "C" {

  void
  cMCMCpaircompare2d(double* sampledata, 
		    const int* samplerow, 
		    const int* samplecol,
		    const unsigned int* MDdata, 
		    const int* MDrow, 
		    const int* MDcol,
		    const int* burnin, 
		    const int* mcmc,  
		    const int* thin,
		    const int *uselecuyer, 
		    const int *seedarray, 
		    const int *lecuyerstream, 
		    const int* verbose, 
		    const double* thetastartdata,
		    const int* thetastartrow, 
		    const int* thetastartcol, 
		    const double* gammastartdata, 
		    const int* gammastartrow, 
		    const int* gammastartcol,
		    const double* tunevalue,
		    const double* thetaeqdata, 
		    const int* thetaeqrow, 
		    const int* thetaeqcol,
		    const double* thetaineqdata, 
		    const int* thetaineqrow, 
		    const int* thetaineqcol, 
		    const int* storegamma, 
		    const int* storetheta,
		    double* gammaacceptrate){


    // put together matrices
    const Matrix<unsigned int> MD(*MDrow, *MDcol, MDdata);
    Matrix<> theta(*thetastartrow, *thetastartcol, thetastartdata);
    Matrix<> gamma(*gammastartrow, *gammastartcol, gammastartdata);
    const Matrix<> theta_eq(*thetaeqrow, *thetaeqcol, thetaeqdata);
    const Matrix<> theta_ineq(*thetaineqrow, *thetaineqcol, thetaineqdata);
    const int samplesize = (*samplerow) * (*samplecol);

    
    MCMCPACK_PASSRNG2MODEL(MCMCpaircompare2d_impl, MD, theta, gamma, 
			   theta_eq, theta_ineq, *tunevalue,
			   *burnin, *mcmc, *thin, 
			   *verbose, *storegamma, *storetheta, 
			   sampledata, samplesize, gammaacceptrate);
  }
}


#endif










