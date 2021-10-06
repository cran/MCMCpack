//////////////////////////////////////////////////////////////////////////
// MCMCpaircompare.cc is C++ code to estimate a pairwise comparison model. 
//
// KQ 3/18/2015
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////

#ifndef MCMCPAIRCOMPARE2DDP_CC
#define MCMCPAIRCOMPARE2DDP_CC

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

/* MCMCpaircompare2dDP implementation. */

// update latent Ystar values in 2-d paired comparisons model
template <typename RNGTYPE>
void paircompare2dDP_Ystar_update (Matrix<>& Ystar,
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
void paircompare2dDP_theta_update (Matrix<>& theta, const Matrix<>& Ystar, 
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

// update cluster_gamma  in 2-d Dirichlet Process paired comparisons model
template <typename RNGTYPE>
void paircompare2dDP_cluster_gamma_update (//Matrix<>& gamma,
					   const Matrix<unsigned int>& gamma_n,
					   const vector< vector < double* > >& gamma_Ystar_ptr,
					   const vector< vector < vector< double* > > >& gamma_theta1_ptr,
					   const vector< vector < vector< double* > > >& gamma_theta2_ptr,
					   const double& tune, // A random walk step in the MH sampler : step \sim runif(-tune, tune). 
					   const unsigned int& clustermcmc,
					   vector<double>& gamma_trial,
					   vector<double>& gamma_accept,				 
					   vector<unsigned int>& judge_cluster_membership, //a vector of length I that stores judge i's cluster_membership
					   vector<double>& cluster_gamma,//a vector of length K that stores the gamma value that charaterizes cluster k				 
					   vector<unsigned int>& cluster_size,
					   rng<RNGTYPE>& stream){
  
  const unsigned int I = judge_cluster_membership.size();  
  const unsigned int K = cluster_gamma.size();
  for (unsigned int k = 0; k < K; ++k){
    // if cluster k has no member for now, we generate the new k from its prior
    if(cluster_size[k]==0){
      //I use 1.5707959999999 to denote PI/2
      cluster_gamma[k]=stream.runif()*1.5707959999999;//directly sampling from unif(0,PI/2)
    }else{ //if cluster k has current members, we generate the new k through HM sampler
      // we update gamma_k with the last value of this smaller HM sampler.
      double gamma_k=cluster_gamma[k];
      for(unsigned int t=0;t<clustermcmc;++t){
	double gamma_k_new=gamma_k+(-2*stream.runif()+1)*tune;
	//I use 1.5707959999999 to denote PI/2
	while(gamma_k_new<0 || gamma_k_new>1.5707959999999){
	  gamma_k_new=gamma_k+(-2*stream.runif()+1)*tune;
	}
	//I assume a uniform prior for gamma, so I only consider likelihood for MH rejection rate.
	double log_lik_old=0;
	double log_lik_new=0;
	double eta_i_j=0;//eta_i_j is the linear function value within the Probit link for judge i's j'th comparison. 
	double eta_i_j_new=0; 
	for(unsigned int i=0;i<I;++i){
	  if(judge_cluster_membership[i]==k){ //only use cluster k's current member judges' likelihood 
	    gamma_trial[i]+=1;//A new value has been proposed for gamma i, if judge i is in cluster k
	    
	    for (unsigned int j = 0; j < gamma_n[i]; ++j){
	      eta_i_j=std::cos(gamma_k)* *gamma_theta1_ptr[i][j][0] + std::sin(gamma_k)**gamma_theta1_ptr[i][j][1] -
		std::cos(gamma_k)* *gamma_theta2_ptr[i][j][0]  - std::sin(gamma_k)* *gamma_theta2_ptr[i][j][1] ;
	      log_lik_old += lndnorm(*gamma_Ystar_ptr[i][j], eta_i_j, 1);//lndnorm is the log scale normal density
	      
	      eta_i_j_new=std::cos(gamma_k_new)* *gamma_theta1_ptr[i][j][0] + std::sin(gamma_k_new)**gamma_theta1_ptr[i][j][1] -
		std::cos(gamma_k_new)* *gamma_theta2_ptr[i][j][0]  - std::sin(gamma_k_new)* *gamma_theta2_ptr[i][j][1] ;
	      log_lik_new += lndnorm(*gamma_Ystar_ptr[i][j], eta_i_j_new, 1);
	    }
	  }
	}
	
	if(stream.runif()<std::exp(log_lik_new-log_lik_old)){
	  gamma_k = gamma_k_new;
	  for(unsigned int i=0;i<I;++i){
	    if(judge_cluster_membership[i]==k){
	      gamma_accept[i]+=1;//A new value has been accepted for gamma i
	    }
	  } 
	}	
      }
      cluster_gamma[k]=gamma_k;
    }//work on cluster k
  }//loop through all clusters
}


// update cluster_weight_log  in 2-d Dirichlet Process paired comparisons model
template <typename RNGTYPE>
void paircompare2dDP_cluster_weight_log_update (
						//vector<double>& cluster_gamma,//a vector of length K that stores the gamma value that charaterizes cluster k
						vector<double>& cluster_weight_log,//on log scale
						vector<unsigned int>& cluster_size,
						const double& alpha,
						const unsigned int& I,
						rng<RNGTYPE>& stream){
  //IMPORTANT: cluster_weight_log is on the log scale
  
  
  //work on cluster 1 through K-1:
  const unsigned int K=cluster_weight_log.size();
  double cluster_cumulative_size= (double) I;
  double weight_log=0;
  double V_log=0;
  double cumulative_one_minue_V_log=0;//must be initialized at zero
  double rbeta_param1=0;
  double rbeta_param2=0;
  double beta_rv=0;
  for(unsigned int k=0;k<K-1;++k){
    rbeta_param1=1+cluster_size[k];
    cluster_cumulative_size -= cluster_size[k];
    rbeta_param2=alpha+cluster_cumulative_size;
    beta_rv=stream.rbeta(rbeta_param1,rbeta_param2);
    while(beta_rv>0.9999){
      beta_rv=stream.rbeta(rbeta_param1,rbeta_param2);
    }
    V_log=std::log(beta_rv);
    weight_log=V_log+cumulative_one_minue_V_log;
    cumulative_one_minue_V_log+=std::log(1-beta_rv);
    cluster_weight_log[k]=weight_log;
  }
  //work on cluster K indexed by K-1
  cluster_weight_log[K-1]=cumulative_one_minue_V_log;
}



// x: n array of original indices running from 0 to (n-1)
// n: length of x 
// sample one element from x with equal likelihood
template <typename RNGTYPE>
unsigned int SampleOneEqualLikelihood(int n, vector<unsigned int>& x, rng<RNGTYPE>& stream) {
	//for (int i = 0; i < n; i++)
	//x[i] = i;
	//for (int i = 0; i < k; i++) {
	double u = stream.runif();//Rcpp::runif(1, 0, 1)[0];//runif(1,0,1) has a NumericVecor output, so I need to index to get the double number out
	int j = static_cast<int>(n * u);
	return x[j];
		// y[i] = x[j];
		// x[j] = x[--n];
	//}
}

// draws 1 element from a const vector<unsigned int> x with prob weights
// in array p
// n is length of x and p 
template <typename RNGTYPE>
int ProbSamp(const vector<unsigned int>& x, const vector< double >& p,
	     const unsigned int& n, rng<RNGTYPE>& stream) {
  const double u = stream.runif();
  double cumprob = p[0];
  for (unsigned int i = 0; i < (n - 1); ++i) {
		if (u <= cumprob) {
		  return x[i];
		}
		else {
		  cumprob += p[i + 1];
		}
  }
  return x[n - 1];
}


// update judge_cluster_membership  in 2-d Dirichlet Process paired comparisons model
template <typename RNGTYPE>
void paircompare2dDP_judge_cluster_membership_update (//Matrix<>& gamma,
	    const Matrix<unsigned int>& gamma_n,
	    const vector< vector < double* > >& gamma_Ystar_ptr,
	    const vector< vector < vector< double* > > >& gamma_theta1_ptr,
	    const vector< vector < vector< double* > > >& gamma_theta2_ptr,
	    vector<unsigned int>& judge_cluster_membership, //a vector of length I that stores judge i's cluster_membership
	    vector<double>& cluster_gamma,//a vector of length K that stores the gamma value that charaterizes cluster k
	    vector<double>& cluster_weight_log,//on log scale
            vector<unsigned int>& cluster_size,
	    const vector<unsigned int>& cluster_labels,
	    unsigned int& unique_cluster,
	    rng<RNGTYPE>& stream){
  
  //IMPORTANT: I use the log-sum-exp trick to control underflow or overflow
  // I store log-scale cluster probabiltiy values  in cluster_probability_log.
  // Moreover, I substract the max log-probabiltiy from every log-probability
  
  const unsigned int I = judge_cluster_membership.size();
  const unsigned int K = cluster_weight_log.size();
  vector< double > cluster_probability;
  vector< double > cluster_probability_log;
  cluster_probability.reserve(K);
  cluster_probability_log.reserve(K);


  
  
  for (unsigned int i = 0; i < I; ++i){
    double max_cluster_prob_log=-100000000000;
    unsigned int cluster_i_new=0;
    double denominator=0;
    for(unsigned int k=0; k< K; ++k){
      //I assume a uniform prior for gamma, so I only consider likelihood for MH rejection rate.
      double log_lik_i_k=0; //the log_lik of judge i being in cluster k
      double eta_i_j=0;//eta_i_j is the linear function value within the Probit link for judge i's j'th comparison. 
      double gamma_i=cluster_gamma[k];//we let gamma_i try on each value in cluster_gamma
      double cluster_prob_log=0;
      for (unsigned int j = 0; j < gamma_n[i]; ++j){
	eta_i_j=std::cos(gamma_i)* *gamma_theta1_ptr[i][j][0] + std::sin(gamma_i)**gamma_theta1_ptr[i][j][1] -
	  std::cos(gamma_i)* *gamma_theta2_ptr[i][j][0]  - std::sin(gamma_i)* *gamma_theta2_ptr[i][j][1] ;
	log_lik_i_k += lndnorm(*gamma_Ystar_ptr[i][j], eta_i_j, 1);//lndnorm is the log scale normal density
	
      }
      cluster_prob_log=cluster_weight_log[k]+log_lik_i_k;//both on log scale
      cluster_probability_log[k]=cluster_prob_log;
      if(cluster_prob_log>max_cluster_prob_log){
	max_cluster_prob_log=cluster_prob_log;
      }
    }//finish loop for computing cluster_prob_log for judge i
    
    //do the loop below to prevent underflow or overflow
    for(unsigned int kk=0; kk< K; ++kk){
      cluster_probability_log[kk]-=max_cluster_prob_log;
    }
    //compute the cluster probability
    
    for(unsigned int kkk=0; kkk< K; ++kkk){
      denominator+=std::exp(cluster_probability_log[kkk]);
    }
    
    
    for(unsigned int kkkk=0; kkkk< K; ++kkkk){
      cluster_probability[kkkk]=std::exp(cluster_probability_log[kkkk])/denominator;
    }
    
    
    
    cluster_i_new = ProbSamp(cluster_labels, cluster_probability, K, stream);
    if(cluster_i_new != judge_cluster_membership[i]){
      
      //bookkeeping for the number of unique clusters:
      if(cluster_size[cluster_i_new]==0 && cluster_size[judge_cluster_membership[i]]>1){
	//new cluster used to be empty, and the current cluster has more than one member
	++unique_cluster;
      }else if(cluster_size[cluster_i_new]>0 && cluster_size[judge_cluster_membership[i]]==1){
	//new cluster already has members, and the current cluster will be empty
	--unique_cluster;
      }//else: new cluster already has members,and the current cluster has more than one member. 
      //Don't need to change unique_cluster
      
      cluster_size[judge_cluster_membership[i]]-=1;
      cluster_size[cluster_i_new]+=1;
      judge_cluster_membership[i]=cluster_i_new;
    }//else: the new cluster label for i is the same as the old one, we don't need to do anything	
  }
  
}




// update alpha  in 2-d Dirichlet Process paired comparisons model
template <typename RNGTYPE>
void paircompare2dDP_alpha_update (
				 vector<double>& cluster_weight_log,//on log scale. I use the last element to update alpha
				 double& alpha,
				 const unsigned int& K,
				 const double& a,
				 const double& b,
				 rng<RNGTYPE>& stream){
	
	// shape=a+K-1;
	// rate=b-cluster_weight_log[K-1];
	// cluster_weight_log is already on the log scale
	// Rprintf("\nshape is: %7.3f", a+K-1, " ");
	// Rprintf("\nrate is:    %7.3f", b-cluster_weight_log[K-1], " ");
	alpha=stream.rgamma(a+K-1, b-cluster_weight_log[K-1]);
	// if(alpha < 0.1){
		// alpha=stream.rgamma(a+K-1, b-cluster_weight_log[K-1]);
	// }

}



template <typename RNGTYPE>
void MCMCpaircompare2dDP_impl (rng<RNGTYPE>& stream,
			     const Matrix<unsigned int>& MD, 
			     Matrix<>& theta, Matrix<>& gamma,
				 Matrix<>& cluster_gamma_mat,
				 Matrix<unsigned int>& judge_cluster_membership_mat,
			     const Matrix<>& theta_eq, 
			     const Matrix<>& theta_ineq,
			     const double tune,
			     const unsigned int burnin,
			     const unsigned int mcmc,
			     const unsigned int clustermcmc,				 
			     const unsigned int thin,
			     const unsigned int verbose,
			     const bool storegamma, const bool storetheta, 
			     double* sampledata,
			     const unsigned int samplesize,
			     double* gammaacceptrate,
				 double alpha,
				 const unsigned int clustermax,
				 const int alphafixed,
				 const double a, const double b){


  // constants
  const unsigned int N = MD.rows();//data matrix
  const unsigned int J = theta.rows();
  const unsigned int I = gamma.rows();
  const unsigned int tot_iter = burnin + mcmc;  
  const unsigned int nsamp = mcmc / thin;
  const unsigned int K = clustermax;//the maximum number of clusters
  
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
  vector<unsigned int> judge_cluster_membership;//a vector of length I that stores judge i's cluster_membership
  vector<double> cluster_gamma;//a vector of length K that stores the gamma value that charaterizes cluster k
  vector<double> cluster_weight_log;
  vector<unsigned int> cluster_size;
  vector<unsigned int> cluster_labels;
  
  
  //vector<double> gamma_accept_rate_all;
  vector<double> gamma_trial;
  vector<double> gamma_accept;
  
  theta_Ystar_ptr.reserve(J);
  theta_gamma_ptr.reserve(J);
  theta_theta_ptr.reserve(J);
  theta_sign.reserve(J);

  gamma_Ystar_ptr.reserve(I);
  gamma_theta1_ptr.reserve(I);
  gamma_theta2_ptr.reserve(I);
  judge_cluster_membership.reserve(I);
  cluster_gamma.reserve(K);
  cluster_weight_log.reserve(K);
  cluster_size.reserve(K);
  
  
  //gamma_trial and gamma_accept record MH sampler proposal and accept for each cluster
  gamma_trial.reserve(I);
  gamma_accept.reserve(I);
  
  
  for(unsigned int k=0; k<K; ++k){
	cluster_size.push_back(0);
	}
  //read in the supplied starting values for judge_cluster_membership
  //record each cluster's size based on the starting values
  for(unsigned int i=0; i<I; ++i){
	gamma_trial.push_back(0);
    gamma_accept.push_back(0);
	judge_cluster_membership.push_back( judge_cluster_membership_mat(i) );
	cluster_size[judge_cluster_membership_mat(i)]+=1;
  }
  
  unsigned int unique_cluster=set<double>( judge_cluster_membership.begin(), judge_cluster_membership.end() ).size();
  
  
  for(unsigned int k=0; k<K; ++k){
    
    cluster_gamma.push_back(cluster_gamma_mat(k));
	cluster_weight_log.push_back(0);
	cluster_labels.push_back(k);
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
	
    theta_temp1.clear();
    theta_temp2.clear();
  }


  
  
  
  // storage matrices (col major order)
  Matrix<> theta_store;
  Matrix<> gamma_store;
  Matrix<unsigned int> gamma_cluster_store;
  Matrix<unsigned int> unique_cluster_store;
  Matrix<> alpha_store;
  if (storetheta)
    theta_store = Matrix<>(nsamp, J*2);
  
  if (storegamma){
	  gamma_store = Matrix<>(nsamp, I);
	  gamma_cluster_store= Matrix<unsigned int>(nsamp, I);
	  unique_cluster_store= Matrix<unsigned int>(nsamp, 1);
  }
  
  if(alphafixed==0){
	  alpha_store= Matrix<>(nsamp, 1);
  }
    

	

  

  unsigned int count = 0;
  // MCMC sampling occurs in this for loop
  for (unsigned int iter = 0; iter < tot_iter; ++iter){

    
    // sample Ystar
    paircompare2dDP_Ystar_update(Ystar, MD, theta, gamma, stream);
  
   
	paircompare2dDP_cluster_gamma_update(
				 gamma_n,
				 gamma_Ystar_ptr,
				 gamma_theta1_ptr,
				 gamma_theta2_ptr,
				 tune, // A random walk step in the MH sampler : step \sim runif(-tune, tune). 
			     clustermcmc,
				 gamma_trial,
				 gamma_accept,
				 judge_cluster_membership, //a vector of length I that stores judge i's cluster_membership
				 cluster_gamma,//a vector of length K that stores the gamma value that charaterizes cluster k
                 cluster_size,
				 stream);
	
	paircompare2dDP_cluster_weight_log_update (
				 cluster_weight_log,
                 cluster_size,
				 alpha,
				 I,
				 stream);
	
	paircompare2dDP_judge_cluster_membership_update (
				 gamma_n,
				 gamma_Ystar_ptr,
				 gamma_theta1_ptr,
				 gamma_theta2_ptr,
				 judge_cluster_membership, //a vector of length I that stores judge i's cluster_membership
				 cluster_gamma,//a vector of length K that stores the gamma value that charaterizes cluster k
				 cluster_weight_log,//on log scale
                 cluster_size,
				 cluster_labels,
				 unique_cluster,
				 stream);
	
	// update the gamma_i's according to the new judge_cluster_membership
	for(unsigned int i=0;i<I;++i){
		gamma(i)=cluster_gamma[judge_cluster_membership[i]];
	}
	

    // sample theta
    paircompare2dDP_theta_update(theta, Ystar, MD, gamma, theta_n, theta_eq,
			       theta_ineq,
			       theta_Ystar_ptr,
			       theta_gamma_ptr,
			       theta_theta_ptr,
			       theta_sign,
			       stream);
	
    
	if (alphafixed == 0){
		paircompare2dDP_alpha_update (
				 cluster_weight_log,//on log scale. I use the last element to update alpha
				 alpha,
				 K,
				 a,
				 b,
				 stream); 
		
	}
	
	// if(cluster_size[K-1]==0){
		// vector<unsigned int> candidates;
		// for(unsigned int k=0;k<K;++k){
			// if(cluster_size[k]!=0){
				// candidates.push_back(k);
			// }
		// }
		// unsigned int candidate=SampleOneEqualLikelihood(candidates.size(), candidates, stream);
		// Rprintf("\n\ncandidate: %i\n", candidate);
		// //switch everything between cluster candidate and cluster K-1
		// double cluster_gamma_temp=cluster_gamma[candidate];
		// cluster_gamma[candidate]=cluster_gamma[K-1];
		// cluster_gamma[K-1]=cluster_gamma_temp;
		// double cluster_weight_log_temp=cluster_weight_log[candidate];
		// cluster_weight_log[candidate]=cluster_weight_log[K-1];
		// cluster_weight_log[K-1]=cluster_weight_log_temp;
		// unsigned int cluster_size_temp=cluster_size[candidate];
		// cluster_size[candidate]=cluster_size[K-1];
		// cluster_size[K-1]=cluster_size_temp;
		// for(unsigned int i=0;i<I;++i){
			// if(judge_cluster_membership[i]==candidate){
				// judge_cluster_membership[i]=K-1;
			// }
		// }
		
		// // // update the gamma_i's according to the new judge_cluster_membership
		// // for(unsigned int i=0;i<I;++i){
			// // gamma(i)=cluster_gamma[judge_cluster_membership[i]];
		// // }
		
	// }
	
    // print results to screen
    if (verbose > 0 && iter % verbose == 0) {
		Rprintf("\n\nMCMCpaircompare2dDP iteration %i of %i \n", (iter+1), tot_iter);
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
        if (storegamma){
			gamma_store(count, _) = gamma;
			for(unsigned int i=0; i<I; ++i){
				gamma_cluster_store(count, i)=judge_cluster_membership[i]+1;//cluster label is indexed from 1 in R. 
			}
			unique_cluster_store(count)=unique_cluster;
		}
		
		if(alphafixed==0){
			alpha_store(count)=alpha;
		}	

		count++;	
    }
    
    R_CheckUserInterrupt(); // allow user interrupts
       
  } // end MCMC sampling loop
 


  
  // put output back into sampledata
  Matrix<> output;
  if(storetheta && ! storegamma) {
	  if(alphafixed==0){
		  output = cbind(theta_store, alpha_store);
	  }else{
		  output = theta_store;
	  }
    
  } else if (storegamma && ! storetheta){
		Matrix<> output_temp1;
		Matrix<> output_temp2;
		output_temp1 = cbind(gamma_cluster_store, unique_cluster_store);
		output_temp2 = cbind(gamma_store, output_temp1);
		
		if(alphafixed==1){
		  output = output_temp2;
		}else{
		  output = cbind(output_temp2, alpha_store);
		}
    
  } else {
	Matrix<> output_temp1;
	Matrix<> output_temp2;
	Matrix<> output_temp3;
	output_temp1 = cbind(theta_store, gamma_store);
	output_temp2 = cbind(gamma_cluster_store, unique_cluster_store);
	output_temp3 = cbind(output_temp1, output_temp2);
	if(alphafixed==1){
		output = output_temp3;//cbind(gamma_store, alpha_store);
	}else{
	    output = cbind(output_temp3, alpha_store);
	}
    
  }

  for (unsigned int i = 0; i < samplesize; ++i)
    sampledata[i] = output[i];

  for (unsigned int i = 0; i < I; ++i)
    gammaacceptrate[i] = gamma_accept[i]/gamma_trial[i];

} // end MCMCpaircompare implementation









extern "C" {

  void
  cMCMCpaircompare2dDP(double* sampledata, 
		       const int* samplerow, 
		       const int* samplecol,
		       const unsigned int* MDdata, 
		       const int* MDrow, 
		       const int* MDcol,
		       const int* burnin, 
		       const int* mcmc,   
		       const int* clustermcmc,  
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
		       const double* clustergammastartdata, 
		       const int* clustergammastartrow, 
		       const int* clustergammastartcol,
		       const int* judgeclustermembershipstartdata, 
		       const int* judgeclustermembershipstartrow, 
		       const int* judgeclustermembershipstartcol,
		       const double* tunevalue,
		       const double* thetaeqdata, 
		       const int* thetaeqrow, 
		       const int* thetaeqcol,
		       const double* thetaineqdata, 
		       const int* thetaineqrow, 
		       const int* thetaineqcol, 
		       const int* storegamma, 
		       const int* storetheta,
		       double* gammaacceptrate,
		       double* alpha,
		       const unsigned int* clustermax, 
		       const int* alphafixed, 
		       const double* a,
		       const double* b){
    

    // put together matrices
    const Matrix<unsigned int> MD(*MDrow, *MDcol, MDdata);
    Matrix<> theta(*thetastartrow, *thetastartcol, thetastartdata);
    Matrix<> gamma(*gammastartrow, *gammastartcol, gammastartdata);
	Matrix<> cluster_gamma_mat(*clustergammastartrow, *clustergammastartcol, clustergammastartdata);
	Matrix<unsigned int> judge_cluster_membership_mat(*judgeclustermembershipstartrow, *judgeclustermembershipstartcol, judgeclustermembershipstartdata);
    const Matrix<> theta_eq(*thetaeqrow, *thetaeqcol, thetaeqdata);
    const Matrix<> theta_ineq(*thetaineqrow, *thetaineqcol, thetaineqdata);
    const int samplesize = (*samplerow) * (*samplecol);

    
    MCMCPACK_PASSRNG2MODEL(MCMCpaircompare2dDP_impl, MD, theta, gamma, cluster_gamma_mat,
			   judge_cluster_membership_mat, theta_eq, theta_ineq, *tunevalue,
			   *burnin, *mcmc, *clustermcmc, *thin, 
			   *verbose, *storegamma, *storetheta, 
			   sampledata, samplesize, gammaacceptrate, *alpha,
			   *clustermax, *alphafixed, *a, *b);
  }
}


#endif










