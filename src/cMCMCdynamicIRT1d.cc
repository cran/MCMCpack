//////////////////////////////////////////////////////////////////////////
// cMCMCdynamicIRT1d.cc is C++ code to estimate a dynamic 1d IRT model
//
// Kevin Quinn
// 1/29/2008
//
// Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
// Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
//    and Jong Hee Park
//////////////////////////////////////////////////////////////////////////


#ifndef CMCMCDYNAMICIRT1D_CC
#define CMCMCDYNAMICIRT1D_CC

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

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace std;
using namespace scythe;


// used to access 1d arrays holding R matrices like a 2d array
#define M(ROW,COL,NROWS) (COL*NROWS+ROW)


// cMCMCdynamicIRT1d implementation
template <typename RNGTYPE>
void MCMCdynamicIRT1d_impl(rng<RNGTYPE>& stream,
			   double* thetadraws, const int* nrowthetadraws,
			   const int* ncolthetadraws,
			   double* alphadraws, const int* nrowalphadraws,
			   const int* ncolalphadraws,
			   double* betadraws, const int* nrowbetadraws,
			   const int* ncolbetadraws,
			   double* tau2draws, const int* nrowtau2draws,
			   const int* ncoltau2draws,
			   const int* nsubj, const int* nitems,
			   const int* ntime,
			   const int* Ydata, const int* nrowYdata,
			   const int* ncolYdata,
			   const int* ITdata, const int* lengthITdata,
			   const int* burnin, const int* mcmc, const int* thin,
			   const int* verbose,
			   const double* thetadata,
			   const int* lengththeta,
			   const int* thetainfodata,
			   const int* nrowthetainfo,
			   const int* ncolthetainfo,
			   double* alphadata,
			   const int* lengthalpha,
			   double* betadata,
			   const int* lengthbeta,
			   double* tau2data,
			   const int* lengthtau2,
			   const double* c0,
			   const int* lengthc0,
			   const double* d0,
			   const int* lengthd0,
			   const double* a0,
			   const double* A0,
			   const double* b0,
			   const double* B0,
			   const double* e0,
			   const double* E0inv,
			   const double* thetaeqdata,
			   const int* nrowthetaeq,
			   const int* ncolthetaeq,
			   const double* thetaineqdata,
			   const int* nrowthetaineq,
			   const int* ncolthetaineq,
			   const int* storeitem,
			   const int* storeability){



  const int tot_iter = *burnin + *mcmc;
// JHP  const int nsamp = *mcmc / *thin;
  //sparse matrix of latent outcome variables
  double** Z;
  Z = new double*[*nsubj];
  for (int s=0; s<*nsubj; ++s){
    Z[s] = new double[*nitems];
    for (int i=0; i<*nitems; ++i){
      Z[s][i] = -999;
    }
  }
  for (int j=0; j<*nrowYdata; ++j){
    const int s = Ydata[M(j, 0, *nrowYdata)]; // subject id
    const int i = Ydata[M(j, 1, *nrowYdata)]; // item id
// JHP    const int y = Ydata[M(j, 2, *nrowYdata)]; // y value
    Z[s][i] = 0.0;
  }


  // stuff to make working with theta easy
  // theta[s][t] gives the tth theta for subject s
  // the actual time period corresponding to theta[s][t] is theta_offset[s]+t
  vector< vector<double> > theta;
  vector<int> theta_offset;
  int count = 0;
  for (int s=0; s<*nsubj; ++s){
    vector<double> holder;
    int ntime_s = thetainfodata[M(s, 1, *nrowthetainfo)];
    theta_offset.push_back(thetainfodata[M(s, 2, *nrowthetainfo)]);
    for (int t=0; t<ntime_s; ++t){
      holder.push_back(thetadata[count]);
      ++count;
    }
    theta.push_back(holder);
  }



  // define constants used in the sampling
  const double A0a0 = *A0 * *a0;
  const double B0b0 = *B0 * *b0;


  // set up mappings
  //
  // IS gives the mapping from items to subjects
  // IS[i] provides a vector of integers corresponding to the subjects
  //    who voted on item i.
  // IS[i][s] gets the subject index of the sth subject to vote on item i
  vector< vector<int> > IS;
  for (int i=0; i<*nitems; ++i){
    vector<int> subjholder;
    for (int j=0; j<*nrowYdata; ++j){
      if (Ydata[M(j,1,*nrowYdata)] == i){
	subjholder.push_back(Ydata[M(j,0,*nrowYdata)]);
      }
    }
    sort(subjholder.begin(), subjholder.end());
    IS.push_back(subjholder);
  }

  // SI gives the mapping from subjects to items
  // SI[s] provides a vector of integers corresponding to the items
  //    voted on by subject s.
  // SI[s][i] gets the item index of the ith item voted on by subject s
  vector< vector<int> > SI;
  for (int i=0; i<*nsubj; ++i){
    vector<int> itemholder;
    for (int j=0; j<*nrowYdata; ++j){
      if (Ydata[M(j,0,*nrowYdata)] == i){
	itemholder.push_back(Ydata[M(j,1,*nrowYdata)]);
      }
    }
    sort(itemholder.begin(), itemholder.end());
    SI.push_back(itemholder);
  }

  // TI gives the mapping from times to items
  // TI[t] provides a vector of integers corresponding to the items
  //    voted on in time period t.
  // TI[t][i] gets the item index of the ith item voted on in time period t
  vector< vector<int> > TI;
  for (int t=0; t<*ntime; ++t){
    vector<int> itemholder;
    for (int i=0; i<*lengthITdata; ++i){
      if (ITdata[i] == t){
	itemholder.push_back(i);
      }
    }
    sort(itemholder.begin(), itemholder.end());
    TI.push_back(itemholder);
  }

  // IT gives the mapping from items to times and is just fancy
  // way of holding the stuff in *ITdata
  vector<int> IT;
  for (int i=0; i<*nitems; ++i){
    IT.push_back(ITdata[i]);
  }


  // ST gives the mapping from subjects to times
  // ST[s] provides a vector of integers corresponding to the times
  //    in which items were voted on by subject s.
  // ST[s][t] gets the time index of the tth time period served by subject s
  /*
  vector <vector <int> > ST;
  for (int s=0; s<*nsubj; ++s){
    vector <int> timeholder;
    for (int iind=0; iind<SI[s].size(); ++iind){
      const int i = SI[s][iind];
      const int t = IT[i];
      if (timeholder.empty()){
	timeholder.push_back(t);
      }
      vector<int>::iterator myiter;
      myiter = find(timeholder.begin(), timeholder.end(), t);
      if (myiter == timeholder.end()){ // t not currently in timeholder
	timeholder.push_back(t);
      }
    }
    sort(timeholder.begin(), timeholder.end());
    ST.push_back(timeholder);
  }
  */

  // STI gives the list of item indices of items voted on by by subject s
  // in the tth time period served by s
  vector< vector < vector <int> > > STI;
  for (int s=0; s<*nsubj; ++s){
    vector < vector <int> > timeitemholder;
    for (unsigned int tt=0; tt<theta[s].size(); ++tt){
      const int t = tt + theta_offset[s];
      vector <int> itemholder;
      for (unsigned int ii=0; ii<TI[t].size(); ++ii){
	const int i = TI[t][ii];
	if (Z[s][i] != -999){
	  itemholder.push_back(i);
	}
      }
      timeitemholder.push_back(itemholder);
    }
    STI.push_back(timeitemholder);
  }





  // MCMC iterations start here
  int sampcount = 0;
  for (int iter=0; iter < tot_iter; ++iter){



    // sample latent Z
    for (int j=0; j<*nrowYdata; ++j){
      const int s = Ydata[M(j, 0, *nrowYdata)]; // subject id
      const int i = Ydata[M(j, 1, *nrowYdata)]; // item id
      const int y = Ydata[M(j, 2, *nrowYdata)]; // y value
      const int t = IT[i];                      // time period
      const double mu = -alphadata[i] +
	betadata[i] * theta[s][t-theta_offset[s]];
      if (y == 1.0){
	Z[s][i] = stream.rtbnorm_combo(mu, 1.0, 0);
      }
      if (y == 0.0){
	Z[s][i] = stream.rtanorm_combo(mu, 1.0, 0);
      }
    }





    // sample item parameters (alpha, beta)
    for (int i=0; i<*nitems; ++i){
      const int nsubj_i = IS[i].size();
      double sumthetasq = 0.0;
      double sumtheta = 0.0;
      double sumz = 0.0;
      double sumthetaz = 0.0;
      const int t = IT[i];

      for (int ss=0; ss<nsubj_i; ++ss){
	const int s = IS[i][ss];
	const double theta_st = theta[s][t - theta_offset[s]];
	sumthetasq += std::pow(theta_st, 2.0);
	sumtheta += -1*theta_st;
	sumz += -1*Z[s][i];
	sumthetaz += Z[s][i] * theta_st;
      } // end ss loop

      const double sumthetaallsq = std::pow(sumtheta, 2.0);
      const double invdet = 1.0 / ((nsubj_i + *A0) *
				   (sumthetasq + *B0) -
				   sumthetaallsq);
      const double v11star = invdet * (sumthetasq + *B0);
      const double v12star = invdet * (-1.0 * sumtheta);
      const double v22star = invdet * (nsubj_i + *B0);
      const double s1star = std::sqrt(v11star);
      const double s2star = std::sqrt(v22star);
      const double rho = v12star / (s1star * s2star);
      const double holder1 = sumz + A0a0;
      const double holder2 = sumthetaz + B0b0;
      const double m1star = v11star * holder1 + v12star * holder2;
      const double m2star = v12star * holder1 + v22star * holder2;
      // (alpha[i], beta[i]) ~
      //     N((m1star, m2star), c(v11star, v12star, v22star) )
      alphadata[i] = stream.rnorm(m1star, s1star);
      const double cond_mean = m2star - m1star * (v12star / v11star) +
	alphadata[i] * (v12star / v11star);
      const double cond_sd = std::sqrt(v22star * ( 1 - std::pow(rho, 2.0)));
      betadata[i] = stream.rnorm(cond_mean, cond_sd);
    } // end i loop





    // sample subject parameters (theta, tau2)
    for (int s=0; s<*nsubj; ++s){
      const int ntime_s = theta[s].size();
      if (thetaeqdata[s] != -999){
	for (int t=0; t<ntime_s; ++t){
	  theta[s][t] = thetaeqdata[s];
	}
      }
      else{

	// time period 0
	vector<double> beta_s0;
	beta_s0.reserve(STI[s][0].size());
	vector<double> zalpha_s0;
	zalpha_s0.reserve(STI[s][0].size());
	for (unsigned int ii=0; ii<STI[s][0].size(); ++ii){
	  const int i = STI[s][0][ii];
	  beta_s0.push_back(betadata[i]);
	  zalpha_s0.push_back(Z[s][i] + alphadata[i]);
	}
	vector<double> a;
	a.reserve(ntime_s);
	a.push_back(e0[s]);
	vector<double> R;
	R.reserve(ntime_s);
	R.push_back(E0inv[s] + tau2data[s]);
	vector <double> f_0 = beta_s0;
	for (unsigned int i=0; i<f_0.size(); ++i){
	  f_0[i] = f_0[i] * a[0];
	}
	Matrix <> Q_0(beta_s0.size(), beta_s0.size());
	for (unsigned int i=0; i<beta_s0.size(); ++i){
	  for (unsigned int j=i; j<beta_s0.size(); ++j){
	    if (i!=j){
	      Q_0(i,j) = Q_0(j,i) = beta_s0[i] * beta_s0[j] * R[0];
	    }
	    else{
	      Q_0(i,j) = beta_s0[i] * beta_s0[j] * R[0] + 1.0;
	    }
	  }
	}
	vector <double> e_0 = zalpha_s0;
	for (unsigned int i=0; i<e_0.size(); ++i){
	  e_0[i] = e_0[i] - f_0[i];
	}
	const Matrix <> Q_0_inv = invpd(Q_0);

	vector <double> A_0;
	A_0.reserve(beta_s0.size());
	for (unsigned int i=0; i<beta_s0.size(); ++i){
	  double sumholder = 0.0;
	  for (unsigned int j=0; j<beta_s0.size(); ++j){
	    sumholder += beta_s0[j] * Q_0_inv(j,i);
	  }
	  A_0.push_back(sumholder * R[0]);
	}
	vector<double> m;
	m.reserve(ntime_s);
	double mhold = a[0];
	for (unsigned int i=0; i<A_0.size(); ++i){
	  mhold += A_0[i] * e_0[i];
	}
	m.push_back(mhold);


	vector<double> C;
	C.reserve(ntime_s);
	double Chold = 0.0;

	for (unsigned int i=0; i<A_0.size(); ++i){
	  double hold2 = 0.0;
	  for (unsigned int j=0; j<A_0.size(); ++j){
	    //cout << "i = " << i << "   j = " << j << endl;
	    hold2 += Q_0(i,j) * A_0[j];
	  }

	  Chold += hold2 * A_0[i];
	}
	C.push_back(R[0] - Chold);



	// THE FORWARD RECURSION
	// time periods 1 to T
	for (int tt=1; tt<ntime_s; ++tt){
	  if (STI[s][tt].size() == 0){
	    a.push_back(m[tt-1]);
	    R.push_back(C[tt-1] + tau2data[s]);
	    m.push_back(a[tt]);
	    C.push_back(R[tt]);
	  }
	  else{
	    vector<double> beta_s;
	    beta_s.reserve(STI[s][tt].size());
	    vector<double> zalpha_s;
	    zalpha_s.reserve(STI[s][tt].size());
	    for (unsigned int ii=0; ii<STI[s][tt].size(); ++ii){
	      const int i = STI[s][tt][ii];
	      beta_s.push_back(betadata[i]);
	      zalpha_s.push_back(Z[s][i] + alphadata[i]);
	    }
	    // a
	    a.push_back(m[tt-1]);
	    // R
	    R.push_back(C[tt-1] + tau2data[s]);
	    vector <double> f = beta_s;
	    for (unsigned int i=0; i<f.size(); ++i){
	      f[i] = f[i] * a[tt];
	    }
	    Matrix <> Q(beta_s.size(), beta_s.size());
	    for (unsigned int i=0; i<beta_s.size(); ++i){
	      for (unsigned int j=i; j<beta_s.size(); ++j){
		if (i!=j){
		  Q(i,j) = Q(j,i) = beta_s[i] * beta_s[j] * R[tt];
		}
		else{
		  Q(i,j) = beta_s[i] * beta_s[j] * R[tt] + 1.0;
		}
	      }
	    }
	    vector <double> e = zalpha_s;
	    for (unsigned int i=0; i<e.size(); ++i){
	      e[i] = e[i] - f[i];
	    }
	    const Matrix <> Q_inv = invpd(Q);

	    vector <double> A;
	    A.reserve(beta_s.size());
	    for (unsigned int i=0; i<beta_s.size(); ++i){
	      double sumholder = 0.0;
	      for (unsigned int j=0; j<beta_s.size(); ++j){
		sumholder += beta_s[j] * Q_inv(j,i);
	      }
	      A.push_back(sumholder * R[tt]);
	    }
	    mhold = a[tt];
	    for (unsigned int i=0; i<A.size(); ++i){
	      mhold += A[i] * e[i];
	    }
	    m.push_back(mhold);

	    Chold = 0.0;
	    for (unsigned int i=0; i<A.size(); ++i){
	      double hold2 = 0.0;
	      for (unsigned int j=0; j<A.size(); ++j){
		hold2 += Q(i,j) * A[j];
	      }
	      Chold += hold2 * A[i];
	    }
	    C.push_back(R[tt] - Chold);
	  }
	} // end tt loop

	if (thetaineqdata[s] == 0){
	theta[s][ntime_s-1] = stream.rnorm(m[ntime_s-1],
					   std::sqrt(C[ntime_s-1]));
	}
	else if (thetaineqdata[s] > 0){
	  theta[s][ntime_s-1] = stream.rtbnorm_combo(m[ntime_s-1],
						     C[ntime_s-1], 0.0);
	}
	else{
	  theta[s][ntime_s-1] = stream.rtanorm_combo(m[ntime_s-1],
						     C[ntime_s-1], 0.0);
	}




	for (int tt=(ntime_s-2); tt>=0; --tt){
	  const double B = C[tt] * (1.0 / R[tt+1]);
	  const double h = m[tt] + B * (theta[s][tt+1] - a[tt+1]);
	  const double H = C[tt] - B * R[tt+1] * B;
	  if (thetaineqdata[s] == 0){
	    theta[s][tt] = stream.rnorm(h, std::sqrt(H));
	  }
	  else if (thetaineqdata[s] > 0){
	    theta[s][tt] = stream.rtbnorm_combo(h, H, 0.0);
	  }
	  else {
	    theta[s][tt] = stream.rtanorm_combo(h, H, 0.0);
	  }
	} // end backwards tt loop

      } // end theteqdataa[s] == -999



      // sample smoothing parameters (tau2)
      if (c0[s] > 0 && d0[s] > 0){
	double SSE = 0.0;
	for (int t=1; t<ntime_s; ++t){
	  SSE += std::pow((theta[s][t] - theta[s][t-1]), 2.0);
	}
	const double param1 = (c0[s] + ntime_s - 1) * 0.5;
	const double param2 = (d0[s] + SSE) * 0.5;
	tau2data[s] = stream.rigamma(param1, param2);
      }




    } // end s loop (end of theta tau2 loop)













    // store draws
    if (iter >= *burnin && (iter % *thin == 0)) {
      if (*storeability == 1){
	int internalcount = 0;
	for (int s=0; s<*nsubj; ++s){
	  for (unsigned int t=0; t<theta[s].size(); ++t){
	    thetadraws[M(sampcount, internalcount, *nrowthetadraws)] =
	      theta[s][t];
	    ++internalcount;
	  }
	}
      }
      if (*storeitem == 1){
	for (int i=0; i<*nitems; ++i){
	  alphadraws[M(sampcount, i, *nrowalphadraws)] = alphadata[i];
	  betadraws[M(sampcount, i, *nrowbetadraws)] = betadata[i];
	}
      }
      for (int s=0; s<*nsubj; ++s){
	tau2draws[M(sampcount, s, *nrowtau2draws)] = tau2data[s];
      }
      ++sampcount;
    }


    // print output to stdout
    // print output to stdout
    if(*verbose > 0 && iter % *verbose == 0){
      Rprintf("\n\nMCMCdynamicIRT1d iteration %i of %i \n", (iter+1), tot_iter);
      /*
      for (int t=0; t<theta[1].size(); ++t){
	Rprintf("%f ", theta[1][t]);
      }
      Rprintf("\n");
      for (int t=0; t<theta[7].size(); ++t){
	Rprintf("%f ", theta[7][t]);
      }
      Rprintf("\n");
      for (int t=0; t<theta[5].size(); ++t){
	Rprintf("%f ", theta[5][t]);
      }
      Rprintf("\n");
      for (int t=0; t<theta[8].size(); ++t){
	Rprintf("%f ", theta[8][t]);
      }
      Rprintf("\n");
      for (int t=0; t<theta[2].size(); ++t){
	Rprintf("%f ", theta[2][t]);
      }
      Rprintf("\n");
      for (int t=0; t<theta[4].size(); ++t){
	Rprintf("%f ", theta[4][t]);
      }
      Rprintf("\n");
      for (int t=0; t<theta[0].size(); ++t){
	Rprintf("%f ", theta[0][t]);
      }
      Rprintf("\n");
      for (int t=0; t<theta[3].size(); ++t){
	Rprintf("%f ", theta[3][t]);
      }
      Rprintf("\n");
      for (int t=0; t<theta[6].size(); ++t){
	Rprintf("%f ", theta[6][t]);
      }
      Rprintf("\n");
      Rprintf("\n\n");
      Rprintf("alpha[0] = %f \n", alphadata[0]);
      Rprintf("beta[0] = %f \n", betadata[0]);
      Rprintf("alpha[0]/beta[0] = %f \n", alphadata[0]/betadata[0]);
      Rprintf("tau2[1] = %f \n", tau2data[1]);
      Rprintf("tau2[7] = %f \n", tau2data[7]);
      Rprintf("tau2[5] = %f \n", tau2data[5]);
      Rprintf("tau2[8] = %f \n", tau2data[8]);
      Rprintf("tau2[2] = %f \n", tau2data[2]);
      Rprintf("tau2[4] = %f \n", tau2data[4]);
      Rprintf("tau2[0] = %f \n", tau2data[0]);
      Rprintf("tau2[3] = %f \n", tau2data[3]);
      Rprintf("tau2[6] = %f \n", tau2data[6]);
      Rprintf("\n\n");
      */
    }

    R_CheckUserInterrupt(); // allow user interrupts



  } // end iter loop

  for (int s=0; s<*nsubj; ++s)
    delete[] Z[s];
  delete[] Z;

} // end MCMCdynamicIRT1d_impl function








extern "C" {

  void cMCMCdynamicIRT1d(double* thetadraws, const int* nrowthetadraws,
			const int* ncolthetadraws,
			double* alphadraws, const int* nrowalphadraws,
			const int* ncolalphadraws,
			double* betadraws, const int* nrowbetadraws,
			const int* ncolbetadraws,
			double* tau2draws, const int* nrowtau2draws,
			const int* ncoltau2draws,
			const int* nsubj, const int* nitems,
			const int* ntime,
 			const int* Ydata, const int* nrowYdata,
			const int* ncolYdata,
			const int* ITdata, const int* lengthITdata,
			const int* burnin, const int* mcmc, const int* thin,
			const int* uselecuyer, const int* seedarray,
			const int* lecuyerstream,
			const int* verbose,
			const double* thetastartdata,
			const int* lengththetastart,
			const int* thetainfodata,
			const int* nrowthetainfo,
			const int* ncolthetainfo,
			double* alphastartdata,
			const int* lengthalphastart,
			double* betastartdata,
			const int* lengthbetastart,
			double* tau2startdata,
			const int* lengthtau2start,
			const double* c0,
			const int* lengthc0,
			const double* d0,
			const int* lengthd0,
			const double* a0,
			const double* A0,
			const double* b0,
			const double* B0,
			const double* e0,
			const double* E0inv,
			const double* thetaeqdata,
			const int* nrowthetaeq,
			const int* ncolthetaeq,
			const double* thetaineqdata,
			const int* nrowthetaineq,
			const int* ncolthetaineq,
			const int* storeitem,
			const int* storeability
			){


    MCMCPACK_PASSRNG2MODEL(MCMCdynamicIRT1d_impl, thetadraws, nrowthetadraws,
			   ncolthetadraws, alphadraws, nrowalphadraws,
			   ncolalphadraws, betadraws, nrowbetadraws,
			   ncolbetadraws, tau2draws, nrowtau2draws,
			   ncoltau2draws, nsubj, nitems, ntime,
			   Ydata, nrowYdata, ncolYdata,
			   ITdata,  lengthITdata, burnin,  mcmc, thin,
			   verbose, thetastartdata, lengththetastart,
			   thetainfodata, nrowthetainfo, ncolthetainfo,
			   alphastartdata, lengthalphastart, betastartdata,
			   lengthbetastart, tau2startdata, lengthtau2start,
			   c0, lengthc0, d0, lengthd0,
			   a0, A0, b0, B0, e0, E0inv, thetaeqdata, nrowthetaeq,
			   ncolthetaeq, thetaineqdata, nrowthetaineq,
			   ncolthetaineq, storeitem, storeability );


  } // end cMCMCdynamicIRT1d

} // end extern "C"


#endif
