// fits Wakefield's baseline model for ecological inference via 
// data augmentation
//
// KQ 2/24/2002
// KQ 10/23/2002 [ported to Scythe0.3 and written for an R interface]
//

#include <iostream> 
#include "Scythe_Matrix.h"
#include "Scythe_Simulate.h"
#include "Scythe_Stat.h"
#include "Scythe_Math.h"
#include "Scythe_LA.h"
#include "Scythe_IDE.h"

extern "C"{
 
using namespace SCYTHE;
using namespace std;

  void baselineDA(double* sample, const int* samrow, const int* samcol,
		  const double* Rr0, const double* Rr1, const double* Rc0,
		  const double* Rc1, const int* Rntables, const int* Rburnin,
		  const int* Rmcmc, const int* Rthin, const double* Ralpha0,
		  const double* Rbeta0, const double* Ralpha1, 
		  const double* Rbeta1, const int* Rverbose, 
		  const int* Rseed){



    // load data
    // table notation is:
    // --------------------
    //   Y0  |     | r0
    // --------------------
    //   Y1  |     | r1
    // --------------------
    //   c0  | c1  | N

  
    // initialize seed (mersenne twister / use default seed unless specified)
    if(Rseed[0]==0) 
      set_mersenne_seed(5489UL);
    else 
      set_mersenne_seed(Rseed[0]);
    
    
    int ntables = *Rntables;
    int verbose = *Rverbose;

    Matrix<double> r0(ntables, 1, Rr0);
    Matrix<double> r1(ntables, 1, Rr1);
    Matrix<double> c0(ntables, 1, Rc0);
    Matrix<double> c1(ntables, 1, Rc1);
    Matrix<double> N = c0 + c1;

    
    // MCMC-related quantities
    int burnin = *Rburnin;
    int mcmc =   *Rmcmc;
    int thin =   *Rthin;
    int tot_iter = burnin + mcmc;

    // prior for p0
    // p0 ~ beta(alpha0, beta0)
    double alpha0 = *Ralpha0;
    double beta0  = *Rbeta0;
    
    // prior for p1
    // p1 ~ beta(alpha1, beta1)
    double alpha1 = *Ralpha1;
    double beta1  = *Rbeta1;

    // storage matrices
    Matrix<double> p0mat(mcmc/thin, ntables);
    Matrix<double> p1mat(mcmc/thin, ntables);
    Matrix<double> y0mat(mcmc/thin, ntables);
    Matrix<double> y1mat(mcmc/thin, ntables);
    int count = 0;
    
    // starting values
    Matrix<double> p0 = ones<double>(ntables,1)*0.5;
    Matrix<double> p1 = ones<double>(ntables,1)*0.5;
    Matrix<double> y0(ntables,1);
    Matrix<double> y1(ntables,1);

    for (int iter=0; iter<tot_iter; ++iter){

      for (int j=0; j<ntables; ++j){
	// sample y0|c1,r0,r1,p0,p1
	double psi = ( p0[j]*(1.0-p1[j]) ) / ( p1[j]*(1.0-p0[j]));
	y0[j] = rnchypgeom(c0[j], r0[j], r1[j], psi);
	y1[j] = c0[j] - y0[j];
      
	// sample (p0,p1)|y0,y1,r0,r1,c0,c1
	p0[j] = rbeta(alpha0+y0[j], beta0+(r0[j]-y0[j]));
	p1[j] = rbeta(alpha1+y1[j], beta1+(r1[j]-y1[j]));
      
	// if after burnin store samples
	if ((iter >= burnin) && ((iter%thin)==0)){
	  p0mat(count,j) = p0[j];
	  p1mat(count,j) = p1[j];
	  y0mat(count,j) = y0[j];
	  y1mat(count,j) = y1[j];
	}
      }
      if ((iter>=burnin) && ((iter%thin)==0)) ++count;
      // print output to screen
      if (verbose==1 && (iter%10000)==0){
	  cout << "MCMCbaselineEI iteration = " << iter << endl;
      }
    
    }

    // return sample
    Matrix<double> storeagem = cbind(p0mat, p1mat);
    storeagem = cbind(storeagem, y0mat);
    storeagem = cbind(storeagem, y1mat);
    int mat_size = samrow[0] * samcol[0];
    for (int i=0; i<mat_size; ++i)
      sample[i] = storeagem[i];


  }



}// extern "C"
