// MCMCdistn.cc contains a number of functions that interface the scythe
// random number generators, density functions, and distribution functions
// with the MCMCpack R package so that these functions can be called 
// directly from R.
//
// Andrew D. Martin
// Dept. of Political Science
// Washington University in St. Louis
// admartin@wustl.edu
//
// Kevin M. Quinn
// Dept. of Government
// Harvard University
// kevin_quinn@harvard.edu
// 
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Copyright (C) 2004 Andrew D. Martin and Kevin M. Quinn
// 

#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"


#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace SCYTHE;
using namespace std;

extern "C"{

  void rtnormFromR(double* sampledata, const int* samplesize,
		   const double* mean, const int* meanlength,
		   const double* var, const int* varlength,
		   const double* below, const int* belowlength,
		   const double* above, const int* abovelength,
		   const int *lecuyer, const int *seedarray, 
		   const int *lecuyerstream){
    
    // initialize rng stream
    rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    // individual indexing to do R-style index recycling
    int m_index = -1;
    int v_index = -1;
    int b_index = -1;
    int a_index = -1;
    
    // the loop where the random number draws take place
    for (int i=0; i<*samplesize; ++i){
     
      // do the R-style index recycling
      ++m_index;
      if (m_index==*meanlength){
	m_index = 0;
      }
      ++v_index;
      if (v_index==*varlength){
	v_index = 0;
      }
      ++b_index;
      if (b_index==*belowlength){
	b_index = 0;
      }
      ++a_index;
      if (a_index==*abovelength){
	a_index = 0;
      }
      
      sampledata[i] = stream->rtnorm(mean[m_index], var[v_index],
				     below[b_index], above[a_index]);
      
    }
    
  } // end rtnormFromR

  //This next function implements the algorithm for rnchypgeom for a
  //list of data by calling double rng::rnchypgeom() for each element.
  void rnchypgeomFromR(const int * number,
		       const double* m1, const double* n1,
		       const double* n2, const double* psi,
		       const double* delta, const int * listLengths,
		       double* result, const int* lecuyer,
		       const int* seedarray, const int* lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //Second, we create some indices for wrap-around. Note that we will
    //be assuming that the listLenths array has five elements, which it
    //should if it was called from the R code.
    int count[5]={-1,-1,-1,-1,-1};

    //Third, we use the function to calculate the return values, while
    //being sure to use R-like vector wraparound for all of our parameters.
    for(int i=0; i < *number; ++i) {

      for(int j=0; j < 5; ++j) {
	count[j] = (count[j] + 1) % listLengths[j];
      }

      result[i] = stream->rnchypgeom(m1[count[0]], n1[count[1]], n2[count[2]],
				     psi[count[3]], delta[count[4]]);
    }

    //Finally, we're done, so we just end.
  }

  //This next function fills the given list from rng::richisq() from the
  //given value of the nu (degrees of freedom).
  void richisqFromR(const int * number, const double * nu,
		    const int * listLength, double * result,
		    const int *lecuyer, const int *seedarray, 
		    const int *lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //We create an index for the wrap-around.
    int count = -1;

    //We fill our array now from our function.
    for(int i = 0; i < *number; ++i) {
      count = (count + 1)%(*listLength);

      result[i] = stream->richisq(nu[count]);
    }

    //Now, we're done, so get out of here.
  }

  //This next function fills the given list from rng::rigamma() from the
  //given value of the alpha and beta.
  void rigammaFromR(const int * number, const double * alpha,
		    const double * beta, const int * listLengths,
                    double * result, const int *lecuyer,
                    const int *seedarray, const int *lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //We create an index for the wrap-around.
    int count[2] = {-1,-1};

    //We fill our array now from our function.
    for(int i = 0; i < *number; ++i) {

      for(int j = 0; j < 2; ++j) {
	count[j] = (count[j] + 1)%(listLengths[j]);
      }

      result[i] = stream->rigamma(alpha[count[0]], beta[count[1]]);
    }

    //Now, we're done, so get out of here.
  }  

  //This next function fills the given list from rng::rbern() from the
  //given value of the p parameter.
  void rbernFromR(const int * number, const double * p,
		    const int * listLength, double * result,
		    const int *lecuyer, const int *seedarray, 
		    const int *lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //We create an index for the wrap-around.
    int count = -1;

    //We fill our array now from our function.
    for(int i = 0; i < *number; ++i) {
      count = (count + 1)%(*listLength);

      result[i] = stream->rbern(p[count]);
    }

    //Now, we're done, so get out of here.
  }

  //This next function implements the algorithm for rtnormcombo for a
  //list of data by calling double rng::rtnorm_combo() for each element.
  void rtnormcomboFromR(const int * number,
		       const double* m, const double* v,
		       const double* below, const double* above,
		       const int * listLengths,
		       double* result, const int* lecuyer,
		       const int* seedarray, const int* lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //Second, we create some indices for wrap-around. Note that we will
    //be assuming that the listLenths array has four elements, which it
    //should if it was called from the R code.
    int count[4]={-1,-1,-1,-1};

    //Third, we use the function to calculate the return values, while
    //being sure to use R-like vector wraparound for all of our parameters.
    for(int i=0; i < *number; ++i) {

      for(int j=0; j < 4; ++j) {
	count[j] = (count[j] + 1) % listLengths[j];
      }

      result[i] = stream->rtnorm_combo(m[count[0]], v[count[1]], below[count[2]],
				       above[count[3]]);
    }

    //Finally, we're done, so we just end.
  }

  //This next function implements the algorithm for rtbnormslice for a
  //list of data by calling double rng::rtbnorm_slice() for each element.
  void rtbnormsliceFromR(const int * number,
		       const double* m, const double* v,
		       const double* below, const int* iter,
		       const int * listLengths,
		       double* result, const int* lecuyer,
		       const int* seedarray, const int* lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //Second, we create some indices for wrap-around. Note that we will
    //be assuming that the listLenths array has four elements, which it
    //should if it was called from the R code.
    int count[4]={-1,-1,-1,-1};

    //Third, we use the function to calculate the return values, while
    //being sure to use R-like vector wraparound for all of our parameters.
    for(int i=0; i < *number; ++i) {

      for(int j=0; j < 4; ++j) {
	count[j] = (count[j] + 1) % listLengths[j];
      }

      result[i] = stream->rtbnorm_slice(m[count[0]], v[count[1]], below[count[2]],
				       iter[count[3]]);
    }

    //Finally, we're done, so we just end.
  }

  //This next function implements the algorithm for rtanormslice for a
  //list of data by calling double rng::rtanorm_slice() for each element.
  void rtanormsliceFromR(const int * number,
		       const double* m, const double* v,
		       const double* above, const int* iter,
		       const int * listLengths,
		       double* result, const int* lecuyer,
		       const int* seedarray, const int* lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //Second, we create some indices for wrap-around. Note that we will
    //be assuming that the listLenths array has four elements, which it
    //should if it was called from the R code.
    int count[4]={-1,-1,-1,-1};

    //Third, we use the function to calculate the return values, while
    //being sure to use R-like vector wraparound for all of our parameters.
    for(int i=0; i < *number; ++i) {

      for(int j=0; j < 4; ++j) {
	count[j] = (count[j] + 1) % listLengths[j];
      }

      result[i] = stream->rtanorm_slice(m[count[0]], v[count[1]], above[count[2]],
				       iter[count[3]]);
    }

    //Finally, we're done, so we just end.
  }

  //This next function implements the algorithm for rtbnormcombo for a
  //list of data by calling double rng::rtbnorm_combo() for each element.
  void rtbnormcomboFromR(const int * number,
		       const double* m, const double* v,
		       const double* below, const int* iter,
		       const int * listLengths,
		       double* result, const int* lecuyer,
		       const int* seedarray, const int* lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //Second, we create some indices for wrap-around. Note that we will
    //be assuming that the listLenths array has four elements, which it
    //should if it was called from the R code.
    int count[4]={-1,-1,-1,-1};

    //Third, we use the function to calculate the return values, while
    //being sure to use R-like vector wraparound for all of our parameters.
    for(int i=0; i < *number; ++i) {

      for(int j=0; j < 4; ++j) {
	count[j] = (count[j] + 1) % listLengths[j];
      }

      result[i] = stream->rtbnorm_combo(m[count[0]], v[count[1]], below[count[2]],
				       iter[count[3]]);
    }

    //Finally, we're done, so we just end.
  }

  //This next function implements the algorithm for rtanormcombo for a
  //list of data by calling double rng::rtanorm_combo() for each element.
  void rtanormcomboFromR(const int * number,
		       const double* m, const double* v,
		       const double* above, const int* iter,
		       const int * listLengths,
		       double* result, const int* lecuyer,
		       const int* seedarray, const int* lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //Second, we create some indices for wrap-around. Note that we will
    //be assuming that the listLenths array has four elements, which it
    //should if it was called from the R code.
    int count[4]={-1,-1,-1,-1};

    //Third, we use the function to calculate the return values, while
    //being sure to use R-like vector wraparound for all of our parameters.
    for(int i=0; i < *number; ++i) {

      for(int j=0; j < 4; ++j) {
	count[j] = (count[j] + 1) % listLengths[j];
      }

      result[i] = stream->rtanorm_combo(m[count[0]], v[count[1]], above[count[2]],
				       iter[count[3]]);
    }

    //Finally, we're done, so we just end.
  }

  //This next function implements the algorithm for rwish for a
  //by calling double rng::rwish() for the whole matrix.
  void rwishFromR(const int * v, const double* s,
		  const int * dimension,
		  double* result, const int* lecuyer,
		  const int* seedarray, const int* lecuyerstream) {

    //First, we initialize our rng stream, being sure to store the parameters
    //so that we can check for them the next time this is run.
    static rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);
    static int store_lecuyer = *lecuyer;
    static int store_seedarray[6] = {0,0,0,0,0,0};
    static int store_lecuyerstream = *lecuyerstream;
    if((*lecuyer!=store_lecuyer)||(*lecuyerstream!=store_lecuyerstream)||
       (store_seedarray[0]!=seedarray[0])||(store_seedarray[1]!=seedarray[1])||
       (store_seedarray[2]!=seedarray[2])||(store_seedarray[3]!=seedarray[3])||
       (store_seedarray[4]!=seedarray[4])||(store_seedarray[5]!=seedarray[5])) {
      for(int i = 0; i < 6; ++i) store_seedarray[i] = seedarray[i];
      store_lecuyer = *lecuyer;
      store_lecuyerstream = *lecuyerstream;
      delete stream;
      stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);
    }

    //Second, we run our rng::rwish() function to get the result
    Matrix<double> * myMat = new Matrix<double>(*dimension,*dimension,s);
    Matrix<double> myMatrix2 = stream->rwish(*v,*myMat); 

    for(int i = 0; i < (*dimension)*(*dimension); ++i) {
      result[i] = myMatrix2[i];
    }

    //So that we don't have any memory leaks:
    delete myMat;

    //Finally, we're done, so we just end.
  }

  //This next function implements the algorithm for rdirich for a
  //by calling double rng::rdirich() for each element.
  void rdirichFromR(const int * n, const double* alpha,
		    const int * vsize, const int * hsize,
		    double* result, const int* lecuyer,
		    const int* seedarray, const int* lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //Second, we run our rng::rdirich() function to get the result
    Matrix<double> alphaMatrix(*hsize,*vsize,alpha);
    Matrix<double> temp, tempVector;

    for(int i=0; i < *n; ++i) {
      tempVector = alphaMatrix(i%(*hsize), _);
      tempVector.resize(*vsize, 1);
      temp = stream->rdirich(tempVector);
      for(int j = 0; j < *vsize; ++j)
	result[i*(*vsize)+j] = temp[j];
    }

    //Finally, we're done, so we just end.
  }

  //This next function implements the algorithm for rmvnorm for a
  //by calling double rng::rmvnorm() for each element.
  void rmvnormFromR(const int * n, const double* mu,
		    const int * vsize, const int * hsize,
		    const double * sigma,
		    double* result, const int* lecuyer,
		    const int* seedarray, const int* lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //Second, we run our rng::rdirich() function to get the result
    Matrix<double> muMatrix(*hsize,*vsize,mu);
    Matrix<double> sigmaMatrix(*vsize,*vsize,sigma);
    Matrix<double> temp, tempVector;

    for(int i=0; i < *n; ++i) {
      tempVector = muMatrix(i%(*hsize), _);
      tempVector.resize(*vsize, 1);
      temp = stream->rmvnorm(tempVector, sigmaMatrix);
      for(int j = 0; j < *vsize; ++j)
	result[i*(*vsize)+j] = temp[j];
    }

    //Finally, we're done, so we just end.
  }

  //This next function implements the algorithm for rmvt for a
  //by calling double rng::rmvt() for each element.
  void rmvtFromR(const int * n, const double* sigma,
		 const int * size, const double * nu, const int * sizeNu,
		 double* result, const int* lecuyer,
		 const int* seedarray, const int* lecuyerstream) {

    //First, we initialize our rng stream.
    rng * stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

    //Second, we run our rng::rdirich() function to get the result
    Matrix<double> sigmaMatrix(*size,*size,sigma);
    Matrix<double> temp;

    for(int i=0; i < *n; ++i) {
      temp = stream->rmvt(sigmaMatrix, nu[(i%(*sizeNu))]);
      for(int j = 0; j < *size; ++j)
	result[i*(*size)+j] = temp[j];
    }

    //Finally, we're done, so we just end.
  }

  // other *FromR functions should be added here


  
} // end extern "C"
