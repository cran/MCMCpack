/* 
 * Scythe Statistical Library
 * Copyright (C) 2000-2002 Andrew D. Martin and Kevin M. Quinn;
 * 2002-2004 Andrew D. Martin, Kevin M. Quinn, and Daniel
 * Pemstein.  All Rights Reserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * under the terms of the GNU General Public License as published by
 * Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.  See the text files COPYING
 * and LICENSE, distributed with this source code, for further
 * information.
 * --------------------------------------------------------------------
 * scythestat/rng.h
 *
 * Provides the class definition for the rng class.  This abstract
 * class forms the foundation of random number generation in Scythe.
 * Specific random number generators should extend this class and
 * implement the virtual void function runif(); this function should
 * take no arguments and return uniformly distributed random numbers
 * on the interval (0, 1).  The rng class provides no interface for
 * seed-setting or intitialization, allowing for maximal flexibility
 * in underlying implementation.  This class does provide
 * implementations of functions that return random numbers from a wide
 * variety of commonly (and not-so-commonly) used distributions by
 * manipulating the uniform deviates returned by runif().
 *
 * The code for many of the RNGs defined in this file and implemented
 * in rng.cc is based on that in the R project, version 1.6.0-1.7.1.
 * This code is available under the terms of the GNU GPL.  Original
 * copyright:
 * 
 * Copyright (C) 1998      Ross Ihaka
 * Copyright (C) 2000-2002 The R Development Core Team
 * Copyright (C) 2003      The R Foundation
 *
 */

#ifndef SCYTHE_RNG_H
#define SCYTHE_RNG_H

#ifdef SCYTHE_COMPILE_DIRECT
#include "matrix.h"
#else
#include "scythestat/matrix.h"
#endif

namespace SCYTHE {

	class rng
	{
		public:

			/* Default (and only) constructor. */
			rng ();

			/* Destructor */
			virtual ~rng();

			/* Returns random uniform numbers on (0, 1).  This function must
			 * be implemented by extending classes */
			virtual double runif () = 0;

			virtual Matrix<double> runif (const int &, const int &);
			
			/**** Random deviates from various distributions */
			
			/* Beta distribution */
			double rbeta (const double &, const double &);

			Matrix<double> rbeta (const int &, const int &,
														const double &, const double &);
			
			/* Non-central hypergeometric distribution */
			double rnchypgeom(const double& m1, const double& n1, 
												const double& n2, const double& psi, 
												const double& delta=1e-14); 
			
			Matrix<double> rnchypgeom(const int &, const int &,
																const double&, const double&, 
																const double&, const double&, 
																const double& delta=1e-14); 

			/* Binomial distribution */
			int rbinom (const int &, const double &);
			
			Matrix<double> rbinom ( const int &, const int &,
															const int &, const double &);

			/* Chi^2 distribution */
			double rchisq (const double &);

			Matrix<double> rchisq ( const int &, const int &,
															const double &);

			/* Exponential distribution */
			double rexp (const double &);

			Matrix<double> rexp ( const int &rows, const int &cols,
														const double &);

			/* f distribution */
			double rf(const double &, const double &);

			Matrix<double> rf(const int &, const int &,
												const double &, const double &);

			/* Gamma distribution */
			double rgamma (const double &, const double &);

			Matrix<double> rgamma ( const int &, const int &,
														const double &, const double &);

			/* Logistic distribution */
			double rlogis (const double &, const double &);

			Matrix<double> rlogis ( const int &, const int &,
															const double &, const double &);
		
			/* Log Normal distribution */
			double rlnorm(const double &logmean = 0.0,
										const double &logsd = 1.0);

			Matrix<double> rlnorm(const int &, const int &,
														const double &logmean = 0.0,
														const double &logsd = 1.0);
			
			/* Negative Binomial distribution */
			double rnbinom(const double & , const double &);

			Matrix<double> rnbinom(const int &, const int &, const double &,
														const double &);

			/* Normal distribution */
			double rnorm (const double &mu=0.0, const double &sigma=1.0);

			Matrix<double> rnorm (const int &rows, const int &cols,
														const double &mu=0.0, 
														const double &sigma=1.0);
			
			/* Poisson distribution */
			int rpois(const double &);

			Matrix<int> rpois (const int &, const int &, const double &);

			/* Student's t distribution */
			double rt (const double &, const double &, const double &);

			Matrix<double> rt ( const int &, const int &,
													const double &, const double &,
													const double &);

			/* Weibull distribution */
			double rweibull (const double &, const double &scale = 1.0);

			Matrix<double> rweibull(const int &, const int &, const double &,
															const double &scale = 1.0);
			
			/* Inverse Chi^2 distribution */
			double richisq (const double &);

			Matrix<double> richisq (const int &, const int &,
															const double &);

			/* Inverse Gamma distribution */
			double rigamma (const double &, const double &);

			Matrix<double> rigamma( const int &, const int &,
															const double &, const double &);

			/* Wishart (Only for Matrix) distribution */ 
			Matrix<double> rwish(const int &, const Matrix<double> &);

			/* Dirichlet distribution */
			Matrix<double> rdirich(const Matrix<double> &);
				
			/* Multivariate Normal distribution */
			Matrix<double> rmvnorm (const Matrix<double> &,
															const Matrix<double> &);

			/* Multivariate t distribution */
			Matrix<double> rmvt (const Matrix<double> &, const double &);

			/* Bernoulli distribution */
			int rbern (const double &);

			Matrix<double> rbern (const int &, const int &,
														const double &);

			/* Beta-Binomial distribution */
			int rbetabin (const int &, const double &, const double &);

			Matrix<double> rbetabin ( const int &, const int &,
																const int &, const double &,
																const double &);

			/* Truncated Normal distribution */
			
			/* Truncated Normal distribution */
			double rtnorm(const double &, const double &, const double &, 
										const double &);
			
			double rtnorm_combo(const double &, const double &,
													const double &, const double &);

			double rtbnorm_slice( const double &, const double &,
														const double &, const int &iter = 10);

			double rtanorm_slice( const double &, const double &,
														const double &, const int &iter = 10);

			double rtbnorm_combo( const double &, const double &,
														const double &, const int &iter=10);

			double rtanorm_combo( const double &, const double &,
														const double &, const int &iter=10);

		protected:

			double rgamma1 (const double &);
  
			double rnorm1();

	};

}

#endif /* SCYTHE_RNG_H */
