// MCMCfcds.h is the header file for MCMCfcds.cc. It contains declarations 
// for a number of functions that produce draws from full conditional 
// distributions. 
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
// KQ 6/10/2004


#ifndef MCMCFCDS_H
#define MCMCFCDS_H

#include "matrix.h"
#include "smath.h"
#include <cfloat>

namespace SCYTHE {

  // linear regression with Gaussian errors beta draw (multivariate Normal 
  // prior)
  Matrix<double>
    NormNormregress_beta_draw (const Matrix<double> &, 
			       const Matrix<double> &,
			       const Matrix<double> &, 
			       const Matrix<double> &,
			       const double &,
			       rng *);
			       
  // linear regression with Gaussian errors sigma2 draw (inverse-Gamma
  // prior)
  double
    NormIGregress_sigma2_draw (const Matrix<double> &, 
			       const Matrix<double> &,
			       const Matrix<double> &, 
			       const double &,
			       const double &,
			       rng *);    
  // factor scores with N(0, F0^-1) prior
  void 
    NormNormfactanal_phi_draw (Matrix<double> &, 
			       const Matrix<double> &, 
			       const Matrix<double> &,
			       const Matrix<double> &,
			       const Matrix<double> &,
			       const int&, const int&,
			       rng *);
  
  // factor loading matrix
  void
    NormNormfactanal_Lambda_draw(Matrix<double> &, 
				 const Matrix<double> &,
				 const Matrix<double> &,
				 const Matrix<double> &,
				 const Matrix<double> &,
				 const Matrix<double> &,
				 const Matrix<double> &,
				 const Matrix<double> &,
				 const int&, const int&,
				 rng *);

  // samples the Psi matrix for a Normal theory factor model with IG 
  // prior on diag elements of Psi
  void 
    NormIGfactanal_Psi_draw(Matrix<double> &, const Matrix<double> &,
			    const Matrix<double> &, 
			    const Matrix<double> &,
			    const Matrix<double> &,
			    const Matrix<double> &,
			    const int&, const int&,
			    rng *);

  // updates for MCMCirt1d
  void irt_Z_update1(Matrix<double> &, const Matrix<int> &, 
		     const Matrix<double> &, const Matrix<double> &, 
		     rng *);
  
  void irt_eta_update1(Matrix<double> &, const Matrix<double> &,
		       const Matrix<double> &, const Matrix<double> &, 
		       const Matrix<double> &, rng *);
  
  void irt_theta_update1(Matrix<double>&, const Matrix<double>&,
			 const Matrix<double>&, const double&, 
			 const double&, const Matrix<double>&, 
			 const Matrix<double>&, rng *); 
}
#endif
