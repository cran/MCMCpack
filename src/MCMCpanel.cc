// MCMCpanel.cc is C++ code to implement a general linear panel
// model.  It is called from R.
// 
// ADM 10/10/2002
// updated to new Scythe spec and new rngs ADM 7/28/2004

#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

extern "C"{

using namespace SCYTHE;
using namespace std;

// simulate from posterior density and return a Gibbs by parameters matrix 
// of the posterior density sample
void
panelpost (double* sample, const int* samrow, const int* samcol,
         const double* obs, const int* obsrow, const int* obscol,         
         const double* Y, const int* Yrow, const int* Ycol,
         const double* X, const int* Xrow, const int* Xcol,
         const double* W, const int* Wrow, const int* Wcol,                
         const int* burnin, const int* gibbs,  const int* thin,
         const int *lecuyer, const int *seedarray,
         const int *lecuyerstream,  const int* verbose,
         const double* bstart, const int* bstartrow, const int* bstartcol,
         const double* sigma2start, const double* Dstart,
         const int* Dstartrow, const int* Dstartcol,
         const double* b0, const int* b0row, const int* b0col,
         const double* B0, const int* B0row, const int* B0col,
         const double* nu0, const double* delta0, const int* eta0,
         const double* R0, const int* R0row, const int* R0col,
         const int* n, const int* k, const int* p, const int* q      
         )
   {
            
   // put together matrices
   Matrix<double> Msample(samcol[0], samrow[0], sample);
   Msample = t(Msample);
   Matrix<double> Mobs(obscol[0], obsrow[0], obs);
   Mobs = t(Mobs);
   Matrix<double> MX(Xcol[0], Xrow[0], X);
   MX = t(MX);
   Matrix<double> MY(Ycol[0], Yrow[0], Y);
   MY = t(MY);
   Matrix<double> MW(Wcol[0], Wrow[0], W);
   MW = t(MW);
   Matrix<double> Mbetastart(bstartcol[0],bstartrow[0], bstart);
   Mbetastart = t(Mbetastart);
   Matrix<double> MDstart(Dstartcol[0],Dstartrow[0], Dstart);
   MDstart = t(MDstart); 
   Matrix<double> Mb0(b0col[0], b0row[0], b0);
   Mb0 = t(Mb0);
   Matrix<double> MB0(B0col[0], B0row[0], B0);
   MB0 = t(MB0);
   Matrix<double> MR0(R0col[0], R0row[0], R0);  
   MR0 = t(MR0);
   
   // redefine constants
   int Mn = n[0];
   int Mk = k[0];
   int Mp = p[0];
   int Mq = q[0];
   int Mgibbs = gibbs[0];
   int Mburnin = burnin[0];
   int Mthin = thin[0];
   int Mtotiter = Mburnin + Mgibbs;

     // initialize rng stream
     rng *stream = MCMCpack_get_rng(*lecuyer, seedarray, *lecuyerstream);

   // create arrays of matrices for data access
   // Matrix<double> Yarr[Mn];
   Matrix<double>* Yarr = new Matrix<double>[Mn];
   for(int i = 0; i < Mn; ++i) {
      Yarr[i] = Matrix<double>(Mk,1);
      int start = i * Mk;
      for(int j = 0; j < Mk; ++j) {
         Yarr[i](j,0) = MY(start + j,0);
         }
      }   
      
   //Matrix<double> Xarr[Mn];
   Matrix<double>* Xarr = new Matrix<double>[Mn];
   for(int i = 0; i < Mn; ++i) {
      Xarr[i] = Matrix<double>(Mk, Mp, 0.0);
      int start = i * Mk;
      for(int l=0; l<Mk; ++l) {
         for(int m=0; m<Mp; ++m) {
            Xarr[i](l,m) = MX(start + l, m);  
            }
         }  
      }
 
   //Matrix<double> Warr[Mn];
   Matrix<double>* Warr = new Matrix<double>[Mn];
   for(int i = 0; i < Mn; ++i) {
      Warr[i] = Matrix<double>(Mk, Mq, 0.0);
      int start = i * Mk;
      for(int l=0; l<Mk; ++l) {
         for(int m=0; m<Mq; ++m) {
            Warr[i](l,m) = MW(start + l, m);  
            }
         }  
      }      

   // holding matrices
   Matrix<double> beta_holder(Mgibbs/Mthin,Mp);
   Matrix<double> D_holder(Mgibbs/Mthin,Mq*Mq);
   Matrix<double> sigma2_holder(Mgibbs/Mthin, 1);
  
   // starting values for sampler
   Matrix<double> beta_sim = Mbetastart;
   double sigma2_sim = sigma2start[0];   
   Matrix<double> D_sim = MDstart; 

   // gibbs loop
   int count = 0;
   for(int g = 0; g < Mtotiter; ++g) {
   
      // simulate \beta | Y, \sigma^2, D
      Matrix<double> beta_var_sum(Mp,Mp);
      Matrix<double> beta_mean_sum(Mp,1);
      for(int i = 0; i < Mn; ++i) {
         Matrix<double> invVi = invpd(sigma2_sim * eye<double>(Mk) + Warr[i] * D_sim * t(Warr[i]));
         beta_var_sum = beta_var_sum + t(Xarr[i]) * invVi * Xarr[i];
         beta_mean_sum = beta_mean_sum + t(Xarr[i]) * invVi * Yarr[i];
         }
      Matrix<double> beta_sim_var = invpd(MB0 + beta_var_sum);
      Matrix<double> beta_sim_mean = beta_sim_var * (MB0 * Mb0 + beta_mean_sum);
      Matrix<double> beta_sim = beta_sim_mean + cholesky(beta_sim_var) *
         stream->rnorm(Mp,1);
  
      // simulate \b_i | Y, \beta, \sigma^2, D 
      // note: these are not stored because of their large number
      Matrix<double> bi(Mn,Mq);
      for(int i = 0; i < Mn; ++i) {
         Matrix<double> b_sim_var = invpd(invpd(D_sim) + (1/sigma2_sim) * t(Warr[i]) * Warr[i]);
         Matrix<double> b_sim_mean = b_sim_var * (1/sigma2_sim * t(Warr[i]) * (Yarr[i] - Xarr[i] * beta_sim));
         Matrix<double> bi_sim = b_sim_mean + cholesky(b_sim_var) *
           stream->rnorm(Mq,1,0,1);
         for(int w = 0; w < Mq; ++w) {
            bi(i,w) = bi_sim(w,0);
         }
      }
   
      // simulate D^-1 | Y, \beta, \b_i, \sigma^2
      Matrix<double> SSB(Mq,Mq);      
      for(int i = 0; i < Mn; ++i) {
         Matrix<double> bi_sim(Mq,1);         
         for(int w = 0; w < Mq; ++w) {
            bi_sim(w,0) = bi(i,w);       
         }       
         SSB = SSB + (bi_sim * t(bi_sim)); 
      }      
      int D_sim_dof = eta0[0] + Mn;
      Matrix<double> D_sim_scale = invpd(invpd(MR0) + SSB);   
      D_sim = invpd(stream->rwish(D_sim_dof,D_sim_scale));   

      // simulate \sigma^2 | Y, \beta, \b_i, D 
      double SSE = 0;
      for(int i=0; i<Mn; ++i) {
         Matrix<double> bi_sim(Mq,1);
         for(int w = 0; w < Mq; ++w) {
            bi_sim(w,0) = bi(i,w);
            }
         Matrix<double> e = t(Yarr[i] - Xarr[i] * beta_sim - Warr[i] * bi_sim) *
                      (Yarr[i] - Xarr[i] * beta_sim - Warr[i] * bi_sim);
         SSE = SSE + e[0];
         }
      double nu_sim = (nu0[0] + Mn * Mk)/2;
      double delta_sim = (delta0[0] + SSE)/2;
      sigma2_sim = 1/stream->rgamma(nu_sim, delta_sim);

      // save values
      if (g >= Mburnin && (g % Mthin == 0)) {
         for(int j = 0; j < Mp; ++j) {
            beta_holder(count,j) = beta_sim(j,0);
            }
         int Dcounter = 0; 
         for(int j = 0; j < Mq; ++j) {
            for(int m = 0; m < Mq; ++m) {
               D_holder(count,Dcounter) = D_sim(j,m);
               ++Dcounter;            
               }
            }
         sigma2_holder(count,0) = sigma2_sim;
         ++count;
         }
 
        // print output to stdout
       if(*verbose > 0 && g % *verbose == 0) {
         Rprintf("\n\nMCMCpanel iteration %i of %i \n",
           (g+1), Mtotiter);
         Rprintf("beta = \n");
         for (int j=0; j<Mp; ++j)
           Rprintf("%10.5f\n", beta_sim[j]);
         Rprintf("sigma2 = %10.5f\n", sigma2_sim);
       }

       R_CheckUserInterrupt(); // allow user interrupt
   }
   
   delete stream; // clean up random number stream
   
  // return posterior denisty sample to R
  Matrix<double> storeagem = cbind(beta_holder, D_holder);
  storeagem = cbind(storeagem, sigma2_holder);
  int loop = samrow[0] * samcol[0];
  for (int i=0; i<loop; ++i) {
      sample[i] = storeagem[i];
      }
    
   delete [] Xarr;
   delete [] Yarr;
   delete [] Warr;
   } // panelpost function

} // extern

