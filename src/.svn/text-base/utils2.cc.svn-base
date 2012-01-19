#include "la.h"
#include "matrix.h" 
#include "stat.h"
#include "ide.h"
#include "mersenne.h"
#include "distributions.h"
#include "rng.h"
#include <R_ext/PrtUtil.h>
//#include <Rmath.h>
//#include <R_ext/RS.h>     /* for Calloc/Free */
//#include <R_ext/Applic.h> /* for dgemm */

using namespace scythe;
using namespace std;

mersenne myrng;


extern "C"{ 
  void demean(const double* Adata, const int* ncol, 
	      const int* nrow, double* out)
  {
    Matrix<double, Col> A(*nrow, *ncol, Adata);
    Matrix<double> Am = meanc(A);
    
    for(unsigned int i = 0; i< A.rows(); ++i)
      for(unsigned int j = 0; j< A.cols(); ++j)
	A(i,j) = A(i,j)-Am(j);

    for(unsigned int i = 0; i< (A.rows()*A.cols()); ++i)
      out[i] = A(i);
  }
}

extern "C"{
  void wildbootr(const double* Xdata, const double* Ydata, 
		 const double* resid,
		 const double* beta, 
		 const double* factor, 
		 const int* wres, const int* lwres,
		 const int* nrow, const int *ncol,
		 const int* nrep, const int* clusstart, 
		 const int* clussize, const int* nc, 
		 const int* wtyp,
		 double* out1, double* out2)
  {

    /* wtyp
       = 1,    radamacher
       = 2,    mtp
       = 3,    mn1
       = 4,    mn2
     */
    
    void robcovf(const int n, const int kv, const int ncl, const int* start, 
		 const int* len, double* u, double* s, double* v, 
		 double* w);
    // Containers
    const Matrix<double> X(*nrow, *ncol, Xdata);
    const Matrix<double> Y(*nrow,     1, Ydata);
    const Matrix<double> XX = invpd(crossprod(X));
    const Matrix<double> bhat(*ncol, 1, beta);
    const Matrix<double> uhat(*nrow, 1, resid);

    // Working Matrices
    Matrix<double> V(*ncol, *ncol);
    Matrix<double> WW(*ncol, *ncol);      // Store WW unrestricted
    Matrix<double> bout(*nrep, *ncol);    // Store betahat_wild
    Matrix<double> SE(*nrep, *ncol);
    Matrix<double> ustar = uhat;          // Store uwild
    Matrix<double> Xs    = X;
    Matrix<double> Yhat  = X*bhat;
    
    // Constant
    const int n        = X.rows();
    const int k        = X.cols();
    const int r        = bout.rows();
    const int ncl      = *nc;
    const int start    = *clusstart;
    const int len      = *clussize;
    const int wr       = *wres;
    const int lwr      = *lwres;
    const double fac   = *factor;
    // Wild bootstrap weights type
    const int wtype = *wtyp;
    // Wild weights constants
    const double tp1    = -(sqrt(5.)-1.)/2.;
    const double tp2    =  (sqrt(5.)+1.)/2.;
    const double tpp    =  (sqrt(5.)+1.)/(2.*sqrt(5.));
    const double delta1 =  sqrt(3./4.+sqrt(17.)/12.);
    const double delta2 =  sqrt(3./4.-sqrt(17.)/12.);
    // weights and arrays to be passed to robcovf()
    double ti, u[n*k], s[k], v[k*k], w[k*k];

    for(unsigned int i = 0; i< r; ++i)
      {
	for(unsigned int j = 0; j<ncl; ++j)
	  {
	    int start = clusstart[j]-1;
	    int end   = clusstart[j]+clussize[j]-1;
	    if(wtype == 1) {
	      ti = myrng();
	      if(ti<0.5) ti = -1; else ti = 1;
	    }
	    else if(wtype == 2) {
	      ti = myrng();
	      if(ti<tpp) ti = tp1; else ti = tp2;
	    }
	    else if(wtype == 3) {
	      ti = myrng.rnorm(0,1);
	      ti = ti/sqrt(2.0)+(pow(ti,2.0)-1.)/2.;
	    }
	    else if(wtype == 4) {
	      ti = myrng.rnorm(0,1);
	      ti = (delta1+ti/sqrt(2.))*(delta2+ti/sqrt(2.))-delta1*delta2;
	    }
	    for(unsigned int jj = start; jj<end; ++jj) 
	      ustar(jj) = uhat(jj)*ti;
	  }

	bout(i,_) = bhat + XX*t(X)*(ustar);
	ustar = (Yhat+ustar)-X*t(bout(i,_));
	for(unsigned int j = 0; j<n; ++j)
	  for(unsigned int i = 0; i<k; ++i)
	    Xs(j,i) = (X(j,i)*ustar(j))*fac;
	
	if(ncl==1)
	  V = XX*crossprod(Xs)*XX;
	else {
	  for(unsigned int j = 0; j<(n*k); ++j)
	    u[j] = Xs(j);
	  // Call robcovf 
	  robcovf(n, k, ncl, clusstart, clussize, u, s, v, w);
	  
	  for(unsigned int uu = 0; uu<(k*k); ++uu)
	    WW(uu) = w[uu];
	  
	  V = XX*WW*XX;
	}

	if(lwr>1)
	  // Which means we are doing the unconstrained
	  for(unsigned int ii = 0; ii<k; ++ii)
	    SE(i, ii) = sqrt(V(ii,ii)); 
	else {
	  out1[i] = bout(i, (wr-1));
	  out2[i] = sqrt(V((wr-1), (wr-1)));
	}
      }
    
    if(lwr>1) {
      for(unsigned int i = 0; i<(r*k); ++i) 
	{
	  out1[i] = bout(i);
	  out2[i] = SE(i);
	}
    }
  }

  void robcovf(const int n, const int kv, const int ncl,  int* start, 
	       int* len, double* u, double* s, double* v, double* w)
  {
    /* System generated locals */
    int i2;
    
    /* Local variables */
    static int i, j, k;
    
    for (i = 0; i < kv; ++i) {
      for (j = 0; j < kv; ++j) {
	w[i + j * kv] = 0.;
      }
    }
    
    for (k = 0; k < ncl; ++k) {
      for (i = 0; i < kv; ++i) {
	s[i] = 0.;
	for (j = 0; j < kv; ++j) {
	  v[i + j * kv] = 0.;
	}
      }
      i2 = start[k] + len[k] - 1;
      for (i = (start[k]-1); i < i2; ++i) {
	for (j = 0; j < kv; ++j) {
	  s[j] += u[i + j * n];
	}
      }
      
      //     for (j = 0; j < kv; ++j) {
      // 	Rprintf("s[%d] = %f \n",  j, s[j]);
      //       }

      for (i = 0; i < kv; ++i) {
	for (j = 0; j < kv; ++j) {
	  v[i + j * kv] += s[i] * s[j];
	}
      }
      
//     for (i = 0; i < kv; ++i) {
//       for (j = 0; j < kv; ++j) {
// 	Rprintf("v[%d,%d] = %f \n", i, j, v[i + j * kv]);
//       }
//     }
      
      for (i = 0; i < kv; ++i) {
	for (j = 0; j < kv; ++j) {
	  w[i + j * kv] += v[i + j * kv];
	  //	Rprintf("%f %f \n", v[i + j * kv], w[i + j * kv]);//
	}
      }
    }
  }
}


