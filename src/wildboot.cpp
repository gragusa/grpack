#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::RNGScope Scope;
using namespace Rcpp;


extern "C" {
  SEXP wb_null2(SEXP Xr, SEXP Rr, SEXP yr, SEXP betar, SEXP resr, SEXP clusstartr,
		SEXP clussizer, SEXP reps, SEXP wbtype) {

  arma::mat X=Rcpp::as<arma::mat>(Xr);
  arma::mat XX=Rcpp::as<arma::mat>(Rr);
  arma::mat y=Rcpp::as<arma::colvec>(yr);
  arma::mat beta=Rcpp::as<arma::colvec>(betar);
  arma::mat uhat=Rcpp::as<arma::colvec>(resr);

  Rcpp::IntegerVector clusstart(clusstartr);
  Rcpp::IntegerVector clussize(clussizer);
  Rcpp::IntegerVector replications(reps);
  Rcpp::IntegerVector type(wbtype);
  
  void robcovf(const int n, const int kv, const int ncl, const int* start, 
	       const int* len, double* u, double* s, double* v, 
	       double* w);
  // int n : 
  // const int kv:
  // const int ncl: number of cluster  
  // const int* start: clus.start
  // const int* len: clus.size
  // double* u:
  // double* s:
  // double* v:
  // double* w

  int obs = X.n_rows;
  int k   = X.n_cols;
  int ncl = clussize.size();
  
  int clstart[ncl];
  int clsize[ncl];

  for(unsigned int j = 0; j<ncl; ++j)
    {
      clstart[j] = (int)clusstart(j);
      clsize[j] =  (int)clussize(j);
    }

  const int *c1 = clstart;
  const int *c2 = clsize;

  double u[obs*k], s[k], v[k*k], w[k*k];

  // Wild weights constants
  const double tp1    = -(sqrt(5.)-1.)/2.;
  const double tp2    =  (sqrt(5.)+1.)/2.;
  const double tpp    =  (sqrt(5.)+1.)/(2.*sqrt(5.));
  const double delta1 =  sqrt(3./4.+sqrt(17.)/12.);
  const double delta2 =  sqrt(3./4.-sqrt(17.)/12.);
  // scaling
  const double adj1 = ((double)(obs-1)/(double)obs)*(((double)ncl/(double)(ncl-1)));
  const double adj3 = ((double)ncl/(double)(ncl-1));
  double ti;

  
  arma::colvec ustar(obs);
  arma::colvec Yhat=X*beta;
  arma::mat bout(replications[0], k);
  arma::mat seout_HC1(replications[0], k);
  arma::mat seout_HC2(replications[0], k);
  arma::mat seout_HC3(replications[0], k);
  arma::mat score(obs, k);
  arma::mat WW(k,k);
  arma::mat V(k,k);
  arma::mat Xp=strans(X);
  int idr, idy;
  for(int i = 0; i<replications[0]; ++i)
    {
      for(unsigned int j = 0; j<ncl; ++j)
	{
	  idr = clusstart[j]-1;
	  idy = (clussize[j]-2)+clusstart[j];
	  if(type[0] == 1) {
	    ti = (runif(1, 0, 1))[0];
	    if(ti<0.5) ti = -1; else ti = 1;
	  }
	  else if(type[0] == 2) {
	    ti = (runif(1,0,1))[0];
	    if(ti<tpp) ti = tp1; else ti = tp2;
	  }
	  else if(type[0] == 3) {
	    ti = (rnorm(1, 0,1))[0];
	    ti = ti/sqrt(2.0)+(pow(ti,2.0)-1.)/2.;
	  }
	  else if(type[0] == 4) {
	    ti = (rnorm(1,0,1))[0];
	    ti = (delta1+ti/sqrt(2.))*(delta2+ti/sqrt(2.))-delta1*delta2;
	  }
	  ustar.rows(idr, idy) = uhat.rows(idr, idy)*ti;
	}
      
      bout.row(i) = (beta + XX*Xp*ustar).st(); 
      // Rcout << "beta\n" << std::endl;
      // Rf_PrintValue(wrap(bout));
      // Rcout << "ustar\n" << std::endl;
      // Rf_PrintValue(wrap(ustar));
      ustar = Yhat+ustar-X*(bout.row(i)).st();
      WW.zeros(); 

      // HC1 type variance
      for(unsigned int j = 0; j<obs; ++j)
	for(unsigned int i = 0; i<k; ++i)
	  score(j,i) = (X(j,i)*ustar(j));
      for(unsigned int j = 0; j<(obs*k); ++j)
	u[j] = score(j);     
      robcovf(obs, k, ncl, c1, c2, u, s, v, w);
      for(unsigned int j = 0; j<(k*k); ++j)
	WW(j)=w[j];
      V = adj1*(XX*WW*XX);
      seout_HC1.row(i) = sqrt((V.diag()).st());

      // HC2 type variance      
      for(unsigned int jj = 0; jj<ncl; ++jj) {
	idr = clusstart[jj]-1;
	idy = (clussize[jj]-2)+clusstart[jj];
	arma::mat Xi  = X.rows(idr, idy);
	arma::mat Identity=arma::eye(Xi.n_rows, Xi.n_rows);
	arma::mat Hgg = arma::chol(Identity - Xi*XX*strans(Xi));
	ustar.rows(idr, idy) = arma::solve(Hgg, ustar.rows(idr, idy));
      }
      for(unsigned int j = 0; j<obs; ++j)
	for(unsigned int i = 0; i<k; ++i)
	  score(j,i) = (X(j,i)*ustar(j));
      for(unsigned int j = 0; j<(obs*k); ++j)
	u[j] = score(j);     
      robcovf(obs, k, ncl, c1, c2, u, s, v, w);
      for(unsigned int j = 0; j<(k*k); ++j)
	WW(j)=w[j];
      V = (XX*WW*XX);
      seout_HC2.row(i) = sqrt((V.diag()).st());
      
      // HC3 Type Variance

      for(unsigned int jj = 0; jj<ncl; ++jj) {
	idr = clusstart[jj]-1;
	idy = (clussize[jj]-2)+clusstart[jj];
	arma::mat Xi  = X.rows(idr, idy);
	arma::mat Identity=arma::eye(Xi.n_rows, Xi.n_rows);
	arma::mat Hgg = Identity - Xi*XX*strans(Xi);
	ustar.rows(idr, idy) = arma::solve(Hgg, ustar.rows(idr, idy));
      }
      for(unsigned int j = 0; j<obs; ++j)
	for(unsigned int i = 0; i<k; ++i)
	  score(j,i) = (X(j,i)*ustar(j));
      for(unsigned int j = 0; j<(obs*k); ++j)
	u[j] = score(j);     
      robcovf(obs, k, ncl, c1, c2, u, s, v, w);
      for(unsigned int j = 0; j<(k*k); ++j)
	WW(j)=w[j];
      V = adj3*XX*WW*XX;      
      seout_HC3.row(i) = sqrt((V.diag()).st());
    }
  
  // return
  return Rcpp::List::create(
			    Rcpp::Named("coef_wb") = bout,
			    Rcpp::Named("sd_wb_HC1") = seout_HC1, 
			    Rcpp::Named("sd_wb_HC2") = seout_HC2,
			    Rcpp::Named("sd_wb_HC3") = seout_HC3
			    );
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

      for (i = 0; i < kv; ++i) {
	for (j = 0; j < kv; ++j) {
	  v[i + j * kv] += s[i] * s[j];
	}
      }
        
      for (i = 0; i < kv; ++i) {
	for (j = 0; j < kv; ++j) {
	  w[i + j * kv] += v[i + j * kv];
	}
      }
    }
  }
}

extern "C" SEXP resHC3(SEXP Xr, SEXP rr, SEXP Rr, SEXP clusstartr, SEXP clussizer) {
  arma::mat X=Rcpp::as<arma::mat>(Xr);
  arma::mat r=Rcpp::as<arma::colvec>(rr);
  arma::mat R=Rcpp::as<arma::mat>(Rr);
  Rcpp::IntegerVector clusstart(clusstartr);
  Rcpp::IntegerVector clussize(clussizer);
  int ncl = clussize.size();
  int obs = X.n_rows;
  arma::mat res(obs,1);
  for(unsigned int jj = 0; jj<ncl; ++jj) {
    int idr = clusstart[jj]-1;
    int idy = (clussize[jj]-2)+clusstart[jj];
    arma::mat Xi  = X.rows(idr, idy);
    arma::mat Identity=arma::eye(Xi.n_rows, Xi.n_rows);
    res.rows(idr, idy) = arma::solve(Identity - Xi*R*Xi.st(), r.rows(idr, idy));
  }
  res = sqrt((double)ncl/(double)(ncl-1))*res;
  return Rcpp::wrap(res);
}

extern "C" SEXP resHC2(SEXP Xr, SEXP rr, SEXP Rr, SEXP clusstartr, SEXP clussizer) {
  arma::mat X=Rcpp::as<arma::mat>(Xr);
  arma::mat r=Rcpp::as<arma::colvec>(rr);
  arma::mat R=Rcpp::as<arma::mat>(Rr);
  Rcpp::IntegerVector clusstart(clusstartr);
  Rcpp::IntegerVector clussize(clussizer);
  int ncl = clussize.size();
  int obs = X.n_rows;
  arma::mat res(obs,1);
  for(unsigned int jj = 0; jj<ncl; ++jj) {
    int idr = clusstart[jj]-1;
    int idy = (clussize[jj]-2)+clusstart[jj];
    arma::mat Xi  = X.rows(idr, idy);
    arma::mat Identity=arma::eye(Xi.n_rows, Xi.n_rows);
    arma::mat Hgg = arma::chol(Identity - Xi*R*Xi.st());
    res.rows(idr, idy) = arma::solve(Hgg, r.rows(idr, idy));
  }
  return Rcpp::wrap(res);
}





