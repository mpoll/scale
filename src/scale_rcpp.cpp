#include <Rcpp.h>
#include <cmath>
#include <cstdlib>

#include "scale_rcpp.hh"

using namespace Rcpp;
using namespace std;

double _eachi(unsigned long j, double s, double t, double m, double xoy, double u) {
  double aum   = abs(u-m);
  double P     = -2.0*aum*j/(t-s);
  double rval  = (2.0*aum*j+(xoy-m))*exp(P*(aum*j+(xoy-m)));

  return(rval);
}

double _eadel(unsigned long n, double s, double t, double x, double y,  double m, double u) {
  unsigned short int xI = (x==m ? 1 : 0);
  unsigned short int yI = (y==m ? 1 : 0);
  unsigned short int delT = (xI | yI) + 1;

  if(t==s)
    return 1;

  if(m > max(x,y)){
     x = -x;
     y = -y;
     u = -u;
     m = -m;
  }

  double rval;
  /*  if((max(m-u,y-u) >= 0) || (max(m-x,m-y) >= 0)) {
    rval = 0;
    }*/
  if(delT == 1) {
    rval = _eagamma(n, s, t, x, y, m, u) / (1.0 - exp(-2.0 * (x-m)*(y-m)/(t-s)));
  }
  if(delT == 2) {
    if(!(xI & yI)) {
      rval = _eadel2(n, s, t, m, max(x,y), u);
    }
    else {
      rval = 0;
    }
  }
  
  return(max(0.0,min(1.0,rval)));
}

double _eadel2(unsigned long n, double s, double t, double m, double xoy, double u) {
  if(max(xoy-u,m-xoy)>=0)
    return 0;

  if(n == 1)
    return 1.0 - _eapsi(1,s,t,m,xoy,u)/(xoy-m);

  double rval = 0;             
  switch(n%2) {
  case 0: //Even
    for(unsigned long j=1 ; j <= (n/2); ++j) {
      rval += _eapsi(j,s,t,m,xoy,u) - _eachi(j,s,t,m,xoy,u);
    }
    rval = 1.0 - rval/(xoy-m);
    break;

  case 1: //Odd
    for(unsigned long j=1 ; j <= ((n-1)/2); ++j) {
      rval += _eapsi(j,s,t,m,xoy,u) - _eachi(j,s,t,m,xoy,u);
    }
    rval += _eapsi((n+1)/2,s,t,m,xoy,u);
    rval = 1.0 - rval/(xoy-m);
    break;
  }

  return(rval);
}

double _eagamma(unsigned long n, double s, double t, double x, double y, double L, double U) {
  return (1.0 - _eazeta(n, s, t, x, y, L, U));
}

double _eapsi(unsigned long j, double s, double t, double m, double xoy, double u) {
  double aum   = abs(u-m);
  double P     = -2.0*aum*j/(t-s);
  double rval  = (2.0*aum*j-(xoy-m))*exp(P*(aum*j-(xoy-m)));

  return(rval);
}

double _eazeta(unsigned long n, double s, double t, double x, double y, double L, double U) {
  if( (max(x-U, y-U) >=0) || (max(L-x, L-y) >= 0))
    return(1); 
  
  double D = U-L;
  double P = -2 / (t-s);
  double Ds = D*D;
    
  double rval = 0;
  unsigned long j, jc = ceil(static_cast<double>(n)/2.0);
  double js;
  for(j=1; j < jc; ++j) {
    js = j*j;
    rval += exp(P*((D*j+L)-x)*((D*j+L)-y))
      + exp(P*((D*j-U)+x)*((D*j-U)+y))
      - exp(P*js*Ds-(P*j*D*(y-x)))
      - exp(P*js*Ds+(P*j*D*(y-x)));
  }
  j = jc;

  rval += exp(P*((D*j+L)-x)*((D*j+L)-y))
        + exp(P*((D*j-U)+x)*((D*j-U)+y));

  if((n%2)) { //i.e. if n is odd
    return(rval);
  }

  js = j*j;
  rval -= exp(P*(js)*Ds-(P*j*D*(y-x)))
        + exp(P*(js)*Ds+(P*j*D*(y-x)));
  
  return(rval);
}


// [[Rcpp::export]]
extern "C" double eachi_cpp(SEXP jS, SEXP sS, SEXP tS, SEXP mS, SEXP xoyS, SEXP uS) {
  unsigned long j = Rcpp::as<unsigned long>(jS);
  double s = Rcpp::as<double>(sS),
    t = Rcpp::as<double>(tS),
    m = Rcpp::as<double>(mS),
    xoy = Rcpp::as<double>(xoyS),
    u = Rcpp::as<double>(uS);
  
  return(_eachi(j,s,t,m,xoy,u));
}

// [[Rcpp::export]]
extern "C" double eadel_cpp(SEXP nS, SEXP sS, SEXP tS, SEXP xS, SEXP yS,  SEXP mS, SEXP uS) {
  unsigned long n = Rcpp::as<unsigned long>(nS);
  double s = Rcpp::as<double>(sS),
    t = Rcpp::as<double>(tS),
    x = Rcpp::as<double>(xS),
    y = Rcpp::as<double>(yS),
    m = Rcpp::as<double>(mS),
    u = Rcpp::as<double>(uS);

return(_eadel(n,s,t,x,y,m,u));
}

// [[Rcpp::export]]
extern "C" Rcpp::NumericVector eadel_pair_cpp(SEXP nS, SEXP sS, SEXP tS, SEXP xS, SEXP yS,  SEXP mS, SEXP uS) {
  unsigned long n = Rcpp::as<unsigned long>(nS);
  double s = Rcpp::as<double>(sS),
    t = Rcpp::as<double>(tS),
    x = Rcpp::as<double>(xS),
    y = Rcpp::as<double>(yS),
    m = Rcpp::as<double>(mS),
    u = Rcpp::as<double>(uS);

  return Rcpp::NumericVector::create(Rcpp::Named("s1")=_eadel(n,s,t,x,y,m,u),
				     Rcpp::Named("s2")=_eadel(n+1,s,t,x,y,m,u));
}

// [[Rcpp::export]]
extern "C" double eadel2_cpp(SEXP nS, SEXP sS, SEXP tS, SEXP mS, SEXP xoyS, SEXP uS) {
  unsigned long n = Rcpp::as<unsigned long>(nS);
  double s = Rcpp::as<double>(sS),
    t = Rcpp::as<double>(tS),
    m = Rcpp::as<double>(mS),
    xoy = Rcpp::as<double>(xoyS),
    u = Rcpp::as<double>(uS);
  
  return(_eadel2(n,s,t,m,xoy,u));  
}

// [[Rcpp::export]]
extern "C" double eagamma_cpp(SEXP nS, SEXP sS, SEXP tS, SEXP xS, SEXP yS, SEXP LS, SEXP US){
  unsigned long n = Rcpp::as<unsigned long>(nS);
  double s = Rcpp::as<double>(sS),
    t = Rcpp::as<double>(tS),
    x = Rcpp::as<double>(xS),
    y = Rcpp::as<double>(yS),
    L = Rcpp::as<double>(LS),
    U = Rcpp::as<double>(US);

  return static_cast<double>(1 - _eazeta(n,s,t,x,y,L,U));  
}

// [[Rcpp::export]]
extern "C" double eapsi_cpp(SEXP jS, SEXP sS, SEXP tS, SEXP mS, SEXP xoyS, SEXP uS) {
  unsigned long j = Rcpp::as<unsigned long>(jS);
  double s = Rcpp::as<double>(sS),
    t = Rcpp::as<double>(tS),
    m = Rcpp::as<double>(mS),
    xoy = Rcpp::as<double>(xoyS),
    u = Rcpp::as<double>(uS);
  
  return(_eapsi(j,s,t,m,xoy,u));
}

// [[Rcpp::export]]
extern "C" double eazeta_cpp(SEXP nS, SEXP sS, SEXP tS, SEXP xS, SEXP yS, SEXP LS, SEXP US){

  unsigned long n = Rcpp::as<unsigned long>(nS);
  double s = Rcpp::as<double>(sS),
    t = Rcpp::as<double>(tS),
    x = Rcpp::as<double>(xS),
    y = Rcpp::as<double>(yS),
    L = Rcpp::as<double>(LS),
    U = Rcpp::as<double>(US);

  return(_eazeta(n,s,t,x,y,L,U));
}
