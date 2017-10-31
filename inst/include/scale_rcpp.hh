double _eachi(unsigned long, double, double, double, double, double);
double _eadel(unsigned long, double, double, double, double, double, double);
double _eadel2(unsigned long, double, double, double, double, double);
double _eagamma(unsigned long, double, double, double, double, double, double);
double _eapsi(unsigned long, double, double, double, double, double);
double _eazeta(unsigned long, double, double, double, double, double, double);

//R Wrappers
extern "C" double eachi_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" double eadel_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" Rcpp::NumericVector eadel_pair_cpp(SEXP, SEXP, SEXP, SEXP, SEXP,  SEXP, SEXP);
extern "C" double eadel2_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" double eagamma_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" double eapsi_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern "C" double eazeta_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
