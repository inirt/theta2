// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// lognorm_dens_dx
double lognorm_dens_dx(double x, double mean, double sd);
RcppExport SEXP _inirt_lognorm_dens_dx(SEXP xSEXP, SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(lognorm_dens_dx(x, mean, sd));
    return rcpp_result_gen;
END_RCPP
}
// logtruncnorm_dens_dx
double logtruncnorm_dens_dx(double x, double mean, double sd, double a, double b);
RcppExport SEXP _inirt_logtruncnorm_dens_dx(SEXP xSEXP, SEXP meanSEXP, SEXP sdSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(logtruncnorm_dens_dx(x, mean, sd, a, b));
    return rcpp_result_gen;
END_RCPP
}
// log_normd_dx
double log_normd_dx(double x, double mean, double sd);
RcppExport SEXP _inirt_log_normd_dx(SEXP xSEXP, SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(log_normd_dx(x, mean, sd));
    return rcpp_result_gen;
END_RCPP
}
// normd_dx
double normd_dx(double x, double mean, double sd);
RcppExport SEXP _inirt_normd_dx(SEXP xSEXP, SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(normd_dx(x, mean, sd));
    return rcpp_result_gen;
END_RCPP
}
// log_unifd_dx
double log_unifd_dx(double x, double l, double u);
RcppExport SEXP _inirt_log_unifd_dx(SEXP xSEXP, SEXP lSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(log_unifd_dx(x, l, u));
    return rcpp_result_gen;
END_RCPP
}
// unifd_dx
double unifd_dx(double x, double l, double u);
RcppExport SEXP _inirt_unifd_dx(SEXP xSEXP, SEXP lSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    Rcpp::traits::input_parameter< double >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(unifd_dx(x, l, u));
    return rcpp_result_gen;
END_RCPP
}
// logadd
double logadd(arma::vec x);
RcppExport SEXP _inirt_logadd(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logadd(x));
    return rcpp_result_gen;
END_RCPP
}
// r8poly_value_horner
double r8poly_value_horner(int m, arma::vec c, double x);
RcppExport SEXP _inirt_r8poly_value_horner(SEXP mSEXP, SEXP cSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(r8poly_value_horner(m, c, x));
    return rcpp_result_gen;
END_RCPP
}
// r8_uniform_01
double r8_uniform_01(int& seed);
RcppExport SEXP _inirt_r8_uniform_01(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int& >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(r8_uniform_01(seed));
    return rcpp_result_gen;
END_RCPP
}
// normal_01_sample
double normal_01_sample(int& seed);
RcppExport SEXP _inirt_normal_01_sample(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int& >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_01_sample(seed));
    return rcpp_result_gen;
END_RCPP
}
// normal_01_cdf
double normal_01_cdf(double x);
RcppExport SEXP _inirt_normal_01_cdf(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_01_cdf(x));
    return rcpp_result_gen;
END_RCPP
}
// normal_01_cdf_inv
double normal_01_cdf_inv(double p);
RcppExport SEXP _inirt_normal_01_cdf_inv(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_01_cdf_inv(p));
    return rcpp_result_gen;
END_RCPP
}
// truncated_normal_ab_sample
double truncated_normal_ab_sample(double mu, double sigma, double a, double b, int& seed);
RcppExport SEXP _inirt_truncated_normal_ab_sample(SEXP muSEXP, SEXP sigmaSEXP, SEXP aSEXP, SEXP bSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int& >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(truncated_normal_ab_sample(mu, sigma, a, b, seed));
    return rcpp_result_gen;
END_RCPP
}
// normal_cdf_inv
double normal_cdf_inv(double cdf, double mu, double sigma);
RcppExport SEXP _inirt_normal_cdf_inv(SEXP cdfSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type cdf(cdfSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_cdf_inv(cdf, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// log_normal_cdf_inv
double log_normal_cdf_inv(double cdf, double mu, double sigma);
RcppExport SEXP _inirt_log_normal_cdf_inv(SEXP cdfSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type cdf(cdfSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(log_normal_cdf_inv(cdf, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// log_normal_sample
double log_normal_sample(double mu, double sigma, int& seed);
RcppExport SEXP _inirt_log_normal_sample(SEXP muSEXP, SEXP sigmaSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int& >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(log_normal_sample(mu, sigma, seed));
    return rcpp_result_gen;
END_RCPP
}
// sis4
arma::vec sis4(arma::mat& data, int n, double tol);
RcppExport SEXP _inirt_sis4(SEXP dataSEXP, SEXP nSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(sis4(data, n, tol));
    return rcpp_result_gen;
END_RCPP
}
// sis5
arma::vec sis5(arma::mat& data, int n_dimensions, arma::vec dimension_start, arma::vec dimension_end, int n_second_order, int n, double tol);
RcppExport SEXP _inirt_sis5(SEXP dataSEXP, SEXP n_dimensionsSEXP, SEXP dimension_startSEXP, SEXP dimension_endSEXP, SEXP n_second_orderSEXP, SEXP nSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type n_dimensions(n_dimensionsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dimension_start(dimension_startSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type dimension_end(dimension_endSEXP);
    Rcpp::traits::input_parameter< int >::type n_second_order(n_second_orderSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(sis5(data, n_dimensions, dimension_start, dimension_end, n_second_order, n, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_inirt_lognorm_dens_dx", (DL_FUNC) &_inirt_lognorm_dens_dx, 3},
    {"_inirt_logtruncnorm_dens_dx", (DL_FUNC) &_inirt_logtruncnorm_dens_dx, 5},
    {"_inirt_log_normd_dx", (DL_FUNC) &_inirt_log_normd_dx, 3},
    {"_inirt_normd_dx", (DL_FUNC) &_inirt_normd_dx, 3},
    {"_inirt_log_unifd_dx", (DL_FUNC) &_inirt_log_unifd_dx, 3},
    {"_inirt_unifd_dx", (DL_FUNC) &_inirt_unifd_dx, 3},
    {"_inirt_logadd", (DL_FUNC) &_inirt_logadd, 1},
    {"_inirt_r8poly_value_horner", (DL_FUNC) &_inirt_r8poly_value_horner, 3},
    {"_inirt_r8_uniform_01", (DL_FUNC) &_inirt_r8_uniform_01, 1},
    {"_inirt_normal_01_sample", (DL_FUNC) &_inirt_normal_01_sample, 1},
    {"_inirt_normal_01_cdf", (DL_FUNC) &_inirt_normal_01_cdf, 1},
    {"_inirt_normal_01_cdf_inv", (DL_FUNC) &_inirt_normal_01_cdf_inv, 1},
    {"_inirt_truncated_normal_ab_sample", (DL_FUNC) &_inirt_truncated_normal_ab_sample, 5},
    {"_inirt_normal_cdf_inv", (DL_FUNC) &_inirt_normal_cdf_inv, 3},
    {"_inirt_log_normal_cdf_inv", (DL_FUNC) &_inirt_log_normal_cdf_inv, 3},
    {"_inirt_log_normal_sample", (DL_FUNC) &_inirt_log_normal_sample, 3},
    {"_inirt_sis4", (DL_FUNC) &_inirt_sis4, 3},
    {"_inirt_sis5", (DL_FUNC) &_inirt_sis5, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_inirt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
