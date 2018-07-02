// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// optimMMTC
List optimMMTC(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd ThetaX, const Eigen::MatrixXd K, const Eigen::MatrixXd A, Eigen::MatrixXd etainit, int iter, double eigvalthresh, int numexcessthresh, bool calcGradHess);
RcppExport SEXP _mongrel_optimMMTC(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaXSEXP, SEXP KSEXP, SEXP ASEXP, SEXP etainitSEXP, SEXP iterSEXP, SEXP eigvalthreshSEXP, SEXP numexcessthreshSEXP, SEXP calcGradHessSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type ThetaX(ThetaXSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type etainit(etainitSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type eigvalthresh(eigvalthreshSEXP);
    Rcpp::traits::input_parameter< int >::type numexcessthresh(numexcessthreshSEXP);
    Rcpp::traits::input_parameter< bool >::type calcGradHess(calcGradHessSEXP);
    rcpp_result_gen = Rcpp::wrap(optimMMTC(Y, upsilon, ThetaX, K, A, etainit, iter, eigvalthresh, numexcessthresh, calcGradHess));
    return rcpp_result_gen;
END_RCPP
}
// loglikMMTC
double loglikMMTC(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd ThetaX, const Eigen::MatrixXd K, const Eigen::MatrixXd A, Eigen::MatrixXd etainit);
RcppExport SEXP _mongrel_loglikMMTC(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaXSEXP, SEXP KSEXP, SEXP ASEXP, SEXP etainitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type ThetaX(ThetaXSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type etainit(etainitSEXP);
    rcpp_result_gen = Rcpp::wrap(loglikMMTC(Y, upsilon, ThetaX, K, A, etainit));
    return rcpp_result_gen;
END_RCPP
}
// gradMMTC
Eigen::MatrixXd gradMMTC(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd ThetaX, const Eigen::MatrixXd K, const Eigen::MatrixXd A, Eigen::MatrixXd etainit);
RcppExport SEXP _mongrel_gradMMTC(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaXSEXP, SEXP KSEXP, SEXP ASEXP, SEXP etainitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type ThetaX(ThetaXSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type etainit(etainitSEXP);
    rcpp_result_gen = Rcpp::wrap(gradMMTC(Y, upsilon, ThetaX, K, A, etainit));
    return rcpp_result_gen;
END_RCPP
}
// hessMMTC
Eigen::MatrixXd hessMMTC(const Eigen::ArrayXXd Y, const double upsilon, const Eigen::MatrixXd ThetaX, const Eigen::MatrixXd K, const Eigen::MatrixXd A, Eigen::MatrixXd etainit);
RcppExport SEXP _mongrel_hessMMTC(SEXP YSEXP, SEXP upsilonSEXP, SEXP ThetaXSEXP, SEXP KSEXP, SEXP ASEXP, SEXP etainitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type ThetaX(ThetaXSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type etainit(etainitSEXP);
    rcpp_result_gen = Rcpp::wrap(hessMMTC(Y, upsilon, ThetaX, K, A, etainit));
    return rcpp_result_gen;
END_RCPP
}
// uncollapseMMTC
List uncollapseMMTC(const Eigen::Map<Eigen::VectorXd> eta, const Eigen::Map<Eigen::MatrixXd> X, const Eigen::Map<Eigen::MatrixXd> Theta, const Eigen::Map<Eigen::MatrixXd> Gamma, const Eigen::Map<Eigen::MatrixXd> Xi, const double upsilon);
RcppExport SEXP _mongrel_uncollapseMMTC(SEXP etaSEXP, SEXP XSEXP, SEXP ThetaSEXP, SEXP GammaSEXP, SEXP XiSEXP, SEXP upsilonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::VectorXd> >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Gamma(GammaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type Xi(XiSEXP);
    Rcpp::traits::input_parameter< const double >::type upsilon(upsilonSEXP);
    rcpp_result_gen = Rcpp::wrap(uncollapseMMTC(eta, X, Theta, Gamma, Xi, upsilon));
    return rcpp_result_gen;
END_RCPP
}
// rMatNormalCholesky_test
Eigen::MatrixXd rMatNormalCholesky_test(Eigen::MatrixXd M, Eigen::MatrixXd LU, Eigen::MatrixXd LV);
RcppExport SEXP _mongrel_rMatNormalCholesky_test(SEXP MSEXP, SEXP LUSEXP, SEXP LVSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type M(MSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type LU(LUSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type LV(LVSEXP);
    rcpp_result_gen = Rcpp::wrap(rMatNormalCholesky_test(M, LU, LV));
    return rcpp_result_gen;
END_RCPP
}
// rInvWishCholesky_test
Eigen::MatrixXd rInvWishCholesky_test(int v, Eigen::MatrixXd Psi);
RcppExport SEXP _mongrel_rInvWishCholesky_test(SEXP vSEXP, SEXP PsiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type v(vSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type Psi(PsiSEXP);
    rcpp_result_gen = Rcpp::wrap(rInvWishCholesky_test(v, Psi));
    return rcpp_result_gen;
END_RCPP
}
// rMatUnitNormal_test
Eigen::MatrixXd rMatUnitNormal_test(int n, int m);
RcppExport SEXP _mongrel_rMatUnitNormal_test(SEXP nSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(rMatUnitNormal_test(n, m));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mongrel_optimMMTC", (DL_FUNC) &_mongrel_optimMMTC, 10},
    {"_mongrel_loglikMMTC", (DL_FUNC) &_mongrel_loglikMMTC, 6},
    {"_mongrel_gradMMTC", (DL_FUNC) &_mongrel_gradMMTC, 6},
    {"_mongrel_hessMMTC", (DL_FUNC) &_mongrel_hessMMTC, 6},
    {"_mongrel_uncollapseMMTC", (DL_FUNC) &_mongrel_uncollapseMMTC, 6},
    {"_mongrel_rMatNormalCholesky_test", (DL_FUNC) &_mongrel_rMatNormalCholesky_test, 3},
    {"_mongrel_rInvWishCholesky_test", (DL_FUNC) &_mongrel_rInvWishCholesky_test, 2},
    {"_mongrel_rMatUnitNormal_test", (DL_FUNC) &_mongrel_rMatUnitNormal_test, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_mongrel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
