// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bal
std::vector<double> bal(const std::vector<double> dbh, const std::vector<double> expansion, bool imperial_units);
RcppExport SEXP _biometrics_utilities_bal(SEXP dbhSEXP, SEXP expansionSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(bal(dbh, expansion, imperial_units));
    return rcpp_result_gen;
END_RCPP
}
// ccfl
std::vector<double> ccfl(const std::vector<double> dbh, const std::vector<double> mcw, const std::vector<double> expansion, bool imperial_units);
RcppExport SEXP _biometrics_utilities_ccfl(SEXP dbhSEXP, SEXP mcwSEXP, SEXP expansionSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type mcw(mcwSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(ccfl(dbh, mcw, expansion, imperial_units));
    return rcpp_result_gen;
END_RCPP
}
// cch
std::vector<double> cch(const std::vector<double> dbh, const std::vector<double> height, const std::vector<double> crown_length, const std::vector<double> dacb, const std::vector<double> lcw, const std::vector<double> expansion, const std::vector<double> parameters, bool imperial_units);
RcppExport SEXP _biometrics_utilities_cch(SEXP dbhSEXP, SEXP heightSEXP, SEXP crown_lengthSEXP, SEXP dacbSEXP, SEXP lcwSEXP, SEXP expansionSEXP, SEXP parametersSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type height(heightSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type crown_length(crown_lengthSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type dacb(dacbSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type lcw(lcwSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(cch(dbh, height, crown_length, dacb, lcw, expansion, parameters, imperial_units));
    return rcpp_result_gen;
END_RCPP
}
// dominant_height
double dominant_height(const std::vector<double> height, const std::vector<double> dbh, const std::vector<double> expansion, const int dominant_cohort_size, const int method);
RcppExport SEXP _biometrics_utilities_dominant_height(SEXP heightSEXP, SEXP dbhSEXP, SEXP expansionSEXP, SEXP dominant_cohort_sizeSEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type height(heightSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< const int >::type dominant_cohort_size(dominant_cohort_sizeSEXP);
    Rcpp::traits::input_parameter< const int >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(dominant_height(height, dbh, expansion, dominant_cohort_size, method));
    return rcpp_result_gen;
END_RCPP
}
// relative_spacing
double relative_spacing(const std::vector<double> expansion, const double dominant_height, bool imperial);
RcppExport SEXP _biometrics_utilities_relative_spacing(SEXP expansionSEXP, SEXP dominant_heightSEXP, SEXP imperialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< const double >::type dominant_height(dominant_heightSEXP);
    Rcpp::traits::input_parameter< bool >::type imperial(imperialSEXP);
    rcpp_result_gen = Rcpp::wrap(relative_spacing(expansion, dominant_height, imperial));
    return rcpp_result_gen;
END_RCPP
}
// curtis_rd
double curtis_rd(const std::vector<double> dbh, const std::vector<double> expansion, bool imperial_units);
RcppExport SEXP _biometrics_utilities_curtis_rd(SEXP dbhSEXP, SEXP expansionSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(curtis_rd(dbh, expansion, imperial_units));
    return rcpp_result_gen;
END_RCPP
}
// reineke_sdi
double reineke_sdi(const std::vector<double> dbh, const std::vector<double> expansion, bool imperial_units);
RcppExport SEXP _biometrics_utilities_reineke_sdi(SEXP dbhSEXP, SEXP expansionSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(reineke_sdi(dbh, expansion, imperial_units));
    return rcpp_result_gen;
END_RCPP
}
// ccf
double ccf(const std::vector<double> crown_width, const std::vector<double> expansion, bool imperial_units);
RcppExport SEXP _biometrics_utilities_ccf(SEXP crown_widthSEXP, SEXP expansionSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type crown_width(crown_widthSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(ccf(crown_width, expansion, imperial_units));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_biometrics_utilities_bal", (DL_FUNC) &_biometrics_utilities_bal, 3},
    {"_biometrics_utilities_ccfl", (DL_FUNC) &_biometrics_utilities_ccfl, 4},
    {"_biometrics_utilities_cch", (DL_FUNC) &_biometrics_utilities_cch, 8},
    {"_biometrics_utilities_dominant_height", (DL_FUNC) &_biometrics_utilities_dominant_height, 5},
    {"_biometrics_utilities_relative_spacing", (DL_FUNC) &_biometrics_utilities_relative_spacing, 3},
    {"_biometrics_utilities_curtis_rd", (DL_FUNC) &_biometrics_utilities_curtis_rd, 3},
    {"_biometrics_utilities_reineke_sdi", (DL_FUNC) &_biometrics_utilities_reineke_sdi, 3},
    {"_biometrics_utilities_ccf", (DL_FUNC) &_biometrics_utilities_ccf, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_biometrics_utilities(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
