// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bal
std::vector<double> bal(const std::vector<double> dbh, const std::vector<double> expansion, const bool imperial_units);
RcppExport SEXP _biometrics_utilities_bal(SEXP dbhSEXP, SEXP expansionSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< const bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(bal(dbh, expansion, imperial_units));
    return rcpp_result_gen;
END_RCPP
}
// ccfl
std::vector<double> ccfl(const std::vector<double> dbh, const std::vector<double> mcw, const std::vector<double> expansion, const bool imperial_units);
RcppExport SEXP _biometrics_utilities_ccfl(SEXP dbhSEXP, SEXP mcwSEXP, SEXP expansionSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type mcw(mcwSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< const bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(ccfl(dbh, mcw, expansion, imperial_units));
    return rcpp_result_gen;
END_RCPP
}
// cch
std::vector<double> cch(const std::vector<int>& species, const std::vector<double>& dbh, const std::vector<double>& height, const std::vector<double>& crown_length, const std::vector<double>& dacb, const std::vector<double>& lcw, const std::vector<double>& expansion, Rcpp::DataFrame& parameters, const bool imperial_units);
RcppExport SEXP _biometrics_utilities_cch(SEXP speciesSEXP, SEXP dbhSEXP, SEXP heightSEXP, SEXP crown_lengthSEXP, SEXP dacbSEXP, SEXP lcwSEXP, SEXP expansionSEXP, SEXP parametersSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type species(speciesSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type height(heightSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type crown_length(crown_lengthSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dacb(dacbSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type lcw(lcwSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(cch(species, dbh, height, crown_length, dacb, lcw, expansion, parameters, imperial_units));
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
// qmd
double qmd(const std::vector<double> dbh, const std::vector<double> expansion);
RcppExport SEXP _biometrics_utilities_qmd(SEXP dbhSEXP, SEXP expansionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    rcpp_result_gen = Rcpp::wrap(qmd(dbh, expansion));
    return rcpp_result_gen;
END_RCPP
}
// relative_spacing
double relative_spacing(const std::vector<double> expansion, const double dominant_height, const bool imperial);
RcppExport SEXP _biometrics_utilities_relative_spacing(SEXP expansionSEXP, SEXP dominant_heightSEXP, SEXP imperialSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< const double >::type dominant_height(dominant_heightSEXP);
    Rcpp::traits::input_parameter< const bool >::type imperial(imperialSEXP);
    rcpp_result_gen = Rcpp::wrap(relative_spacing(expansion, dominant_height, imperial));
    return rcpp_result_gen;
END_RCPP
}
// curtis_rd
double curtis_rd(const std::vector<double> dbh, const std::vector<double> expansion, const bool imperial_units);
RcppExport SEXP _biometrics_utilities_curtis_rd(SEXP dbhSEXP, SEXP expansionSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< const bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(curtis_rd(dbh, expansion, imperial_units));
    return rcpp_result_gen;
END_RCPP
}
// reineke_sdi
double reineke_sdi(const std::vector<double> dbh, const std::vector<double> expansion, const bool imperial_units);
RcppExport SEXP _biometrics_utilities_reineke_sdi(SEXP dbhSEXP, SEXP expansionSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< const bool >::type imperial_units(imperial_unitsSEXP);
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
// Clark_Evans_R
double Clark_Evans_R(const std::vector<double>& x, const std::vector<double>& y, double plotarea, const std::vector<double>& poly_x, const std::vector<double>& poly_y);
RcppExport SEXP _biometrics_utilities_Clark_Evans_R(SEXP xSEXP, SEXP ySEXP, SEXP plotareaSEXP, SEXP poly_xSEXP, SEXP poly_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type plotarea(plotareaSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type poly_x(poly_xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type poly_y(poly_ySEXP);
    rcpp_result_gen = Rcpp::wrap(Clark_Evans_R(x, y, plotarea, poly_x, poly_y));
    return rcpp_result_gen;
END_RCPP
}
// Clark_Evans_R_circle
double Clark_Evans_R_circle(const std::vector<double>& x, const std::vector<double>& y, double plotarea, const double plot_center_x, const double plot_center_y, const double plot_radius);
RcppExport SEXP _biometrics_utilities_Clark_Evans_R_circle(SEXP xSEXP, SEXP ySEXP, SEXP plotareaSEXP, SEXP plot_center_xSEXP, SEXP plot_center_ySEXP, SEXP plot_radiusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type plotarea(plotareaSEXP);
    Rcpp::traits::input_parameter< const double >::type plot_center_x(plot_center_xSEXP);
    Rcpp::traits::input_parameter< const double >::type plot_center_y(plot_center_ySEXP);
    Rcpp::traits::input_parameter< const double >::type plot_radius(plot_radiusSEXP);
    rcpp_result_gen = Rcpp::wrap(Clark_Evans_R_circle(x, y, plotarea, plot_center_x, plot_center_y, plot_radius));
    return rcpp_result_gen;
END_RCPP
}
// Hegyi
std::vector<double> Hegyi(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double> dbh, const std::vector<double>& poly_x, const std::vector<double>& poly_y, const bool imperial_units);
RcppExport SEXP _biometrics_utilities_Hegyi(SEXP xSEXP, SEXP ySEXP, SEXP dbhSEXP, SEXP poly_xSEXP, SEXP poly_ySEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::vector<double> >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type poly_x(poly_xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type poly_y(poly_ySEXP);
    Rcpp::traits::input_parameter< const bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(Hegyi(x, y, dbh, poly_x, poly_y, imperial_units));
    return rcpp_result_gen;
END_RCPP
}
// Arney_CSI
std::vector<double> Arney_CSI(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& dbh, const std::vector<double>& mcw);
RcppExport SEXP _biometrics_utilities_Arney_CSI(SEXP xSEXP, SEXP ySEXP, SEXP dbhSEXP, SEXP mcwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type mcw(mcwSEXP);
    rcpp_result_gen = Rcpp::wrap(Arney_CSI(x, y, dbh, mcw));
    return rcpp_result_gen;
END_RCPP
}
// mcw
std::vector<double> mcw(const std::vector<int>& fia, const std::vector<double>& dbh, const bool imperial_units, const int default_fia);
RcppExport SEXP _biometrics_utilities_mcw(SEXP fiaSEXP, SEXP dbhSEXP, SEXP imperial_unitsSEXP, SEXP default_fiaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type fia(fiaSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const bool >::type imperial_units(imperial_unitsSEXP);
    Rcpp::traits::input_parameter< const int >::type default_fia(default_fiaSEXP);
    rcpp_result_gen = Rcpp::wrap(mcw(fia, dbh, imperial_units, default_fia));
    return rcpp_result_gen;
END_RCPP
}
// mcw_species
Rcpp::DataFrame mcw_species();
RcppExport SEXP _biometrics_utilities_mcw_species() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(mcw_species());
    return rcpp_result_gen;
END_RCPP
}
// hd_fit
Rcpp::DataFrame hd_fit(const std::vector<int>& fia, const std::vector<double>& dbh, const std::vector<double>& height, const double bh);
RcppExport SEXP _biometrics_utilities_hd_fit(SEXP fiaSEXP, SEXP dbhSEXP, SEXP heightSEXP, SEXP bhSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type fia(fiaSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type height(heightSEXP);
    Rcpp::traits::input_parameter< const double >::type bh(bhSEXP);
    rcpp_result_gen = Rcpp::wrap(hd_fit(fia, dbh, height, bh));
    return rcpp_result_gen;
END_RCPP
}
// hd_predict
std::vector<double> hd_predict(Rcpp::DataFrame& hd_parameters, const std::vector<int>& fia, const std::vector<double>& dbh, const double bh);
RcppExport SEXP _biometrics_utilities_hd_predict(SEXP hd_parametersSEXP, SEXP fiaSEXP, SEXP dbhSEXP, SEXP bhSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame& >::type hd_parameters(hd_parametersSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type fia(fiaSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const double >::type bh(bhSEXP);
    rcpp_result_gen = Rcpp::wrap(hd_predict(hd_parameters, fia, dbh, bh));
    return rcpp_result_gen;
END_RCPP
}
// Glover_Hool
std::vector<double> Glover_Hool(const std::vector<double>& dbh, const std::vector<double>& expansion, const bool use_arithmetic, const bool imperial_units);
RcppExport SEXP _biometrics_utilities_Glover_Hool(SEXP dbhSEXP, SEXP expansionSEXP, SEXP use_arithmeticSEXP, SEXP imperial_unitsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type expansion(expansionSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_arithmetic(use_arithmeticSEXP);
    Rcpp::traits::input_parameter< const bool >::type imperial_units(imperial_unitsSEXP);
    rcpp_result_gen = Rcpp::wrap(Glover_Hool(dbh, expansion, use_arithmetic, imperial_units));
    return rcpp_result_gen;
END_RCPP
}
// APA
std::vector<double> APA(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& dbh, const std::vector<double>& plot_x, const std::vector<double>& plot_y, const bool weighted);
RcppExport SEXP _biometrics_utilities_APA(SEXP xSEXP, SEXP ySEXP, SEXP dbhSEXP, SEXP plot_xSEXP, SEXP plot_ySEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type plot_x(plot_xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type plot_y(plot_ySEXP);
    Rcpp::traits::input_parameter< const bool >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(APA(x, y, dbh, plot_x, plot_y, weighted));
    return rcpp_result_gen;
END_RCPP
}
// APA_Polygons
Rcpp::DataFrame APA_Polygons(const std::vector<int>& tree_id, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& dbh, const std::vector<double>& plot_x, const std::vector<double>& plot_y, const bool weighted);
RcppExport SEXP _biometrics_utilities_APA_Polygons(SEXP tree_idSEXP, SEXP xSEXP, SEXP ySEXP, SEXP dbhSEXP, SEXP plot_xSEXP, SEXP plot_ySEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type tree_id(tree_idSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type plot_x(plot_xSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type plot_y(plot_ySEXP);
    Rcpp::traits::input_parameter< const bool >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(APA_Polygons(tree_id, x, y, dbh, plot_x, plot_y, weighted));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_biometrics_utilities_bal", (DL_FUNC) &_biometrics_utilities_bal, 3},
    {"_biometrics_utilities_ccfl", (DL_FUNC) &_biometrics_utilities_ccfl, 4},
    {"_biometrics_utilities_cch", (DL_FUNC) &_biometrics_utilities_cch, 9},
    {"_biometrics_utilities_dominant_height", (DL_FUNC) &_biometrics_utilities_dominant_height, 5},
    {"_biometrics_utilities_qmd", (DL_FUNC) &_biometrics_utilities_qmd, 2},
    {"_biometrics_utilities_relative_spacing", (DL_FUNC) &_biometrics_utilities_relative_spacing, 3},
    {"_biometrics_utilities_curtis_rd", (DL_FUNC) &_biometrics_utilities_curtis_rd, 3},
    {"_biometrics_utilities_reineke_sdi", (DL_FUNC) &_biometrics_utilities_reineke_sdi, 3},
    {"_biometrics_utilities_ccf", (DL_FUNC) &_biometrics_utilities_ccf, 3},
    {"_biometrics_utilities_Clark_Evans_R", (DL_FUNC) &_biometrics_utilities_Clark_Evans_R, 5},
    {"_biometrics_utilities_Clark_Evans_R_circle", (DL_FUNC) &_biometrics_utilities_Clark_Evans_R_circle, 6},
    {"_biometrics_utilities_Hegyi", (DL_FUNC) &_biometrics_utilities_Hegyi, 6},
    {"_biometrics_utilities_Arney_CSI", (DL_FUNC) &_biometrics_utilities_Arney_CSI, 4},
    {"_biometrics_utilities_mcw", (DL_FUNC) &_biometrics_utilities_mcw, 4},
    {"_biometrics_utilities_mcw_species", (DL_FUNC) &_biometrics_utilities_mcw_species, 0},
    {"_biometrics_utilities_hd_fit", (DL_FUNC) &_biometrics_utilities_hd_fit, 4},
    {"_biometrics_utilities_hd_predict", (DL_FUNC) &_biometrics_utilities_hd_predict, 4},
    {"_biometrics_utilities_Glover_Hool", (DL_FUNC) &_biometrics_utilities_Glover_Hool, 4},
    {"_biometrics_utilities_APA", (DL_FUNC) &_biometrics_utilities_APA, 6},
    {"_biometrics_utilities_APA_Polygons", (DL_FUNC) &_biometrics_utilities_APA_Polygons, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_biometrics_utilities(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
