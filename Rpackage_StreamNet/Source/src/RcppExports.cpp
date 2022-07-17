// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getDistance
double getDistance(CharacterMatrix flow_path_nodeIDs, NumericVector stream_path_dist, CharacterMatrix total_edgelist);
RcppExport SEXP _StreamNet_getDistance(SEXP flow_path_nodeIDsSEXP, SEXP stream_path_distSEXP, SEXP total_edgelistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterMatrix >::type flow_path_nodeIDs(flow_path_nodeIDsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type stream_path_dist(stream_path_distSEXP);
    Rcpp::traits::input_parameter< CharacterMatrix >::type total_edgelist(total_edgelistSEXP);
    rcpp_result_gen = Rcpp::wrap(getDistance(flow_path_nodeIDs, stream_path_dist, total_edgelist));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StreamNet_getDistance", (DL_FUNC) &_StreamNet_getDistance, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_StreamNet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
