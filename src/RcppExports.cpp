// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// simulate_single_pd_tree_cpp
Rcpp::NumericMatrix simulate_single_pd_tree_cpp(Rcpp::NumericVector pars, float max_t, float max_N, int max_tries);
RcppExport SEXP _emphasis_simulate_single_pd_tree_cpp(SEXP parsSEXP, SEXP max_tSEXP, SEXP max_NSEXP, SEXP max_triesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< float >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< float >::type max_N(max_NSEXP);
    Rcpp::traits::input_parameter< int >::type max_tries(max_triesSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_single_pd_tree_cpp(pars, max_t, max_N, max_tries));
    return rcpp_result_gen;
END_RCPP
}
// simulate_pd_trees_cpp
Rcpp::NumericMatrix simulate_pd_trees_cpp(Rcpp::NumericVector pars, float max_t, size_t repl, float max_N);
RcppExport SEXP _emphasis_simulate_pd_trees_cpp(SEXP parsSEXP, SEXP max_tSEXP, SEXP replSEXP, SEXP max_NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< float >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< size_t >::type repl(replSEXP);
    Rcpp::traits::input_parameter< float >::type max_N(max_NSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_pd_trees_cpp(pars, max_t, repl, max_N));
    return rcpp_result_gen;
END_RCPP
}
// explore_grid_cpp
Rcpp::NumericMatrix explore_grid_cpp(Rcpp::NumericVector par1, Rcpp::NumericVector par2, Rcpp::NumericVector par3, Rcpp::NumericVector par4, float max_t, int num_repl, int max_N);
RcppExport SEXP _emphasis_explore_grid_cpp(SEXP par1SEXP, SEXP par2SEXP, SEXP par3SEXP, SEXP par4SEXP, SEXP max_tSEXP, SEXP num_replSEXP, SEXP max_NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par1(par1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par2(par2SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par3(par3SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type par4(par4SEXP);
    Rcpp::traits::input_parameter< float >::type max_t(max_tSEXP);
    Rcpp::traits::input_parameter< int >::type num_repl(num_replSEXP);
    Rcpp::traits::input_parameter< int >::type max_N(max_NSEXP);
    rcpp_result_gen = Rcpp::wrap(explore_grid_cpp(par1, par2, par3, par4, max_t, num_repl, max_N));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_mce
List rcpp_mce(const std::vector<double>& brts, const std::vector<double>& init_pars, int sample_size, int maxN, const std::string& plugin, int soc, int max_missing, double max_lambda, const std::vector<double>& lower_bound, const std::vector<double>& upper_bound, double xtol_rel, int num_threads);
RcppExport SEXP _emphasis_rcpp_mce(SEXP brtsSEXP, SEXP init_parsSEXP, SEXP sample_sizeSEXP, SEXP maxNSEXP, SEXP pluginSEXP, SEXP socSEXP, SEXP max_missingSEXP, SEXP max_lambdaSEXP, SEXP lower_boundSEXP, SEXP upper_boundSEXP, SEXP xtol_relSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type brts(brtsSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type init_pars(init_parsSEXP);
    Rcpp::traits::input_parameter< int >::type sample_size(sample_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type maxN(maxNSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type plugin(pluginSEXP);
    Rcpp::traits::input_parameter< int >::type soc(socSEXP);
    Rcpp::traits::input_parameter< int >::type max_missing(max_missingSEXP);
    Rcpp::traits::input_parameter< double >::type max_lambda(max_lambdaSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type lower_bound(lower_boundSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type upper_bound(upper_boundSEXP);
    Rcpp::traits::input_parameter< double >::type xtol_rel(xtol_relSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_mce(brts, init_pars, sample_size, maxN, plugin, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_mcem
List rcpp_mcem(const std::vector<double>& brts, const std::vector<double>& init_pars, int sample_size, int maxN, const std::string& plugin, int soc, int max_missing, double max_lambda, const std::vector<double>& lower_bound, const std::vector<double>& upper_bound, double xtol_rel, int num_threads, bool copy_trees, Nullable<Function> rconditional);
RcppExport SEXP _emphasis_rcpp_mcem(SEXP brtsSEXP, SEXP init_parsSEXP, SEXP sample_sizeSEXP, SEXP maxNSEXP, SEXP pluginSEXP, SEXP socSEXP, SEXP max_missingSEXP, SEXP max_lambdaSEXP, SEXP lower_boundSEXP, SEXP upper_boundSEXP, SEXP xtol_relSEXP, SEXP num_threadsSEXP, SEXP copy_treesSEXP, SEXP rconditionalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type brts(brtsSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type init_pars(init_parsSEXP);
    Rcpp::traits::input_parameter< int >::type sample_size(sample_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type maxN(maxNSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type plugin(pluginSEXP);
    Rcpp::traits::input_parameter< int >::type soc(socSEXP);
    Rcpp::traits::input_parameter< int >::type max_missing(max_missingSEXP);
    Rcpp::traits::input_parameter< double >::type max_lambda(max_lambdaSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type lower_bound(lower_boundSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type upper_bound(upper_boundSEXP);
    Rcpp::traits::input_parameter< double >::type xtol_rel(xtol_relSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type copy_trees(copy_treesSEXP);
    Rcpp::traits::input_parameter< Nullable<Function> >::type rconditional(rconditionalSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_mcem(brts, init_pars, sample_size, maxN, plugin, soc, max_missing, max_lambda, lower_bound, upper_bound, xtol_rel, num_threads, copy_trees, rconditional));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_mcm
List rcpp_mcm(List e_step, const std::vector<double>& init_pars, const std::string& plugin, const std::vector<double>& lower_bound, const std::vector<double>& upper_bound, double xtol_rel, int num_threads, Nullable<Function> rconditional);
RcppExport SEXP _emphasis_rcpp_mcm(SEXP e_stepSEXP, SEXP init_parsSEXP, SEXP pluginSEXP, SEXP lower_boundSEXP, SEXP upper_boundSEXP, SEXP xtol_relSEXP, SEXP num_threadsSEXP, SEXP rconditionalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type e_step(e_stepSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type init_pars(init_parsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type plugin(pluginSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type lower_bound(lower_boundSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type upper_bound(upper_boundSEXP);
    Rcpp::traits::input_parameter< double >::type xtol_rel(xtol_relSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    Rcpp::traits::input_parameter< Nullable<Function> >::type rconditional(rconditionalSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_mcm(e_step, init_pars, plugin, lower_bound, upper_bound, xtol_rel, num_threads, rconditional));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_emphasis_simulate_single_pd_tree_cpp", (DL_FUNC) &_emphasis_simulate_single_pd_tree_cpp, 4},
    {"_emphasis_simulate_pd_trees_cpp", (DL_FUNC) &_emphasis_simulate_pd_trees_cpp, 4},
    {"_emphasis_explore_grid_cpp", (DL_FUNC) &_emphasis_explore_grid_cpp, 7},
    {"_emphasis_rcpp_mce", (DL_FUNC) &_emphasis_rcpp_mce, 12},
    {"_emphasis_rcpp_mcem", (DL_FUNC) &_emphasis_rcpp_mcem, 14},
    {"_emphasis_rcpp_mcm", (DL_FUNC) &_emphasis_rcpp_mcm, 8},
    {NULL, NULL, 0}
};

void remphasis_init(DllInfo *dll);
RcppExport void R_init_emphasis(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    remphasis_init(dll);
}
