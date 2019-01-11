// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// pairWiseMSA_cpp
Rcpp::List pairWiseMSA_cpp(CharacterVector r_sread, CharacterVector r_ref_seq, CharacterVector r_seq_names);
RcppExport SEXP _MotifBinner2_pairWiseMSA_cpp(SEXP r_sreadSEXP, SEXP r_ref_seqSEXP, SEXP r_seq_namesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type r_sread(r_sreadSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_ref_seq(r_ref_seqSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_seq_names(r_seq_namesSEXP);
    rcpp_result_gen = Rcpp::wrap(pairWiseMSA_cpp(r_sread, r_ref_seq, r_seq_names));
    return rcpp_result_gen;
END_RCPP
}
// trimFront_cpp
Rcpp::List trimFront_cpp(CharacterVector r_sread, CharacterVector r_qual, CharacterVector r_primer, std::vector<int> prefix_lens);
RcppExport SEXP _MotifBinner2_trimFront_cpp(SEXP r_sreadSEXP, SEXP r_qualSEXP, SEXP r_primerSEXP, SEXP prefix_lensSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type r_sread(r_sreadSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_qual(r_qualSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_primer(r_primerSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type prefix_lens(prefix_lensSEXP);
    rcpp_result_gen = Rcpp::wrap(trimFront_cpp(r_sread, r_qual, r_primer, prefix_lens));
    return rcpp_result_gen;
END_RCPP
}
// map_reads_no_ins_cpp
Rcpp::List map_reads_no_ins_cpp(CharacterVector r_profile, CharacterVector r_reads, CharacterVector r_quals);
RcppExport SEXP _MotifBinner2_map_reads_no_ins_cpp(SEXP r_profileSEXP, SEXP r_readsSEXP, SEXP r_qualsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type r_profile(r_profileSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_reads(r_readsSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_quals(r_qualsSEXP);
    rcpp_result_gen = Rcpp::wrap(map_reads_no_ins_cpp(r_profile, r_reads, r_quals));
    return rcpp_result_gen;
END_RCPP
}
// transfer_gaps_cpp
Rcpp::List transfer_gaps_cpp(CharacterVector aligned_read, CharacterVector r_qual, NumericVector gap_only_cols);
RcppExport SEXP _MotifBinner2_transfer_gaps_cpp(SEXP aligned_readSEXP, SEXP r_qualSEXP, SEXP gap_only_colsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type aligned_read(aligned_readSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_qual(r_qualSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gap_only_cols(gap_only_colsSEXP);
    rcpp_result_gen = Rcpp::wrap(transfer_gaps_cpp(aligned_read, r_qual, gap_only_cols));
    return rcpp_result_gen;
END_RCPP
}
// gapQualityTweaker_ol_cpp
Rcpp::List gapQualityTweaker_ol_cpp(CharacterVector reads, NumericMatrix q_mat);
RcppExport SEXP _MotifBinner2_gapQualityTweaker_ol_cpp(SEXP readsSEXP, SEXP q_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type reads(readsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type q_mat(q_matSEXP);
    rcpp_result_gen = Rcpp::wrap(gapQualityTweaker_ol_cpp(reads, q_mat));
    return rcpp_result_gen;
END_RCPP
}
// gapQualityTweaker_non_ol_cpp
Rcpp::List gapQualityTweaker_non_ol_cpp(CharacterVector reads, NumericMatrix q_mat, std::string which_pair, NumericVector avg_quals);
RcppExport SEXP _MotifBinner2_gapQualityTweaker_non_ol_cpp(SEXP readsSEXP, SEXP q_matSEXP, SEXP which_pairSEXP, SEXP avg_qualsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type reads(readsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type q_mat(q_matSEXP);
    Rcpp::traits::input_parameter< std::string >::type which_pair(which_pairSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type avg_quals(avg_qualsSEXP);
    rcpp_result_gen = Rcpp::wrap(gapQualityTweaker_non_ol_cpp(reads, q_mat, which_pair, avg_quals));
    return rcpp_result_gen;
END_RCPP
}
// scoreAlignmentPositions_cpp
Rcpp::List scoreAlignmentPositions_cpp(CharacterVector reads, NumericMatrix q_mat);
RcppExport SEXP _MotifBinner2_scoreAlignmentPositions_cpp(SEXP readsSEXP, SEXP q_matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type reads(readsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type q_mat(q_matSEXP);
    rcpp_result_gen = Rcpp::wrap(scoreAlignmentPositions_cpp(reads, q_mat));
    return rcpp_result_gen;
END_RCPP
}
// buildConsensus_cpp
Rcpp::List buildConsensus_cpp(NumericMatrix score_mat, double required_dominance, double minimum_score);
RcppExport SEXP _MotifBinner2_buildConsensus_cpp(SEXP score_matSEXP, SEXP required_dominanceSEXP, SEXP minimum_scoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type score_mat(score_matSEXP);
    Rcpp::traits::input_parameter< double >::type required_dominance(required_dominanceSEXP);
    Rcpp::traits::input_parameter< double >::type minimum_score(minimum_scoreSEXP);
    rcpp_result_gen = Rcpp::wrap(buildConsensus_cpp(score_mat, required_dominance, minimum_score));
    return rcpp_result_gen;
END_RCPP
}
// tallyPrimerSeqErrors_cpp
std::map<char, std::map<char, std::map<char, int> > > tallyPrimerSeqErrors_cpp(CharacterVector r_sread, CharacterVector r_primer, CharacterVector r_qual);
RcppExport SEXP _MotifBinner2_tallyPrimerSeqErrors_cpp(SEXP r_sreadSEXP, SEXP r_primerSEXP, SEXP r_qualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type r_sread(r_sreadSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_primer(r_primerSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_qual(r_qualSEXP);
    rcpp_result_gen = Rcpp::wrap(tallyPrimerSeqErrors_cpp(r_sread, r_primer, r_qual));
    return rcpp_result_gen;
END_RCPP
}
// regionSplit_cpp
Rcpp::List regionSplit_cpp(CharacterVector mapped_read, CharacterVector profile, CharacterVector region_map, CharacterVector mapped_qual);
RcppExport SEXP _MotifBinner2_regionSplit_cpp(SEXP mapped_readSEXP, SEXP profileSEXP, SEXP region_mapSEXP, SEXP mapped_qualSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type mapped_read(mapped_readSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type profile(profileSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type region_map(region_mapSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type mapped_qual(mapped_qualSEXP);
    rcpp_result_gen = Rcpp::wrap(regionSplit_cpp(mapped_read, profile, region_map, mapped_qual));
    return rcpp_result_gen;
END_RCPP
}
// removeChars_cpp
Rcpp::List removeChars_cpp(CharacterVector r_sread, CharacterVector r_qual, std::string char_to_strip);
RcppExport SEXP _MotifBinner2_removeChars_cpp(SEXP r_sreadSEXP, SEXP r_qualSEXP, SEXP char_to_stripSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type r_sread(r_sreadSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type r_qual(r_qualSEXP);
    Rcpp::traits::input_parameter< std::string >::type char_to_strip(char_to_stripSEXP);
    rcpp_result_gen = Rcpp::wrap(removeChars_cpp(r_sread, r_qual, char_to_strip));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MotifBinner2_pairWiseMSA_cpp", (DL_FUNC) &_MotifBinner2_pairWiseMSA_cpp, 3},
    {"_MotifBinner2_trimFront_cpp", (DL_FUNC) &_MotifBinner2_trimFront_cpp, 4},
    {"_MotifBinner2_map_reads_no_ins_cpp", (DL_FUNC) &_MotifBinner2_map_reads_no_ins_cpp, 3},
    {"_MotifBinner2_transfer_gaps_cpp", (DL_FUNC) &_MotifBinner2_transfer_gaps_cpp, 3},
    {"_MotifBinner2_gapQualityTweaker_ol_cpp", (DL_FUNC) &_MotifBinner2_gapQualityTweaker_ol_cpp, 2},
    {"_MotifBinner2_gapQualityTweaker_non_ol_cpp", (DL_FUNC) &_MotifBinner2_gapQualityTweaker_non_ol_cpp, 4},
    {"_MotifBinner2_scoreAlignmentPositions_cpp", (DL_FUNC) &_MotifBinner2_scoreAlignmentPositions_cpp, 2},
    {"_MotifBinner2_buildConsensus_cpp", (DL_FUNC) &_MotifBinner2_buildConsensus_cpp, 3},
    {"_MotifBinner2_tallyPrimerSeqErrors_cpp", (DL_FUNC) &_MotifBinner2_tallyPrimerSeqErrors_cpp, 3},
    {"_MotifBinner2_regionSplit_cpp", (DL_FUNC) &_MotifBinner2_regionSplit_cpp, 4},
    {"_MotifBinner2_removeChars_cpp", (DL_FUNC) &_MotifBinner2_removeChars_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_MotifBinner2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
