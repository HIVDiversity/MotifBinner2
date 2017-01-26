#' Builds the overlapping or non-overlapping production configs
#' @export

buildConfig <- function(overlapping,
                        fwd_file, fwd_primer_seq, fwd_primer_lens, fwd_min_score,
                        rev_file, rev_primer_seq, rev_primer_lens, rev_min_score,
                        fwd_pid_in_which_fragment, rev_pid_in_which_fragment,
                        min_read_length = 295,
                        pattern_to_chop_from_names = ' [0-9]:N:[0-9]*:[0-9]*$',
                        output_dir = "/fridge/data/MotifBinner2_test",
                        base_for_names = "CAP129_2040_009wpi_C2C3",
                        intermediate_reports = TRUE,
                        erase_history = FALSE,
                        verbosity = 3,
                        report_type = c('html'),
                        ncpu = 4,
                        bins_to_process = Inf)
{
  if (overlapping)
  {
    buildConfig_ol_prod(fwd_file = fwd_file, 
                        fwd_primer_seq = fwd_primer_seq, 
                        fwd_primer_lens = fwd_primer_lens, 
                        fwd_min_score = fwd_min_score,
                        rev_file = rev_file, 
                        rev_primer_seq = rev_primer_seq, 
                        rev_primer_lens = rev_primer_lens, 
                        rev_min_score = rev_min_score,
                        fwd_pid_in_which_fragment = fwd_pid_in_which_fragment, 
                        rev_pid_in_which_fragment = rev_pid_in_which_fragment,
                        min_read_length = min_read_length,
                        pattern_to_chop_from_names = pattern_to_chop_from_names,
                        output_dir = output_dir,
                        base_for_names = base_for_names,
                        intermediate_reports = intermediate_reports,
                        erase_history = erase_history,
                        verbosity = verbosity,
                        report_type = report_type,
                        ncpu = ncpu,
                        bins_to_process = bins_to_process)
  } else {
    buildConfig_nol_test(fwd_file = fwd_file, 
                         fwd_primer_seq = fwd_primer_seq, 
                         fwd_primer_lens = fwd_primer_lens, 
                         fwd_min_score = fwd_min_score,
                         rev_file = rev_file, 
                         rev_primer_seq = rev_primer_seq, 
                         rev_primer_lens = rev_primer_lens, 
                         rev_min_score = rev_min_score,
                         fwd_pid_in_which_fragment = fwd_pid_in_which_fragment, 
                         rev_pid_in_which_fragment = rev_pid_in_which_fragment,
                         min_read_length = min_read_length,
                         pattern_to_chop_from_names = pattern_to_chop_from_names,
                         output_dir = output_dir,
                         base_for_names = base_for_names,
                         intermediate_reports = intermediate_reports,
                         erase_history = erase_history,
                         verbosity = verbosity,
                         report_type = report_type,
                         ncpu = ncpu,
                         bins_to_process = bins_to_process)
  }
}

#' Builds the config for overlapping reads - stable version
#' @export

buildConfig_ol_prod <- function(fwd_file, fwd_primer_seq, fwd_primer_lens, fwd_min_score,
                                rev_file, rev_primer_seq, rev_primer_lens, rev_min_score,
                                fwd_pid_in_which_fragment, rev_pid_in_which_fragment,
                                min_read_length = 295,
                                pattern_to_chop_from_names = ' [0-9]:N:[0-9]*:[0-9]*$',
                                output_dir = "/fridge/data/MotifBinner2_test",
                                base_for_names = "CAP129_2040_009wpi_C2C3",
                                intermediate_reports = TRUE,
                                erase_history = FALSE,
                                verbosity = 3,
                                report_type = c('html'),
                                ncpu = 4,
                                bins_to_process = Inf
                                )
{
  if (fwd_pid_in_which_fragment == "NULL"){fwd_pid_in_which_fragment <- NULL}
  if (rev_pid_in_which_fragment == "NULL"){rev_pid_in_which_fragment <- NULL}
  operation_list = list(
    'n001' = 
      list(name = 'fwd_loadData',
        op = 'loadData',
        data_source = fwd_file,
        cache_data = TRUE),
    'n002' =
      list(name = 'fwd_basicQC',
        op = 'basicQC',
        data_source = "n001",
        cache_data = FALSE),
    'n003' =
      list(name = 'fwd_ambigSeqs',
        op = 'ambigSeqs',
        data_source = "n001",
        threshold = 0.02,
        cache_data = TRUE),
    'n004' =
      list(name = 'fwd_primerDimer',
        op = 'primerDimer',
        data_source = "n003",
        threshold = 80,
        cache_data = TRUE),
    'n005' =
      list(name = 'fwd_seqLength',
        op = 'seqLength',
        data_source = "n004",
        threshold = min_read_length,
        cache_data = TRUE),
    'n006' =
      list(name = 'fwd_qualTrim',
        op = 'qualTrim',
        data_source = "n005",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.15,
        cache_data = TRUE),
    'n007' =
      list(name = 'fwd_trimAffixes',
        op = 'trimAffixes',
        data_source = "n006",
        
        primer_seq = fwd_primer_seq,
        primer_lens = fwd_primer_lens,
        min_score = fwd_min_score,

        primer_location = 'front',
        front_gaps_allowed = 0,
        cache_data = TRUE),
    'n008' = 
      list(name = 'rev_loadData',
        op = 'loadData',
        data_source = rev_file,
        cache_data = TRUE),
    'n009' =
      list(name = 'rev_basicQC',
        op = 'basicQC',
        data_source = "n008",
        cache_data = FALSE),
    'n010' =
      list(name = 'rev_ambigSeqs',
        op = 'ambigSeqs',
        data_source = "n008",
        threshold = 0.02,
        cache_data = TRUE),
    'n011' =
      list(name = 'rev_primerDimer',
        op = 'primerDimer',
        data_source = "n010",
        threshold = 80,
        cache_data = TRUE),
    'n012' =
      list(name = 'rev_seqLength',
        op = 'seqLength',
        data_source = "n011",
        threshold = min_read_length,
        cache_data = TRUE),
    'n013' =
      list(name = 'rev_qualTrim',
        op = 'qualTrim',
        data_source = "n012",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.15,
        cache_data = TRUE),
    'n014' =
      list(name = 'rev_trimAffixes',
        op = 'trimAffixes',
        data_source = "n013",
        
        primer_seq = rev_primer_seq,
        primer_lens = rev_primer_lens,
        min_score = rev_min_score,
        
        primer_location = 'front',
        front_gaps_allowed = 0,
        cache_data = TRUE),
    'n015' =
      list(name = 'fwd_extractPIDs',
        op = 'extractPIDs',
        data_source = "n007",
        pid_in_which_fragment = fwd_pid_in_which_fragment,
        pattern_to_chop_from_names = pattern_to_chop_from_names,
        pid_gaps_allowed = 0,
        cache_data = TRUE),
    'n016' =
      list(name = 'rev_extractPIDs',
        op = 'extractPIDs',
        data_source = "n014",
        pid_in_which_fragment = rev_pid_in_which_fragment,
        pattern_to_chop_from_names = pattern_to_chop_from_names,
        pid_gaps_allowed = 0,
        cache_data = TRUE),
    'n017' =
      list(name = 'matchPairs',
        op = 'matchPairs',
        data_source = c("fwd" = "n015", "rev" = "n016"),
        cache_data = TRUE),
    'n018' =
      list(name = 'processBadPIDs',
        op = 'processBadPIDs',
        data_source = "n017",
        cache_data = TRUE),
    'n019' =
      list(name = 'mergePEAR',
        op = 'mergePEAR',
        data_source = "n018",
        cache_data = TRUE),
    'n020' =
      list(name = 'merge_qualTrim',
        op = 'qualTrim',
        data_source = "n019",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.05,
        cache_data = TRUE),
    'n021' =
      list(name = 'binSizeCheck',
        op = 'binSizeCheck',
        data_source = "n020",
        min_bin_size = 3,
        cache_data = TRUE),
    'n022' =
      list(name = 'alignBinsMSA',
        op = 'alignBinsMSA',
        bins_to_process = bins_to_process,
        data_source = "n021",
        cache_data = TRUE),
    'n023' =
      list(name = 'buildConsensus',
        op = 'buildConsensus',
        data_source = "n022",
        cache_data = TRUE),
    'n024' =
      list(name = 'primerSeqErr',
        op = 'primerSeqErr',
        data_source = c("fwd" = "n007", "rev" = "n014"),
        cache_data = FALSE),
    'n025' =
      list(name = 'binSeqErr',
        op = 'binSeqErr',
        data_source = c("bin_msa_merged" = "n022", "cons_merged" = "n023", "primer_err" = "n024"),
        cache_data = FALSE),
    'n026' = 
      list(name = 'removeGaps',
        op = 'removeChars',
        data_source = "n023",
        char_to_remove = "-",
        cache_data = TRUE
           ),
    'n027' =
      list(name = 'vsearchCluster',
        op = 'vsearchCluster98',
        data_source = "n026",
        id = 0.98,
        min_clus_size = 1,
        cache_data = TRUE),
    'n028' =
      list(name = 'vsearchCluster',
        op = 'vsearchCluster95',
        data_source = "n026",
        id = 0.95,
        min_clus_size = 1,
        cache_data = TRUE),
    'n100' =
      list(name = 'dataTracing',
        op = 'dataTracing',
        data_source = c(
          "fwdReads.1" = "n001", "fwdReads.2" = "n003", "fwdReads.3" = "n004",
          "fwdReads.4" = "n005", "fwdReads.5" = "n006", "fwdReads.6" = "n007",

          "revReads.1" = "n008", "revReads.2" = "n010", "revReads.3" = "n011",
          "revReads.4" = "n012", "revReads.5" = "n013", "revReads.6" = "n014",

          "mergeReads.1" = "n017", "mergeReads.2" = "n018", "mergeReads.3" = "n019",
          "mergeReads.4" = "n020", "mergeReads.5" = "n021", "mergeReads.6" = "n022",
          "mergeReads.7" = "n023", "mergeReads.8" = "n027", "mergeReads.9" = "n028"),
      cache_data = FALSE)
  )

  return(list(operation_list = operation_list,
              output_dir = output_dir,
              base_for_names = base_for_names,
              intermediate_reports = intermediate_reports,
              verbosity = verbosity,
              erase_history = erase_history,
              report_type = report_type,
              ncpu = ncpu))
}

#' Test config builder for overlapping reads
#' @export

buildConfig_ol_test <- function(fwd_file, fwd_primer_seq, fwd_primer_lens, fwd_min_score,
                        rev_file, rev_primer_seq, rev_primer_lens, rev_min_score,
                        fwd_pid_in_which_fragment, rev_pid_in_which_fragment,
                        min_read_length = 295,
                        pattern_to_chop_from_names = ' [0-9]:N:[0-9]*:[0-9]*$',
                        output_dir = "/fridge/data/MotifBinner2_test",
                        base_for_names = "CAP129_2040_009wpi_C2C3",
                        intermediate_reports = TRUE,
                        erase_history = TRUE,
                        verbosity = 3,
                        report_type = c('html'),
                        ncpu = 4,
                        bins_to_process = Inf
                        )
{
  if (!is.null(fwd_pid_in_which_fragment)){
    if (fwd_pid_in_which_fragment == "NULL"){fwd_pid_in_which_fragment <- NULL}
  }
  if (!is.null(fwd_pid_in_which_fragment)){
    if (rev_pid_in_which_fragment == "NULL"){rev_pid_in_which_fragment <- NULL}
  }
  operation_list = list(
    'n001' = 
      list(name = 'fwd_loadData',
        op = 'loadData',
        data_source = fwd_file,
        cache_data = TRUE),
    'n002' =
      list(name = 'fwd_basicQC',
        op = 'basicQC',
        data_source = "n001",
        cache_data = FALSE),
    'n003' =
      list(name = 'fwd_ambigSeqs',
        op = 'ambigSeqs',
        data_source = "n001",
        threshold = 0.02,
        cache_data = TRUE),
    'n004' =
      list(name = 'fwd_primerDimer',
        op = 'primerDimer',
        data_source = "n003",
        threshold = 80,
        cache_data = TRUE),
    'n005' =
      list(name = 'fwd_seqLength',
        op = 'seqLength',
        data_source = "n004",
        threshold = min_read_length,
        cache_data = TRUE),
    'n006' =
      list(name = 'fwd_qualTrim',
        op = 'qualTrim',
        data_source = "n005",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.10,
        cache_data = TRUE),
    'n007' =
      list(name = 'fwd_trimAffixes',
        op = 'trimAffixes',
        data_source = "n006",
        
        primer_seq = fwd_primer_seq,
        primer_lens = fwd_primer_lens,
        min_score = fwd_min_score,

        primer_location = 'front',
        front_gaps_allowed = 0,
        cache_data = TRUE),
    'n008' = 
      list(name = 'rev_loadData',
        op = 'loadData',
        data_source = rev_file,
        cache_data = TRUE),
    'n009' =
      list(name = 'rev_basicQC',
        op = 'basicQC',
        data_source = "n008",
        cache_data = FALSE),
    'n010' =
      list(name = 'rev_ambigSeqs',
        op = 'ambigSeqs',
        data_source = "n008",
        threshold = 0.02,
        cache_data = TRUE),
    'n011' =
      list(name = 'rev_primerDimer',
        op = 'primerDimer',
        data_source = "n010",
        threshold = 80,
        cache_data = TRUE),
    'n012' =
      list(name = 'rev_seqLength',
        op = 'seqLength',
        data_source = "n011",
        threshold = min_read_length,
        cache_data = TRUE),
    'n013' =
      list(name = 'rev_qualTrim',
        op = 'qualTrim',
        data_source = "n012",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.20,
        cache_data = TRUE),
    'n014' =
      list(name = 'rev_trimAffixes',
        op = 'trimAffixes',
        data_source = "n013",
        
        primer_seq = rev_primer_seq,
        primer_lens = rev_primer_lens,
        min_score = rev_min_score,
        
        primer_location = 'front',
        front_gaps_allowed = 0,
        cache_data = TRUE),
    'n015' =
      list(name = 'fwd_extractPIDs',
        op = 'extractPIDs',
        data_source = "n007",
        pid_in_which_fragment = fwd_pid_in_which_fragment,
        pattern_to_chop_from_names = pattern_to_chop_from_names,
        pid_gaps_allowed = 0,
        cache_data = TRUE),
    'n016' =
      list(name = 'rev_extractPIDs',
        op = 'extractPIDs',
        data_source = "n014",
        pid_in_which_fragment = rev_pid_in_which_fragment,
        pattern_to_chop_from_names = pattern_to_chop_from_names,
        pid_gaps_allowed = 0,
        cache_data = TRUE),
    'n017' =
      list(name = 'matchPairs',
        op = 'matchPairs',
        data_source = c("fwd" = "n015", "rev" = "n016"),
        cache_data = TRUE),
    'n018' =
      list(name = 'processBadPIDs',
        op = 'processBadPIDs',
        data_source = "n017",
        cache_data = TRUE),
    'n019' =
      list(name = 'mergePEAR',
        op = 'mergePEAR',
        data_source = "n018",
        cache_data = TRUE),
    'n020' =
      list(name = 'merge_qualTrim',
        op = 'qualTrim',
        data_source = "n019",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.05,
        cache_data = TRUE),
    'n021' =
      list(name = 'binSizeCheck',
        op = 'binSizeCheck',
        data_source = "n020",
        min_bin_size = 3,
        cache_data = TRUE),
    'n022' =
      list(name = 'alignBinsMSA',
        op = 'alignBinsMSA',
        bins_to_process = bins_to_process,
        data_source = "n021",
        cache_data = TRUE),
    'n023' =
      list(name = 'buildConsensus',
        op = 'buildConsensus',
        data_source = "n022",
        cache_data = TRUE),
    'n024' =
      list(name = 'primerSeqErr',
        op = 'primerSeqErr',
        data_source = c("fwd" = "n007", "rev" = "n014"),
        cache_data = FALSE),
    'n025' =
      list(name = 'binSeqErr',
        op = 'binSeqErr',
        data_source = c("bin_msa_merged" = "n022", "cons_merged" = "n023", "primer_err" = "n024"),
        cache_data = FALSE),
    'n026' = 
      list(name = 'removeGaps',
        op = 'removeChars',
        data_source = "n023",
        char_to_remove = "-",
        cache_data = TRUE
           ),
    'n040' =
      list(name = 'fwd_extractReads',
        op = 'extractData',
        data_source = "n018",
        extract_levels = c("seq_dat", "fwd"),
        cache_data = TRUE),
    'n100' =
      list(name = 'dataTracing',
        op = 'dataTracing',
        data_source = c(
          "fwdReads.1" = "n001", "fwdReads.2" = "n003", "fwdReads.3" = "n004",
          "fwdReads.4" = "n005", "fwdReads.5" = "n006", "fwdReads.6" = "n007",

          "revReads.1" = "n008", "revReads.2" = "n010", "revReads.3" = "n011",
          "revReads.4" = "n012", "revReads.5" = "n013", "revReads.6" = "n014",

          "mergeReads.1" = "n017", "mergeReads.2" = "n018", "mergeReads.3" = "n019",
          "mergeReads.4" = "n020", "mergeReads.5" = "n021", "mergeReads.6" = "n022",
          "mergeReads.7" = "n023"),
      cache_data = FALSE)
    )

  config <- list(operation_list = operation_list,
                 output_dir = output_dir,
                 base_for_names = base_for_names,
                 intermediate_reports = intermediate_reports,
                 verbosity = verbosity,
                 erase_history = erase_history,
                 report_type = report_type,
                 ncpu = ncpu)
}

#' Test config builder for non-overlapping reads
#' @export

buildConfig_nol_test <- function(fwd_file, fwd_primer_seq, fwd_primer_lens, fwd_min_score,
                        rev_file, rev_primer_seq, rev_primer_lens, rev_min_score,
                        fwd_pid_in_which_fragment, rev_pid_in_which_fragment,
                        min_read_length = 295,
                        pattern_to_chop_from_names = ' [0-9]:N:[0-9]*:[0-9]*$',
                        output_dir = "/fridge/data/MotifBinner2_test",
                        base_for_names = "CAP129_2040_009wpi_C2C3",
                        intermediate_reports = TRUE,
                        erase_history = TRUE,
                        verbosity = 3,
                        report_type = c('html'),
                        ncpu = 4,
                        bins_to_process = Inf
                        )
{
  if (!is.null(fwd_pid_in_which_fragment)){
    if (fwd_pid_in_which_fragment == "NULL"){fwd_pid_in_which_fragment <- NULL}
  }
  if (!is.null(fwd_pid_in_which_fragment)){
    if (rev_pid_in_which_fragment == "NULL"){rev_pid_in_which_fragment <- NULL}
  }
  operation_list = list(
    'n001' = 
      list(name = 'fwd_loadData',
        op = 'loadData',
        data_source = fwd_file,
        cache_data = TRUE),
    'n002' =
      list(name = 'fwd_basicQC',
        op = 'basicQC',
        data_source = "n001",
        cache_data = FALSE),
    'n003' =
      list(name = 'fwd_ambigSeqs',
        op = 'ambigSeqs',
        data_source = "n001",
        threshold = 0.02,
        cache_data = TRUE),
    'n004' =
      list(name = 'fwd_primerDimer',
        op = 'primerDimer',
        data_source = "n003",
        threshold = 80,
        cache_data = TRUE),
    'n005' =
      list(name = 'fwd_seqLength',
        op = 'seqLength',
        data_source = "n004",
        threshold = min_read_length,
        cache_data = TRUE),
    'n006' =
      list(name = 'fwd_qualTrim',
        op = 'qualTrim',
        data_source = "n005",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.15,
        cache_data = TRUE),
    'n007' =
      list(name = 'fwd_trimAffixes',
        op = 'trimAffixes',
        data_source = "n006",
        
        primer_seq = fwd_primer_seq,
        primer_lens = fwd_primer_lens,
        min_score = fwd_min_score,

        primer_location = 'front',
        front_gaps_allowed = 0,
        cache_data = TRUE),
    'n008' = 
      list(name = 'rev_loadData',
        op = 'loadData',
        data_source = rev_file,
        cache_data = TRUE),
    'n009' =
      list(name = 'rev_basicQC',
        op = 'basicQC',
        data_source = "n008",
        cache_data = FALSE),
    'n010' =
      list(name = 'rev_ambigSeqs',
        op = 'ambigSeqs',
        data_source = "n008",
        threshold = 0.02,
        cache_data = TRUE),
    'n011' =
      list(name = 'rev_primerDimer',
        op = 'primerDimer',
        data_source = "n010",
        threshold = 80,
        cache_data = TRUE),
    'n012' =
      list(name = 'rev_seqLength',
        op = 'seqLength',
        data_source = "n011",
        threshold = min_read_length,
        cache_data = TRUE),
    'n013' =
      list(name = 'rev_qualTrim',
        op = 'qualTrim',
        data_source = "n012",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.15,
        cache_data = TRUE),
    'n014' =
      list(name = 'rev_trimAffixes',
        op = 'trimAffixes',
        data_source = "n013",
        
        primer_seq = rev_primer_seq,
        primer_lens = rev_primer_lens,
        min_score = rev_min_score,
        
        primer_location = 'front',
        front_gaps_allowed = 0,
        cache_data = TRUE),
    'n015' =
      list(name = 'fwd_extractPIDs',
        op = 'extractPIDs',
        data_source = "n007",
        pid_in_which_fragment = fwd_pid_in_which_fragment,
        pattern_to_chop_from_names = pattern_to_chop_from_names,
        pid_gaps_allowed = 0,
        cache_data = TRUE),
    'n016' =
      list(name = 'rev_extractPIDs',
        op = 'extractPIDs',
        data_source = "n014",
        pid_in_which_fragment = rev_pid_in_which_fragment,
        pattern_to_chop_from_names = pattern_to_chop_from_names,
        pid_gaps_allowed = 0,
        cache_data = TRUE),
    'n017' =
      list(name = 'matchPairs',
        op = 'matchPairs',
        data_source = c("fwd" = "n015", "rev" = "n016"),
        cache_data = TRUE),
    'n018' =
      list(name = 'processBadPIDs',
        op = 'processBadPIDs',
        data_source = "n017",
        cache_data = TRUE),

    'n019' =
      list(name = 'fwd_extractReads',
        op = 'extractData',
        data_source = "n018",
        extract_levels = c("seq_dat", "fwd"),
        cache_data = TRUE),
    'n020' =
      list(name = 'fwd_alignBinsMSA',
        op = 'alignBinsMSA',
        bins_to_process = bins_to_process,
        data_source = "n019",
        cache_data = TRUE),
    'n021' =
      list(name = 'fwd_buildConsensus',
        op = 'buildConsensus',
        data_source = "n020",
        cache_data = TRUE),
    'n022' = 
      list(name = 'fwd_removeGaps',
        op = 'removeChars',
        data_source = "n021",
        char_to_remove = "-",
        cache_data = TRUE
           ),
    'n023' =
      list(name = 'rev_extractReads',
        op = 'extractData',
        data_source = "n018",
        extract_levels = c("seq_dat", "rev"),
        cache_data = TRUE),
    'n024' =
      list(name = 'rev_alignBinsMSA',
        op = 'alignBinsMSA',
        bins_to_process = bins_to_process,
        data_source = "n023",
        cache_data = TRUE),
    'n025' =
      list(name = 'rev_buildConsensus',
        op = 'buildConsensus',
        data_source = "n024",
        cache_data = TRUE),
    'n026' = 
      list(name = 'rev_removeGaps',
        op = 'removeChars',
        data_source = "n025",
        char_to_remove = "-",
        cache_data = TRUE
           ),
    'n030' =
      list(name = 'primerSeqErr',
        op = 'primerSeqErr',
        data_source = c("fwd" = "n007", "rev" = "n014"),
        cache_data = FALSE),
    'n031' =
      list(name = 'binSeqErr',
        op = 'binSeqErr',
        data_source = c("bin_msa_fwd" = "n020", "bin_msa_rev" = "n024",
                        "cons_fwd" = "n021", "cons_rev" = "n025",
                        "primer_err" = "n030"),
        cache_data = FALSE),

    # todo binSeqErr for non-overlapping reads

    'n100' =
      list(name = 'dataTracing',
        op = 'dataTracing',
        data_source = c(
          "fwdReads.01" = "n001", "fwdReads.02" = "n003", "fwdReads.03" = "n004",
          "fwdReads.04" = "n005", "fwdReads.05" = "n006", "fwdReads.06" = "n007",
          "fwdReads.07" = "n017", "fwdReads.08" = "n018", "fwdReads.09" = "n020",
          "fwdReads.10" = "n021",

          "revReads.01" = "n008", "revReads.02" = "n010", "revReads.03" = "n011",
          "revReads.04" = "n012", "revReads.05" = "n013", "revReads.06" = "n014",
          "revReads.07" = "n017", "revReads.08" = "n018", "revReads.09" = "n024",
          "revReads.10" = "n025"),
      cache_data = FALSE)
    )

  config <- list(operation_list = operation_list,
                 output_dir = output_dir,
                 base_for_names = base_for_names,
                 intermediate_reports = intermediate_reports,
                 verbosity = verbosity,
                 erase_history = erase_history,
                 report_type = report_type,
                 ncpu = ncpu)
}

