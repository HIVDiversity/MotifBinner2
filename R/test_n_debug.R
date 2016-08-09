buildConfig <- function(fwd_file, fwd_primer_seq, fwd_primer_lens, fwd_min_score,
                        rev_file, rev_primer_seq, rev_primer_lens, rev_min_score,
                        fwd_pid_in_which_fragment, rev_pid_in_which_fragment,
                        profile_file,
                        pattern_to_chop_from_names = ' [0-9]:N:[0-9]*:[0-9]*$',
                        output_dir = "/fridge/data/MotifBinner2_test",
                        base_for_names = "CAP129_2040_009wpi_C2C3",
                        intermediate_reports = TRUE,
                        erase_history = TRUE,
                        verbosity = 3,
                        report_type = c('html'),
                        ncpu = 4
                        )
{
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
        threshold = 295,
        cache_data = TRUE),
    'n006' =
      list(name = 'fwd_qualTrim',
        op = 'qualTrim',
        data_source = "n005",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.05,
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
        threshold = 295,
        cache_data = TRUE),
    'n013' =
      list(name = 'rev_qualTrim',
        op = 'qualTrim',
        data_source = "n012",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.05,
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
      list(name = 'alignBins',
        op = 'alignBins',
        bins_to_process = Inf,
        data_source = "n018",
        profile_file = profile_file,
        cache_data = TRUE),
    'n020' =
      list(name = 'buildConsensus',
        op = 'buildConsensus',
        data_source = "n019",
        cache_data = TRUE),
    'n021' =
      list(name = 'primerSeqErr',
        op = 'primerSeqErr',
        data_source = c("fwd" = "n007", "rev" = "n014"),
        cache_data = FALSE),
    'n022' =
      list(name = 'binSeqErr',
        op = 'binSeqErr',
        data_source = c("bin_msa" = "n019", "cons" = "n020", "primer_err" = "n021"),
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

store_configs <- function()
{
#buildConfig <- function(fwd_file, fwd_primer_seq, fwd_primer_lens, fwd_min_score,
#                        rev_file, rev_primer_seq, rev_primer_lens, rev_min_score,
#                        fwd_pid_in_which_fragment, rev_pid_in_which_fragment,
#                        profile_file, output_dir, base_for_names,
#                        intermediate_reports = TRUE,
#                        verbosity = 3,
#                        report_type = c('html')
#                        )
#  buildConfig(fwd_file = 
#              fwd_primer_seq = 
#              fwd_primer_lens = 
#              fwd_min_score = 
#              rev_file = 
#              rev_primer_seq = 
#              rev_primer_lens = 
#              rev_min_score = 
#              fwd_pid_in_which_fragment = 
#              rev_pid_in_which_fragment = 
#              profile_file = 
#              output_dir = 
#              base_for_names = 
#              )
  config <-
  buildConfig(fwd_file = "/fridge/data/MotifBinner2_test/raw/CAP129_2040_009wpi_C2C3_R1.fastq",
              fwd_primer_seq = 'CTCTTTTGACCCAATTCCTATACATTATTG',
              fwd_primer_lens = 30,
              fwd_min_score = 20,
              rev_file = "/fridge/data/MotifBinner2_test/raw/CAP129_2040_009wpi_C2C3_R2.fastq",
              rev_primer_seq = 'NNNNNNNNNNNNNNTGCAATAGAAAAATTCTCCTCTACAATT',
              rev_primer_lens = c(14, 28),
              rev_min_score = 28,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              profile_file = "/fridge/data/MotifBinner2_test/c2c3_mol_clock_profile_1.fasta",
              output_dir = "/fridge/data/MotifBinner2_test",
              base_for_names = "CAP129_2040_009wpi_C2C3"
              )
  config <-
  buildConfig(fwd_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq",
              fwd_primer_seq = 'TATGGGAYSAAAGYCTMAARCCATGTG',
              fwd_primer_lens = 27,
              fwd_min_score = 18,
              rev_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R2.fastq",
              rev_primer_seq = 'CACACGCTCAGNNNNNNNNNATTCCATGTGTACATTGTACTGTRCTG',
              rev_primer_lens = c(11, 9, 27),
              rev_min_score = 30,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 2,
              profile_file = "/fridge/data/MotifBinner2_test/v1v2_profile1.fasta",
              output_dir = "/fridge/data/MotifBinner2_test",
              base_for_names = "CAP256_3100_030wpi_v1v2"
              )
  config <-    ## Zhou et al 2015 env region
  buildConfig(fwd_file = "/fridge/data/zhou2015/data/SRR1761912/SRR1761912_1.fastq",
              fwd_primer_seq = "NNNNTTATGGGATCAAAGCCTAAAGCCATGTGTA",
              fwd_primer_lens = 34,
              fwd_min_score = 23,
              rev_file = "/fridge/data/zhou2015/data/SRR1761912/SRR1761912_2.fastq",
              rev_primer_seq = "NNNNNNNNNCAGTCCATTTTGCTCTACTAATGTTACAATGTGC",
              rev_primer_lens = c(9, 34),
              rev_min_score = 29,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              profile_file = "/fridge/data/zhou2015/NA_aligned.fasta",
              pattern_to_chop_from_names = "\\.[12] [0-9]+ length=[0-9].*",
              output_dir = "/fridge/data/zhou2015/binned",
              base_for_names = "zhou_2015_v1v3_1761912",
              erase_history = FALSE,
              ncpu = 2
              )
}

dummy_test_debug <- function()
{
  getwd()
  library(devtools)
  setwd('~/projects/MotifBinner2/code/MotifBinner2')
  load_all(quiet=TRUE)
  unlink(file.path(config$output_dir, config$base_for_names), recursive=T)
#  load_all()
  all_results <- list()
  class(all_results) <- 'allResults'

  all_results <- applyOperation(all_results, config, op_number = 'n001') # fwd_loadData
  all_results <- applyOperation(all_results, config, op_number = 'n002') # fwd_basicQC
  all_results <- applyOperation(all_results, config, op_number = 'n003') # fwd_ambigSeqs
  all_results <- applyOperation(all_results, config, op_number = 'n004') # fwd_primerDimer
  all_results <- applyOperation(all_results, config, op_number = 'n005') # fwd_seqLength
  all_results <- applyOperation(all_results, config, op_number = 'n006') # fwd_qualTrim
  all_results <- applyOperation(all_results, config, op_number = 'n007') # fwd_trimAffixes
  all_results <- applyOperation(all_results, config, op_number = 'n008') # rev_loadData
  all_results <- applyOperation(all_results, config, op_number = 'n009') # rev_basicQC
  all_results <- applyOperation(all_results, config, op_number = 'n010') # rev_ambigSeqs
  all_results <- applyOperation(all_results, config, op_number = 'n011') # rev_primerDimer
  all_results <- applyOperation(all_results, config, op_number = 'n012') # rev_seqLength
  all_results <- applyOperation(all_results, config, op_number = 'n013') # rev_qualTrim
  all_results <- applyOperation(all_results, config, op_number = 'n014') # rev_trimAffixes
  all_results <- applyOperation(all_results, config, op_number = 'n015') # fwdExtractPIDs
  all_results <- applyOperation(all_results, config, op_number = 'n016') # revExtractPIDs
  all_results <- applyOperation(all_results, config, op_number = 'n017') # matchPairs
  all_results <- applyOperation(all_results, config, op_number = 'n018') # processBadPIDs
  all_results <- applyOperation(all_results, config, op_number = 'n019') # alignBins
  all_results <- applyOperation(all_results, config, op_number = 'n020') # buildConsensus
  all_results <- applyOperation(all_results, config, op_number = 'n021') # primerSeqErr
  all_results <- applyOperation(all_results, config, op_number = 'n022') # binSeqErr

  timing <- list()
  ptm <- proc.time()
  op_number <- 'n022'
  config$current_op_number <- op_number
  op <- get(config$operation_list[[op_number]]$op)
  result <- op(all_results, config)
  seq_dat <- result$input_dat
  result$input_dat <- NULL
  timing$main <- proc.time() - ptm
  
  ptm <- proc.time()
  result <- genSummary(result, config, seq_dat = seq_dat)
  timing$summary <- proc.time() - ptm
  
  ptm <- proc.time()
  result <- computeMetrics(result, config, seq_dat = seq_dat)
  timing$metrics <- proc.time() - ptm
  
  ptm <- proc.time()
  result$timing <- timing
  result <- genReport(result, config)
  timing$report <- proc.time() - ptm
  
  ptm <- proc.time()
  result <- saveToDisk(result, config, seq_dat)
  timing$saveToDisk <- proc.time() - ptm

  ptm <- proc.time()

  print('bye')

#  config <- list(fwd_reads_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_20k_R1.fastq",
#                 rev_reads_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_20k_R2.fastq",
#                 output_dir = "/fridge/data/MotifBinner2_test",
#                 prefix_for_names = "",
#                 fwd_primer = 'TATGGGAYSAAAGYCTMAARCCATGTG',
#                 rev_primer = 'CACACGCTCAGNNNNNNNNNATTCCATGTGTACATTGTACTGTRCTG',
#                 operation_list = c('loadData', 'basicQC', 'ambigSeqs', 'primerDimer', 'seqLength'),
#                 ambigSeqs = list(max_ambig = 5),
#                 primerDimer = list(primer_dimer_len = 80),
#                 seqLength = list(short_seq_len = 295,
#                                  long_seq_len = 305),
#                 intermediate_reports = TRUE,
#                 verbosity = 3,
#                 report_type = c('html'))
#
#  x <- do.call(processPrimers, config)
#  all_results <- x
#  all_results$summary <- NULL
#  
#  operation_function <- trimAffixes
#  config$operation_number <- length(all_results)
#  result <- operation_function(all_results, config)
#
#  str(result$fwd_result$kept$trim_stats)
#  print(result$fwd_result$kept$seq_dat@sread)
#
#  all_results <- list()
#  class(all_results) <- 'allResults'
#
#  dir.create(file.path(config$output_dir, config$prefix_for_names),
#             showWarnings = FALSE, recursive = TRUE)
#
#  all_results <- applyOperation('loadData', all_results, config)
#  all_results <- applyOperation('basicQC', all_results, config)
#  all_results <- applyOperation('ambigSeqs', all_results, config)
#  all_results <- applyOperation('primerDimer', all_results, config)
#  all_results <- applyOperation('seqLength', all_results, config)
#  all_results <- applyOperation('trimEnds', all_results, config)
#
#  setwd('~/projects/MotifBinner2/code/MotifBinner2')
#  unlink('/fridge/data/MotifBinner2_test/CAP256_3100_030wpi_v1v2_20k/n005_seqLength', recursive=T)
#  load_all(quiet=TRUE)
#
#  ptm <- proc.time()
#  timing <- list()
#  operation_function <- trimAffixes
#  config$operation_number <- length(all_results)
#
#  result <- operation_function(all_results, config)
#  timing$main <- proc.time() - ptm
#  ptm <- proc.time()
#
#  result <- saveToDisk(result, config)
#  timing$saveToDisk <- proc.time() - ptm
#
#  ptm <- proc.time()
#
#  result <- genSummary(result, config)
#  timing$genSummary <- proc.time() - ptm
#  ptm <- proc.time()
#
#  result <- computeMetrics(result, config)
#  timing$computeMetrics <- proc.time() - ptm
#  result$timing <- timing
#  ptm <- proc.time()
#
#  result <- genReport(result, config)
#  timing$genReport <- proc.time() - ptm
#  ptm <- proc.time()
#
#  result <- print(result, config)
#  timing$print <- proc.time() - ptm
#  result$timing <- timing
#
#  all_results[[basename(result$op_dir)]] <- result
#
#  return(all_results)
#
}

alignment_tests <- function()
{
  library(devtools)
  setwd('~/projects/MotifBinner2/code/MotifBinner2')
  Rcpp::compileAttributes()
#  load_all(quiet=TRUE)
  cat('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n'); load_all()
  seq_dat <- readFastq('~/projects/MotifBinner2/code/MotifBinner2/inst/test_dat.fastq')
  seq_dat@sread

  prefix <- 'TATGGGAYSAAAGYCTMAARCCATGTG'

  trimmed <- trimFront_cpp(as.character(seq_dat@sread),
                          as.character(seq_dat@quality@quality),
                          prefix, c(13, 7, 7))



  trimmed$sread
}

memory_optimization <- function()
{
  for (i in 1:length(all_results)){
    res <- all_results[[i]]
    cat(paste(names(all_results)[i], ':\n', sep = ''))
    for (j in names(res)){
      size_rep <- round(pryr::object_size(res[j])/1000000,0)
      if (size_rep > 1){
        cat(paste(' - ', j, ' : ', size_rep, '\n', sep = ''))
      }
    }
  }
}

build_dependency_map <- function(operation_list){
  dependencies <- list()
  for (i in 1:length(operation_list)){
    dependencies[[names(operation_list)[i]]] <- NULL
    for (j in operation_list[[i]]$data_source){
      if (grepl("^n[0-9]*$", j)){
        dependencies[[j]] <- c(dependencies[[j]], names(operation_list)[i])
      }
    }
  }
  ops_performed <- gsub("_.*$", "", names(all_results))
  for (purge_step in names(dependencies)){
    if (all(dependencies[[purge_step]] %in% ops_performed)){
      all_results$seq_dat <- NULL
      all_results$trim_dat <- NULL
    } else {
    }
  }
}

purge_all_results <- function(all_results){
}
