dummy_test_debug <- function()
{
  unlink('/fridge/data/MotifBinner2_test/CAP256_3100_030wpi_v1v2_20k', recursive=T)
  getwd()
  library(devtools)
  setwd('~/projects/MotifBinner2/code/MotifBinner2')
  #Rcpp::compileAttributes()
  load_all(quiet=TRUE)

  # arguments passed in:
  operation_list = list(
    'n001' = 
      list(name = 'fwd_loadData',
        op = 'loadData',
        #data_source = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_80k_R1.fastq",
        data_source = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq",
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
        primer_seq = 'TATGGGAYSAAAGYCTMAARCCATGTG',
        primer_lens = 27,
        primer_location = 'front',
        min_score = 18,
        front_gaps_allowed = 0,
        cache_data = TRUE),
    'n008' = 
      list(name = 'rev_loadData',
        op = 'loadData',
        #data_source = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_80k_R2.fastq",
        data_source = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R2.fastq",
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
        primer_seq = 'CACACGCTCAGNNNNNNNNNATTCCATGTGTACATTGTACTGTRCTG',
        primer_lens = c(11, 9, 27),
        primer_location = 'front',
        min_score = 30,
        front_gaps_allowed = 3,
        cache_data = TRUE),
    'n015' =
      list(name = 'fwd_extractPIDs',
        op = 'extractPIDs',
        data_source = "n007",
        pid_in_which_fragment = NULL,
        pid_gaps_allowed = 0,
        cache_data = TRUE),
    'n016' =
      list(name = 'rev_extractPIDs',
        op = 'extractPIDs',
        data_source = "n014",
        pid_in_which_fragment = 2,
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
        profile_file = "/fridge/data/MotifBinner2_test/v1v2_profile1.fasta",
        cache_data = TRUE),
    'n020' =
      list(name = 'buildConsensus',
        op = 'buildConsensus',
        data_source = "n019",
        cache_data = TRUE)
    )

  output_dir = "/fridge/data/MotifBinner2_test"
  base_for_names = "CAP256_3100_030wpi_v1v2_20k"
  intermediate_reports = TRUE
  verbosity = 3
  report_type = c('html')

  config <- list(operation_list = operation_list,
                 output_dir = output_dir,
                 base_for_names = base_for_names,
                 intermediate_reports = TRUE,
                 verbosity = 3,
                 report_type = c('html'))

  all_results <- list()
  class(all_results) <- 'allResults'

  all_results <- applyOperation(all_results, config, op_number = 'n001')
  all_results <- applyOperation(all_results, config, op_number = 'n002')
  all_results <- applyOperation(all_results, config, op_number = 'n003')
  all_results <- applyOperation(all_results, config, op_number = 'n004')
  all_results <- applyOperation(all_results, config, op_number = 'n005')
  all_results <- applyOperation(all_results, config, op_number = 'n006')
  all_results <- applyOperation(all_results, config, op_number = 'n007')
  all_results <- applyOperation(all_results, config, op_number = 'n008')
  all_results <- applyOperation(all_results, config, op_number = 'n009')
  all_results <- applyOperation(all_results, config, op_number = 'n010')
  all_results <- applyOperation(all_results, config, op_number = 'n011')
  all_results <- applyOperation(all_results, config, op_number = 'n012')
  all_results <- applyOperation(all_results, config, op_number = 'n013')
  all_results <- applyOperation(all_results, config, op_number = 'n014')
  all_results <- applyOperation(all_results, config, op_number = 'n015')
  all_results <- applyOperation(all_results, config, op_number = 'n016')
  all_results <- applyOperation(all_results, config, op_number = 'n017')
  all_results <- applyOperation(all_results, config, op_number = 'n018')
  all_results <- applyOperation(all_results, config, op_number = 'n019')
  all_results <- applyOperation(all_results, config, op_number = 'n020')

  timing <- list()
  ptm <- proc.time()
  op_number <- 'n020'
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
