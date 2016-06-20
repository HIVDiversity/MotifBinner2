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
        data_source = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_20k_R1.fastq",
        cache_data = TRUE),
    'n002' =
      list(name = 'fwd_basicQC',
        op = 'basicQC',
        data_source = "n001",
        cache_data = FALSE)
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

  all_results <- applyOperation(all_results, config, operation = 'loadData')

  timing <- list()
  ptm <- proc.time()
  op_number <- 'n002'
  config$current_op_number <- op_number
  op <- get(config$operation_list[[op_number]]$op)
  result <- op(all_results, config)
  if ('seq_dat' %in% names(result))
  {
    seq_dat <- result$seq_dat
  } else {
    seq_dat <- result$tmp
    result$tmp <- NULL
  }
  timing$main <- proc.time() - ptm
  
  ptm <- proc.time()
  result <- genSummary(result, config, seq_dat = seq_dat)
  timing$summary <- proc.time() - ptm
  
  ptm <- proc.time()
  result <- computeMetrics(result, config)
  timing$metrics <- proc.time() - ptm
  
  ptm <- proc.time()
  result$timing <- timing
  result <- genReport(result, config)
  timing$report <- proc.time() - ptm
  
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
