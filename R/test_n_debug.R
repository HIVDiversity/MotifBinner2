dummy_test_debug <- function()
{
  getwd()
  library(devtools)
  setwd('~/projects/MotifBinner2/code/MotifBinner2')
  load_all(quiet=TRUE)

  config <- list(fwd_reads_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_20k_R1.fastq",
                 rev_reads_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_20k_R2.fastq",
                 output_dir = "/fridge/data/MotifBinner2_test",
                 prefix_for_names = "CAP256_3100_030wpi_v1v2_20k",
                 operation_list = c('loadData', 'basicQC', 'ambigSeqs'),
                 ambigSeqs = list(max_ambig = 5),
                 intermediate_reports = TRUE,
                 verbosity = 3,
                 report_type = c('html', 'pdf'))

#  x <- do.call(processPrimers, config)

  all_results <- list()
  class(all_results) <- 'allResults'
  
  dir.create(file.path(config$output_dir, config$prefix_for_names), 
             showWarnings = FALSE, recursive = TRUE)

  all_results <- applyOperation('loadData', all_results, config)
  all_results <- applyOperation('basicQC', all_results, config)

  ptm <- proc.time()
  timing <- list()
  operation_function <- ambigSeqs
  config$operation_number <- length(all_results)

  result <- operation_function(all_results, config)
  timing$main <- proc.time() - ptm
  ptm <- proc.time()

  result <- saveToDisk(result, config)
  timing$saveToDisk <- proc.time() - ptm
  ptm <- proc.time()

  result <- genSummary(result, config)
  timing$genSummary <- proc.time() - ptm
  ptm <- proc.time()

  result <- computeMetrics(result, config)
  timing$computeMetrics <- proc.time() - ptm
  result$timing <- timing
  ptm <- proc.time()

  result <- genReport(result, config)
  timing$genReport <- proc.time() - ptm
  ptm <- proc.time()

  result <- print(result, config)
  timing$print <- proc.time() - ptm
  result$timing <- timing

  all_results[[basename(result$op_dir)]] <- result
  
  return(all_results)

}
