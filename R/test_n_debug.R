dummy_test_debug <- function()
{
  getwd()
  library(devtools)
  setwd('~/projects/MotifBinner2/code/MotifBinner2')
  load_all()

  all_results <- list()
  class(all_results) <- 'all_results'
  
  dir.create(file.path(config$output_dir, config$prefix_for_names), 
             showWarnings = FALSE, recursive = TRUE)

  config <- list(fwd_reads_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq",
                 rev_reads_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R2.fastq",
                 output_dir = "/fridge/data/MotifBinner2_test",
                 prefix_for_names = "CAP256_3100_030wpi_v1v2",
                 intermediate_reports = TRUE,
                 verbosity = 3,
                           report_type = c('html', 'pdf'))

  result <- loadData(all_results, config)
  result <- genSummary(result, config)
  result <- genReport(result, config)
}
