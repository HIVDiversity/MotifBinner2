dummy_test_debug <- function()
{
  library(devtools)
  load_all()

  all_results <- list()
  class(all_results) <- 'all_results'

  config <- list(fwd_reads_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq",
                 rev_reads_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R2.fastq",
                 output_dir = "/fridge/data/MotifBinner2_test",
                 prefix_for_names = "CAP256_3100_030wpi_v1v2",
                 intermediate_reports = TRUE,
                 verbosity = 3)

  result <- loadData(all_results, config)
}
