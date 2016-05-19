#' Loads the fwd and rev fastq reads
#' @inheritParams applyOperation
#' @export

loadData <- function(all_results, config)
{
  if (!is.null(config$fwd_reads_file))
  {
    fwd_reads <- readFastq(config$fwd_reads_file)
  }
  if (!is.null(config$rev_reads_file))
  {
    rev_reads <- readFastq(config$rev_reads_file)
  }
  result <- list(fwd_reads = fwd_reads,
                 rev_reads = rev_reads)
  class(result) <- 'loadData'
  return(result)
}

saveToDisk.loadData <- function(result, config)
{
  return(NULL)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
