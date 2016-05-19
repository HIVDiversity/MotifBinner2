#' Applies a single operation
#'
#' @param operation The name of the operation to apply. A function by this name
#' with the arguments all_results and config must exists. Additionally, the
#' functions saveToDisk, genReport, genSummary and print must be specialized
#' for objects of the class matching the name of the operation.
#' @param all_results A list of all results given the class 'all_results'.
#' @param config A list of configuration options
#' @export

applyOperation <- function(operation, all_results, config)
{
  operation_function <- get(operation)
  config$operation_number <- length(all_results)

  result <- operation_function(all_results, config)
  result <- saveToDisk(result, config)
  result <- genReport(result, config)
  result <- genSummary(result, config)
  result <- print(result, config)
  all_results[[operation]] <- result
  
  return(all_results)
}

#' Processes a dataset with Primer IDs
#'
#' Applies a series of operations to the input dataset that generates consensus
#' sequences from raw (or preprocessed) datasets containing Primer IDs.
#'
#' @param fwd_reads The name of the fastq file with the forward reads.
#' @param rev_reads The name of the fastq file with the reverse reads.
#' @param output_dir The directory in which the output must be produced.
#' @param operation_list The list of operations to apply to the input data.
#' @param intermediate_reports Should intermediate reports be produced after
#' each operations? (Useful for debugging)
#' @param verbosity The output level. 0 = no output. 1 = info about which
#' operation is currently running. 2 = info about progress of current
#' operation. 3 = super verbose debugging output.
#' @export

processPrimers <- function(fwd_reads = NULL, rev_reads = NULL, output_dir = NULL, 
                           operation_list = c('loadData', 'basicQualityPlots'),
                           intermediate_reports = TRUE, verbosity = 0)
{
## Process arguments
  if (is.null(fwd_reads) & is.null(rev_reads))
  {
    stop('At least one of fwd_reads or rev_reads must be non-null')
  }

## Setup data structures
  all_results <- list()
  class(all_results) <- 'all_results'

  config <- list(fwd_reads = fwd_reads,
                 rev_reads = rev_reads,
                 output_dir = output_dir,
                 intermediate_reports = intermediate_reports,
                 verbosity = verbosity)

## Perform operations
  all_results <- applyOperation('loadData', all_results, config)  
  all_results <- applyOperation('basicQualityPlots', all_results, config)  

## Finalize
  all_results <- saveToDisk(all_results, config)
  all_results <- genReport(all_results, config)
  all_results <- genSummary(all_results, config)
  all_results <- print(all_results, config)

  return(all_results)
}
