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
#' @param fwd_reads_file The name of the fastq file with the forward reads.
#' @param rev_reads_file The name of the fastq file with the reverse reads.
#' @param output_dir The directory in which the output must be produced.
#' @param prefix_for_names The basename to use for naming output.
#' @param operation_list The list of operations to apply to the input data.
#' @param intermediate_reports Should intermediate reports be produced after
#' each operations? (Useful for debugging)
#' @param verbosity The output level. 0 = no output. 1 = info about which
#' operation is currently running. 2 = info about progress of current
#' operation. 3 = super verbose debugging output.
#' @param report_type vector of types of reports to procude. Valid options:
#' 'html', 'pdf'.
#' @export

processPrimers <- function(fwd_reads_file = NULL, rev_reads_file = NULL, 
                           output_dir = NULL, prefix_for_names = NULL,
                           operation_list = c('loadData'),
                           intermediate_reports = TRUE, verbosity = 0,
                           report_type = c('html', 'pdf'))
{
## Process arguments
  if (is.null(fwd_reads_file) & is.null(rev_reads_file))
  {
    stop('At least one of fwd_reads or rev_reads must be non-null')
  }

## Setup data structures
  all_results <- list()
  class(all_results) <- 'all_results'

  config <- list(fwd_reads_file = fwd_reads_file,
                 rev_reads_file = rev_reads_file,
                 output_dir = output_dir,
                 prefix_for_names = prefix_for_names,
                 intermediate_reports = intermediate_reports,
                 verbosity = verbosity,
                 report_type = report_type)

## Perform operations
  for (operation in operation_list)
  {
    all_results <- applyOperation(operation, all_results, config)
  }

## Finalize
  all_results <- saveToDisk(all_results, config)
  all_results <- genReport(all_results, config)
  all_results <- genSummary(all_results, config)
  all_results <- print(all_results, config)

  return(all_results)
}
