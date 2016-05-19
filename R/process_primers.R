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

processPrimers <- function(fwd_reads, rev_reads, output_dir, intermediate_reports = TRUE, verbosity = 0)
{
## Process arguments

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
