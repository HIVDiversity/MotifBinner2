#' Loads the fwd and rev fastq reads
#' @inheritParams applyOperation
#' @export

loadData <- function(all_results, config)
{
  dir.create(file.path(config$output_dir, config$prefix_for_names, 
                       paste('n', sprintf("%03d", length(all_results)), 'loadData', 
                             sep = '')),
             showWarnings = FALSE, recursive = TRUE)
  if (!is.null(config$fwd_reads_file))
  {
    fwd_reads <- readFastq(config$fwd_reads_file)
  }
  if (!is.null(config$rev_reads_file))
  {
    rev_reads <- readFastq(config$rev_reads_file)
  }
  final <- list(fwd_reads = fwd_reads,
                rev_reads = rev_reads)
  result <- list(final = final,
                 step_num = length(all_results))
  class(result) <- 'loadData'
  return(result)
}

saveToDisk.loadData <- function(result, config)
{
  return(NULL)
}

genSummary.loadData <- function(result, config)
{
  summary_tab <- rbind(
    genSummary_internal(operation = 'loadData',
                        parameters = 'fwd_reads',
                        good_seq_dat = result$final$fwd_reads,
                        bad_seq_dat = DNAStringSet(NULL)),
    genSummary_internal(operation = 'loadData',
                        parameters = 'rev_reads',
                        good_seq_dat = result$final$rev_reads,
                        bad_seq_dat = DNAStringSet(NULL)))
  result$summary <- summary_tab
  return(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
