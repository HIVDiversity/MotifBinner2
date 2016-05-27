#' Loads the fwd and rev fastq reads
#' @inheritParams applyOperation
#' @export

loadData <- function(all_results, config)
{
  op_dir <- file.path(config$output_dir, config$prefix_for_names,
                      paste('n', sprintf("%03d", length(all_results)+1), '_loadData', sep = ''))
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(config$fwd_reads_file))
  {
    fwd_reads <- readFastq(config$fwd_reads_file)
  }
  if (!is.null(config$rev_reads_file))
  {
    rev_reads <- readFastq(config$rev_reads_file)
  }
  kept <- list(fwd_reads = fwd_reads,
                rev_reads = rev_reads)
  result <- list(kept = kept,
                 step_num = length(all_results)+1,
                 op_dir = op_dir)
  class(result) <- 'loadData'
  return(result)
}

saveToDisk.loadData <- function(result, config)
{
  return(result)
}

genSummary.loadData <- function(result, config)
{
  summary_tab <- rbind(
    genSummary_internal(operation = 'loadData',
                        parameters = 'fwd_reads',
                        kept_seq_dat = result$kept$fwd_reads,
                        trimmed_seq_dat = DNAStringSet(NULL)),
    genSummary_internal(operation = 'loadData',
                        parameters = 'rev_reads',
                        kept_seq_dat = result$kept$rev_reads,
                        trimmed_seq_dat = DNAStringSet(NULL)))
  result$summary <- summary_tab
  write.csv(summary_tab, file.path(result$op_dir, 'loadData_summary.csv'), row.names=FALSE)
  return(result)
}

computeMetrics.loadData <- function(result, config)
{
  return(result)
}

print.loadData <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: loadData')
  cat('\n-------------------')
  cat('\nLoaded Sequences:\n')
  print(result$summary[,c('parameters', 'seqs_kept', 'mean_length_kept', 'mean_qual_kept')])
  return(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
