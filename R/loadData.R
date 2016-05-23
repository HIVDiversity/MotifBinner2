#' Loads the fwd and rev fastq reads
#' @inheritParams applyOperation
#' @export

loadData <- function(all_results, config)
{
  opdir <- file.path(config$output_dir, config$prefix_for_names,
                      paste('n', sprintf("%03d", length(all_results)+1), '_loadData', sep = ''))
  dir.create(opdir, showWarnings = FALSE, recursive = TRUE)
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
                 step_num = length(all_results)+1,
                 opdir = opdir)
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
                        kept_seq_dat = result$final$fwd_reads,
                        trimmed_seq_dat = DNAStringSet(NULL)),
    genSummary_internal(operation = 'loadData',
                        parameters = 'rev_reads',
                        kept_seq_dat = result$final$rev_reads,
                        trimmed_seq_dat = DNAStringSet(NULL)))
  result$summary <- summary_tab
  write.csv(summary_tab, file.path(result$opdir, 'loadData_summary.csv'), row.names=FALSE)
  return(result)
}

print.loadData <- function(result, config)
{
  cat('\nloadData\n')
  print(result$summary)
  return(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
