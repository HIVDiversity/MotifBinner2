#' Loads the fwd and rev fastq reads
#' @inheritParams applyOperation
#' @export

loadData <- function(all_results, config)
{
  op_name <- paste('n', sprintf("%03d", length(all_results)+1), '_loadData', sep = '')
  op_dir <- file.path(config$output_dir, config$prefix_for_names, op_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(config$fwd_reads_file))
  {
    fwd_reads <- readFastq(config$fwd_reads_file)
  }
  if (!is.null(config$rev_reads_file))
  {
    rev_reads <- readFastq(config$rev_reads_file)
  }
  result <- list(
    fwd_reads = list(
      kept = list(
        seq_dat = fwd_reads
                  )
                     ),
    rev_reads = list(
      kept = list(
        seq_dat = rev_reads)
      ),
    step_num = length(all_results)+1,
    op_dir = op_dir,
    op_name = op_name
                 )
  class(result) <- 'loadData'
  return(result)
}

saveToDisk.loadData <- function(result, config)
{
  return(result)
}

genSummary.loadData <- function(result, config)
{
  summary_tabs <- list()
  for (orient in c('fwd_reads', 'rev_reads')){
    if (orient %in% names(result)){
      for (status in c('kept', 'trimmed')){
        if (status %in% names(result[[orient]])){
          status_dat <- result[[orient]][[status]]
          result[[orient]][[status]][['summary']] <- genSummary_internal_one(result[[orient]][[status]][['seq_dat']])
          summary_tabs[[paste(orient, status, sep = '_')]] <- result[[orient]][[status]][['summary']]
        }
      }
    }
  }

  for (i in c('fwd_reads_kept', 'fwd_reads_trimmed',
              'rev_reads_kept', 'rev_reads_trimmed'))
  {

  }

  summary_tab <- rbind(
    genSummary_internal(operation = class(result),
                        parameters = 'fwd_reads',
                        kept_seq_dat = result$kept$fwd_reads,
                        trimmed_seq_dat = DNAStringSet(NULL)),
    genSummary_internal(operation = class(result),
                        parameters = 'rev_reads',
                        kept_seq_dat = result$kept$rev_reads,
                        trimmed_seq_dat = DNAStringSet(NULL)))
  result$summary <- summary_tab
  write.csv(summary_tab, file.path(result$op_dir, 'loadData_summary.csv'), row.names=FALSE)
  return(result)
}

genSummary_internal <- function(operation, parameters, kept_seq_dat, trimmed_seq_dat)
{
  kept_tab <- genSummary_internal_one(kept_seq_dat)
  names(kept_tab) <- paste(names(kept_tab), 'kept', sep = '_')
  trimmed_tab <- genSummary_internal_one(trimmed_seq_dat)
  names(trimmed_tab) <- paste(names(trimmed_tab), 'trimmed', sep = '_')
  cbind(data.frame(operation = operation,
                   parameters = parameters),
        kept_tab, trimmed_tab)
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
