#' Removes sequences with too many ambig bases
#' @inheritParams applyOperation
#' @export

ambigSeqs <- function(all_results, config)
{
  op_dir <- file.path(config$output_dir, config$prefix_for_names,
                      paste('n', sprintf("%03d", length(all_results)+1), '_ambigSeqs', sep = ''))
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)

  for (data_set_name in names(all_results[[length(all_results)]]$final)){
    print(data_set_name)
  }

  final <- list(fwd_reads = 1,
                rev_reads = 2)

  result <- list(final = final,
                 step_num = length(all_results)+1,
                 op_dir = op_dir)
  class(result) <- 'ambigSeqs'
  return(result)
}

ambigSeqs_internal <- function(seq_dat, max_ambig)
{
  counts <- alphabetFrequency(seq_dat)
  tmp <- data.frame(counts)
  names(tmp) <- paste('c',attr(counts, 'dimnames')[[2]],sep='')
  counts <- tmp
  rm(tmp)

  ambigCols <- !(gsub('^c','', names(counts)) %in% c('A','C','G','T','-'))
  counts$ambig <- apply(counts[,ambigCols], 1, sum)
  counts$perc_ambig <- counts$ambig/(apply(counts, 1, sum) - counts$ambig)
}

saveToDisk.ambigSeqs <- function(result, config)
{
  return(result)
}

genSummary.ambigSeqs <- function(result, config)
{
  summary_tab <- rbind(
    genSummary_internal(operation = 'ambigSeqs',
                        parameters = 'fwd_reads',
                        kept_seq_dat = result$final$fwd_reads,
                        trimmed_seq_dat = DNAStringSet(NULL)),
    genSummary_internal(operation = 'ambigSeqs',
                        parameters = 'rev_reads',
                        kept_seq_dat = result$final$rev_reads,
                        trimmed_seq_dat = DNAStringSet(NULL)))
  result$summary <- summary_tab
  write.csv(summary_tab, file.path(result$op_dir, 'ambigSeqs_summary.csv'), row.names=FALSE)
  return(result)
}

computeMetrics.ambigSeqs <- function(result, config)
{
  return(result)
}

print.ambigSeqs <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: ambigSeqs')
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
