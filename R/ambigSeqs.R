#' Removes sequences with too many ambig bases
#' @inheritParams applyOperation
#' @export

ambigSeqs <- function(all_results, config)
{
  op_dir <- file.path(config$output_dir, config$prefix_for_names,
                      paste('n', sprintf("%03d", length(all_results)+1), '_ambigSeqs', sep = ''))
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)

  kept <- list()
  trimmed <- list()
  for (data_set_name in names(all_results[[length(all_results)]]$kept)){
    seq_dat <- all_results[[length(all_results)]]$kept[[data_set_name]]
    if (length(seq_dat) > 0)
    {
      tmp <- ambigSeqs_internal(seq_dat, config$ambigSeqs$max_ambig)
      kept[[data_set_name]] <- tmp$kept
      trimmed[[data_set_name]] <- tmp$trimmed
      rm(tmp)
    }
  }

  result <- list(kept = kept,
                 trimmed = trimmed,
                 step_num = length(all_results)+1,
                 op_dir = op_dir)
  class(result) <- 'ambigSeqs'
  return(result)
}

ambigSeqs_internal <- function(seq_dat, max_ambig)
{
  counts <- alphabetFrequency(seq_dat@sread)
  tmp <- data.frame(counts)
  names(tmp) <- paste('c',attr(counts, 'dimnames')[[2]],sep='')
  counts <- tmp
  rm(tmp)

  ambigCols <- !(gsub('^c','', names(counts)) %in% c('A','C','G','T','-'))
  counts$ambig <- apply(counts[,ambigCols], 1, sum)
  counts$perc_ambig <- counts$ambig/(apply(counts, 1, sum) - counts$ambig)
  if (max_ambig < 1)
  {
    kept_list <- counts$perc_ambig <= max_ambig
  } else {
    kept_list <- counts$ambig <= max_ambig
  }
  return(list(kept = seq_dat[kept_list],
              trimmed = seq_dat[(!kept_list)]))
}

saveToDisk.ambigSeqs <- function(result, config)
{
  for (data_set_name in names(result$kept)){
    seq_dat <- result$kept[[data_set_name]]
    if (length(seq_dat) > 0)
    {
      writeFastq(seq_dat, file.path(result$op_dir, 
        paste(config$prefix_for_names, '_kept_', data_set_name, '.fastq', sep = '')), compress=F)
    }
  }
  for (data_set_name in names(result$trimmed)){
    seq_dat <- result$trimmed[[data_set_name]]
    if (length(seq_dat) > 0)
    {
      writeFastq(seq_dat, file.path(result$op_dir, 
        paste(config$prefix_for_names, '_trimmed_', data_set_name, '.fastq', sep = '')), compress=F)
    }
  }

}

genSummary.ambigSeqs <- function(result, config)
{
  summary_tab <- rbind(
    genSummary_internal(operation = 'ambigSeqs',
                        parameters = 'fwd_reads',
                        kept_seq_dat = result$kept$fwd_reads,
                        trimmed_seq_dat = result$trimmed$fwd_reads),
    genSummary_internal(operation = 'ambigSeqs',
                        parameters = 'rev_reads',
                        kept_seq_dat = result$kept$rev_reads,
                        trimmed_seq_dat = result$trimmed$rev_reads))
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
