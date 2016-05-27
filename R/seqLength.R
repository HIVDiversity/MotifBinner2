#' Removes perimer-dimer sequences (based on sequence length)
#' @inheritParams applyOperation
#' @export

seqLength <- function(all_results, config)
{
  op_dir <- file.path(config$output_dir, config$prefix_for_names,
                      paste('n', sprintf("%03d", length(all_results)+1), '_seqLength', sep = ''))
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)

  kept <- list()
  trimmed <- list()
  for (data_set_name in names(all_results[[length(all_results)]]$kept)){
    seq_dat <- all_results[[length(all_results)]]$kept[[data_set_name]]
    if (length(seq_dat) > 0)
    {
      kept_list <- width(seq_dat) > config$seqLength$short_seq_len &
                   width(seq_dat) < config$seqLength$long_seq_len
      kept[[data_set_name]] <- seq_dat[kept_list]
      trimmed[[data_set_name]] <- seq_dat[!kept_list]
    }
  }

  result <- list(kept = kept,
                 trimmed = trimmed,
                 step_num = length(all_results)+1,
                 op_dir = op_dir)
  class(result) <- 'seqLength'
  return(result)
}

saveToDisk.seqLength <- function(result, config)
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
  return(result)
}

genSummary.seqLength <- function(result, config)
{
  summary_tab <- rbind(
    genSummary_internal(operation = 'seqLength',
                        parameters = 'fwd_reads',
                        kept_seq_dat = result$kept$fwd_reads,
                        trimmed_seq_dat = result$trimmed$fwd_reads),
    genSummary_internal(operation = 'seqLength',
                        parameters = 'rev_reads',
                        kept_seq_dat = result$kept$rev_reads,
                        trimmed_seq_dat = result$trimmed$rev_reads))
  result$summary <- summary_tab
  write.csv(summary_tab, file.path(result$op_dir, 'seqLength_summary.csv'), row.names=FALSE)
  return(result)
}

computeMetrics.seqLength <- function(result, config)
{
  return(result)
}

print.seqLength <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: seqLength')
  cat('\n-------------------')
  cat('\nKept Sequences:\n')
  print(result$summary[,c('parameters', 'seqs_kept', 'mean_length_kept', 'mean_qual_kept')])
  cat('\n-------------------')
  cat('\nTrimmed Sequences:\n')
  print(result$summary[,c('parameters', 'seqs_trimmed', 'mean_length_trimmed', 'mean_qual_trimmed')])
  return(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
