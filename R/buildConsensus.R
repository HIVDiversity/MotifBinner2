#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

buildConsensus <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  seq_dat <- all_results[[data_source_indx]]$seq_dat
  
  per_read_metrics <- data.frame(read_name = as.character(seq_dat@id),
                                 stringsAsFactors = F)
  per_read_metrics$pid <- gsub("_(fwd)|_(rev)$", "", gsub("^.*_PID:" , "", per_read_metrics$read_name))
  per_read_metrics$clean_pid <- gsub("_" , "", per_read_metrics$pid)

  pid <- unique(per_read_metrics$clean_pid)[1]
  for (pid in unique(per_read_metrics$clean_pid)){
    bin_seq_indx <- which(per_read_metrics$clean_pid == pid)
    bin_seqs <- seq_dat[bin_seq_indx]
#    writeFastq(bin_seqs, '/tmp/bin.fastq', compress=F)

    qual_mat <- as(FastqQuality(quality(quality(bin_seqs))), 'matrix')
    x <- score_alignment_positions(as.character(seq_dat@sread),
                                   qual_mat)
  }


  return(result)
}

saveToDisk.buildConsensus <- function(result, config, seq_dat)
{
  kept <- getKept(result, seq_dat)
  trimmed <- getTrimmed(seq_dat = seq_dat, kept_dat = kept)

  if (length(kept) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_kept_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(kept, tmp_name, compress=F)
  }
  if (length(trimmed) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_trimmed_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(trimmed, tmp_name, compress=F)
  }
  return(result)
}

computeMetrics.buildConsensus <- function(result, config, seq_dat)
{
  return(result)
}

print.buildConsensus <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: buildConsensus')
  cat('\n-------------------')
  cat('\nKept Sequences:\n')
  print(result$summary[,c('parameter', 'k_seqs', 'k_mean_length', 'k_mean_qual')])
  cat('\n-------------------')
  cat('\nTrimmed Sequences:\n')
  print(result$summary[,c('parameter', 't_seqs', 't_mean_length', 't_mean_qual')])
  invisible(result)
}


