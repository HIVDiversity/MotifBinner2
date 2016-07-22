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

  required_dominance <- .05
  
  per_read_metrics <- data.frame(read_name = as.character(seq_dat@id),
                                 stringsAsFactors = F)
  per_read_metrics$pid <- gsub("_(fwd)|_(rev)$", "", gsub("^.*_PID:" , "", per_read_metrics$read_name))
  per_read_metrics$clean_pid <- gsub("_" , "", per_read_metrics$pid)
  full_alphabet <- row.names(consensusMatrix(DNAStringSet('A')))

  all_consensuses <- NULL
  pid <- unique(per_read_metrics$clean_pid)[1]
  uniq_pids <- unique(per_read_metrics$clean_pid)
  tmp_x <- foreach(pid = uniq_pids, .combine = "c") %dopar% {
#  for (pid in unique(per_read_metrics$clean_pid)){
    bin_seq_indx <- which(per_read_metrics$clean_pid == pid)
    bin_seqs <- seq_dat[bin_seq_indx]
#    writeFastq(bin_seqs, '/tmp/bin.fastq', compress=F)

    qual_mat <- as(FastqQuality(quality(quality(bin_seqs))), 'matrix')
    tweaked_qual_mat <- gapQualityTweaker_cpp(as.character(bin_seqs@sread),
                                   qual_mat)
    x <- scoreAlignmentPositions_cpp(as.character(bin_seqs@sread),
                                   tweaked_qual_mat)
    gsub('-', '', buildConsensus_cpp(x$score_mat, required_dominance)$consensus)
  }

  fake_qualities <- NULL
  for (i in 1:length(tmp_x)){
    fake_qualities <- c(fake_qualities, paste(rep('G', nchar(tmp_x[i])), collapse = ''))
  }

  consensuses <-
  ShortReadQ(sread = DNAStringSet(tmp_x),
             quality = BStringSet(fake_qualities),
             id = BStringSet(uniq_pids))
  
  per_read_metrics <- data.frame('read_exists' = rep(1, length(consensuses)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'buildConsensus'
  if (op_args$cache){
    result$seq_dat <- consensuses
  }
  result$input_dat <- consensuses
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
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


