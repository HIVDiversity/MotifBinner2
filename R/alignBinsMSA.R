#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

alignBinsMSA <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(op_dir, 'bins'), showWarnings = FALSE, recursive = TRUE)
  
  data_source_indx <- grep(op_args$data_source, names(all_results))

  stopifnot(length(data_source_indx) == 1)
  seq_dat <- all_results[[data_source_indx]]$seq_dat

  bins_to_process <- op_args$bins_to_process

  per_read_metrics <- data.frame(read_name = as.character(seq_dat@id),
                                 stringsAsFactors = F)
  per_read_metrics$pid <- gsub("^.*_PID:" , "", per_read_metrics$read_name)
  per_read_metrics$clean_pid <- gsub("_" , "", per_read_metrics$pid)

  registerDoMC(cores = config$ncpu)
  uniq_pids <- unique(per_read_metrics$clean_pid)
  if (is.null(bins_to_process)){
    bins_to_process <- length(uniq_pids)
  } else {
    bins_to_process <- min(length(uniq_pids), bins_to_process)
  }
  pid <- uniq_pids[583]
  working_dir <- op_dir
  tmp_x <- foreach(pid = uniq_pids[1:bins_to_process], .combine = "c") %dopar% {
    cur_seq_indxs <- which(per_read_metrics$clean_pid == pid)
    cur_seq_names <- per_read_metrics$read_name[cur_seq_indxs]
    stopifnot(all(cur_seq_names == as.character(seq_dat@id)[cur_seq_indxs]))

    cur_seqs <- seq_dat[cur_seq_indxs]

    if (length(cur_seqs) == 1){
      stop('cannot align single sequences. Bin should NOT be size 1')
    } else {
      aligned_with_qual <- alignBinsMSA_internal(cur_seqs, op_dir, pid)
      aligned_with_qual
    }
  }
  all_bins_aligned_with_qual <- shortReadQ_forced_append(tmp_x)
  rm(tmp_x)

  per_read_metrics <- data.frame('read_exists' = rep(1, length(all_bins_aligned_with_qual)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'alignBinsMSA'
  if (op_args$cache){
    result$seq_dat <- all_bins_aligned_with_qual
  }
  result$input_dat <- all_bins_aligned_with_qual
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

#' aligns fwd and rev sequences to a profile
#'
#' a custom guide tree is constructed to ensure that the alignment of the fwd
#' reads does not interfere with the alignment of the rev reads and vice versa.
#' @export

alignBinsMSA_internal <- function(cur_seqs, working_dir, pid)
{
  bin_seqs_file_name <- file.path(working_dir, 'bins', paste(pid, '_bin', '.fasta', sep = ''))
  bin_seqs <- cur_seqs@sread
  names(bin_seqs) <- as.character(cur_seqs@id)
  writeXStringSet(bin_seqs,
                  bin_seqs_file_name,
                  width=20000)
  aligned_file_name <- file.path(working_dir, 'bins', paste(pid, '_aligned', '.fasta', sep = ''))
  system(paste('mafft --quiet ', bin_seqs_file_name, 
               ' > ', aligned_file_name, sep = ''))
  stopifnot(file.exists(aligned_file_name))
  aligned_seqs <- readDNAStringSet(aligned_file_name)
  
  bin_seqs_quals <- cur_seqs@quality@quality
  names(bin_seqs_quals) <- paste(as.character(cur_seqs@id), sep = '_')
  
  gap_only_cols_cpp_indexing <-
  which(consensusMatrix(aligned_seqs)['-',] == length(aligned_seqs)) - 1
  
  reads_and_qual <- transfer_gaps_cpp(as.character(aligned_seqs),
                                      as.character(bin_seqs_quals), 
                                      gap_only_cols_cpp_indexing)
  
  qual_mat <- as(FastqQuality(reads_and_qual$quals), 'matrix')
  avg_quals <- trunc(apply(qual_mat, 1, (function(x) {mean(x[x>0])})))
  tweaked_qual_mat <- gapQualityTweaker_non_ol_cpp(reads_and_qual$reads, qual_mat, which_pair = 'fwd', avg_quals)

  aligned_with_qual <-
  ShortReadQ(sread = DNAStringSet(tweaked_qual_mat$reads),
             quality = BStringSet(tweaked_qual_mat$quals),
             id = BStringSet(names(aligned_seqs)))
  return(aligned_with_qual)
}

saveToDisk.alignBinsMSA <- function(result, config, seq_dat)
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

computeMetrics.alignBinsMSA <- function(result, config, seq_dat)
{
  return(result)
}

print.alignBinsMSA <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: alignBinsMSA')
  cat('\n-------------------')
  cat('\nKept Sequences:\n')
  print(result$summary[,c('parameter', 'k_seqs', 'k_mean_length', 'k_mean_qual')])
  cat('\n-------------------')
  cat('\nTrimmed Sequences:\n')
  print(result$summary[,c('parameter', 't_seqs', 't_mean_length', 't_mean_qual')])
  invisible(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
