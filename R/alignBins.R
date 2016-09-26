#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

alignBins <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(op_dir, 'bins'), showWarnings = FALSE, recursive = TRUE)
  
  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  stopifnot(all(sort(names(all_results[[data_source_indx]]$seq_dat)) == c("fwd", "rev")))
  seq_dat_fwd <- all_results[[data_source_indx]]$seq_dat$fwd
  seq_dat_rev <- all_results[[data_source_indx]]$seq_dat$rev
  stopifnot(all(seq_dat_fwd@id == seq_dat_rev@id))

  bins_to_process <- op_args$bins_to_process

  per_read_metrics <- data.frame(read_name = as.character(seq_dat_fwd@id),
                                 stringsAsFactors = F)
  per_read_metrics$pid <- gsub("^.*_PID:" , "", per_read_metrics$read_name)
  per_read_metrics$clean_pid <- gsub("_" , "", per_read_metrics$pid)

#  profile_seqs <- readDNAStringSet(op_args$profile_file)

  registerDoMC(cores = config$ncpu)
  uniq_pids <- unique(per_read_metrics$clean_pid)
  if (is.null(bins_to_process)){
    bins_to_process <- length(uniq_pids)
  } else {
    bins_to_process <- min(length(uniq_pids), bins_to_process)
  }
  tmp_x <- foreach(pid = uniq_pids[1:bins_to_process], .combine = "c") %dopar% {
    task_indicator_file_name <- file.path(op_dir, 'bins', paste('zz_running_', pid, '.txt', sep = ''))
    file.create(task_indicator_file_name)
    cur_seqs <- which(per_read_metrics$clean_pid == pid)
    cur_seq_names <- per_read_metrics$read_name[cur_seqs]
    stopifnot(all(cur_seq_names == as.character(seq_dat_fwd@id)[cur_seqs]))
    stopifnot(all(cur_seq_names == as.character(seq_dat_rev@id)[cur_seqs]))

    cur_fwd_seqs <- seq_dat_fwd[cur_seqs]
    cur_rev_seqs <- seq_dat_rev[cur_seqs]

#    aligned_with_qual <- alignBins_internal(cur_fwd_seqs, cur_rev_seqs, profile_seqs, op_dir, pid)
    aligned_with_qual <- alignBins_internal(cur_fwd_seqs, cur_rev_seqs, op_args$profile_file, op_dir, pid)
    file.remove(task_indicator_file_name)
    aligned_with_qual
  }
  all_bins_aligned_with_qual <- shortReadQ_forced_append(tmp_x)
  rm(tmp_x)

  per_read_metrics <- data.frame('read_exists' = rep(1, length(all_bins_aligned_with_qual)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'alignBins'
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

map_reads_with_mafft <- function(interleaved_seqs, working_dir, pid, profile_file, n_profile_seqs)
{
#  guide_tree_size <- length(interleaved_seqs)
#  if (guide_tree_size < n_profile_seqs){
#    guide_tree_size <- guide_tree_size + n_profile_seqs
#  }
#  mafft_guide_tree <- data.frame(r1 = 1,
#                                 r2 = 2:guide_tree_size,
#                                 r3 = 0.01,
#                                 r4 = 0.01)
#  gt_file_name <- file.path(working_dir, 'bins', paste(pid, '_guide_tree.txt', sep = ''))
#  write.table(mafft_guide_tree,
#              gt_file_name,
#              col.names = FALSE, row.names = FALSE, sep = '\t')
  interleaved_file_name <- file.path(working_dir, 'bins', paste(pid, '_interleaved', '.fasta', sep = ''))
  writeXStringSet(interleaved_seqs,
                  interleaved_file_name,
                  width=20000)
  aligned_file_name <- file.path(working_dir, 'bins', paste(pid, '_aligned', '.fasta', sep = ''))
  if (is.null(profile_file)){
    system(paste('mafft --quiet --op 2.5 --ep 0.4 --retree 1 --treein ', gt_file_name, ' ', interleaved_file_name, 
                 ' > ', aligned_file_name, sep = ''))
  } else if (file.exists(profile_file)){
#    system(paste('mafft --quiet --ep 0.2 --retree 1 --treein ', gt_file_name, 
#                 ' --addfragments ', interleaved_file_name, 
#                 ' ', profile_file,
#                 ' > ', aligned_file_name, sep = ''))
    system(paste('mafft --quiet --op 2.5 --ep 0.4 ', 
                 ' --addfragments ', interleaved_file_name, 
                 ' ', profile_file,
                 ' > ', aligned_file_name, sep = ''))
  } else {
    stop('profile specification issues')
  }
  stopifnot(file.exists(aligned_file_name))
  aligned_seqs <- readDNAStringSet(aligned_file_name)
  aligned_seqs
}

#' aligns fwd and rev sequences to a profile
#'
#' a custom guide tree is constructed to ensure that the alignment of the fwd
#' reads does not interfere with the alignment of the rev reads and vice versa.
#' @export

#alignBins_internal <- function(cur_fwd_seqs, cur_rev_seqs, profile_seqs, working_dir, pid)
alignBins_internal <- function(cur_fwd_seqs, cur_rev_seqs, profile_file, working_dir, pid)
{
  profile_seqs <- readDNAStringSet(profile_file)
  n_profile_seqs <- length(profile_seqs)
  interleaving_vector <- NULL
  for (fwd_indx in 1:length(cur_fwd_seqs)){
    rev_indx <- length(cur_fwd_seqs) + fwd_indx
    interleaving_vector <- c(interleaving_vector, fwd_indx, rev_indx)
  }

  interleaved_seqs <- c(cur_fwd_seqs@sread, reverseComplement(cur_rev_seqs@sread))[interleaving_vector]
  names(interleaved_seqs) <- c(paste(as.character(cur_fwd_seqs@id), 'fwd', sep = '_'),
                               paste(as.character(cur_fwd_seqs@id), 'rev', sep = '_'))[interleaving_vector]
#  interleaved_seqs <- c(profile_seqs, interleaved_seqs)

  aligned_seqs <- map_reads_with_mafft(interleaved_seqs, working_dir, pid, profile_file, n_profile_seqs)
  if (length(aligned_seqs) < length(interleaved_seqs)){
    warning (paste('PID ', pid, ' mapping FAILED'))
    file.create(file.path(working_dir, 'bins', paste('zzz_', pid, '_mapping_failed.txt', sep = '')))
    aligned_seqs <- map_reads_with_mafft(c(profile_seqs, interleaved_seqs), working_dir, 
                                         pid, profile_file = NULL, n_profile_seqs)
  }
  aligned_seqs <- aligned_seqs[(names(aligned_seqs) %in% names(interleaved_seqs))]

#  aligned_seqs <- aligned_seqs[!(names(aligned_seqs) %in% names(profile_seqs))]
  interleaved_quals <- c(cur_fwd_seqs@quality@quality, reverse(cur_rev_seqs@quality@quality))[interleaving_vector]
  names(interleaved_quals) <- c(paste(as.character(cur_fwd_seqs@id), 'fwd', sep = '_'),
                                paste(as.character(cur_fwd_seqs@id), 'rev', sep = '_'))[interleaving_vector]

  gap_only_cols_cpp_indexing <-
  which(consensusMatrix(aligned_seqs)['-',] == length(aligned_seqs)) - 1
  
  reads_and_qual <- transfer_gaps_cpp(as.character(aligned_seqs),
                                      as.character(interleaved_quals), 
                                      gap_only_cols_cpp_indexing)
  
  qual_mat <- as(FastqQuality(reads_and_qual$quals), 'matrix')
  tweaked_qual_mat <- gapQualityTweaker_ol_cpp(reads_and_qual$reads, qual_mat)

  aligned_with_qual <-
  ShortReadQ(sread = DNAStringSet(tweaked_qual_mat$reads),
             quality = BStringSet(tweaked_qual_mat$quals),
             id = BStringSet(names(aligned_seqs)))
  return(aligned_with_qual)
}

saveToDisk.alignBins <- function(result, config, seq_dat)
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

computeMetrics.alignBins <- function(result, config, seq_dat)
{
  return(result)
}

print.alignBins <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: alignBins')
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
