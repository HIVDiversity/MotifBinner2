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

  profile_seqs <- readDNAStringSet(op_args$profile_file)

  pid <- unique(per_read_metrics$clean_pid)[1]
#  counter <- 1
#  loop_length <- length(unique(per_read_metrics$clean_pid))
#  ptm <- proc.time()
  registerDoMC(cores = 8)
  uniq_pids <- unique(per_read_metrics$clean_pid)
  if (is.null(bins_to_process)){
    bins_to_process <- length(uniq_pids)
  } else {
    bins_to_process <- min(length(uniq_pids), bins_to_process)
  }
  tmp_x <- foreach(pid = uniq_pids[1:bins_to_process], .combine = "c") %dopar% {
#  for (pid in unique(per_read_metrics$clean_pid)){
#    elapsed_time <- proc.time() - ptm
#    print (c(counter, loop_length, round(counter/loop_length,3), elapsed_time[3], elapsed_time[3]/round(counter/loop_length,3)  ))
#    counter <- counter + 1
    cur_seqs <- which(per_read_metrics$clean_pid == pid)
    cur_seq_names <- per_read_metrics$read_name[cur_seqs]
    stopifnot(all(cur_seq_names == as.character(seq_dat_fwd@id)[cur_seqs]))
    stopifnot(all(cur_seq_names == as.character(seq_dat_rev@id)[cur_seqs]))
    cur_fwd_seqs <- seq_dat_fwd[cur_seqs]
    cur_rev_seqs <- seq_dat_rev[cur_seqs]

    interleaving_vector <- NULL
    for (fwd_indx in 1:length(cur_fwd_seqs)){
      rev_indx <- length(cur_fwd_seqs) + fwd_indx
      interleaving_vector <- c(interleaving_vector, fwd_indx, rev_indx)
    }

    interleaved_seqs <- c(cur_fwd_seqs@sread, reverseComplement(cur_rev_seqs@sread))[interleaving_vector]
    names(interleaved_seqs) <- c(paste(as.character(cur_fwd_seqs@id), 'fwd', sep = '_'),
                                 paste(as.character(cur_fwd_seqs@id), 'rev', sep = '_'))[interleaving_vector]
    interleaved_seqs <- c(profile_seqs, interleaved_seqs)
#    if (!all(cur_fwd_seqs@sread %in% interleaved_seqs)){
#        stop('not all fwd reads in interleaved seqs - something is fishy')
#    }
#    if (!all(reverseComplement(cur_rev_seqs@sread) %in% interleaved_seqs)){
#        stop('not all rev reads in interleaved seqs - time to handle this edge case now...')
#    }
    mafft_guide_tree <- data.frame(r1 = 1,
                                   r2 = 2:length(interleaved_seqs),
                                   r3 = 0.01,
                                   r4 = 0.01)
    gt_file_name <- file.path(op_dir, 'bins', paste(pid, '_guide_tree.txt', sep = ''))
    write.table(mafft_guide_tree,
                gt_file_name,
                col.names = FALSE, row.names = FALSE, sep = '\t')
    interleaved_file_name <- file.path(op_dir, 'bins', paste(pid, '_interleaved', '.fasta', sep = ''))
    writeXStringSet(interleaved_seqs,
                    interleaved_file_name,
                    width=20000)
    aligned_file_name <- file.path(op_dir, 'bins', paste(pid, '_aligned', '.fasta', sep = ''))
    system(paste('mafft --quiet --retree 1 --treein ', gt_file_name, ' ', interleaved_file_name, 
                 ' > ', aligned_file_name, sep = ''))
    stopifnot(file.exists(aligned_file_name))
    aligned_seqs <- readDNAStringSet(aligned_file_name)
    aligned_seqs <- aligned_seqs[!(names(aligned_seqs) %in% names(profile_seqs))]
    interleaved_quals <- c(cur_fwd_seqs@quality@quality, reverse(cur_rev_seqs@quality@quality))[interleaving_vector]
    names(interleaved_quals) <- c(paste(as.character(cur_fwd_seqs@id), 'fwd', sep = '_'),
                                  paste(as.character(cur_fwd_seqs@id), 'rev', sep = '_'))[interleaving_vector]

#    output for c++ developemnt and testing
#    x <- paste('!', interleaved_quals, sep = '')
#    x <- BStringSet(x)
#    names(x) <- names(interleaved_quals)
#    writeXStringSet(x, gsub('fasta', 'bstring', gsub('interleaved', 'interleaved_quals', interleaved_file_name)), compress=F, width=20000)
  
     gap_only_cols_cpp_indexing <-
     which(consensusMatrix(aligned_seqs)['-',] == length(aligned_seqs)) - 1
     
     reads_and_qual <- transfer_gaps_cpp(as.character(aligned_seqs),
                                         as.character(interleaved_quals), 
                                         gap_only_cols_cpp_indexing)
     
     qual_mat <- as(FastqQuality(reads_and_qual$quals), 'matrix')
     cat('\n\n\n\n\n'); Rcpp::sourceCpp('src/trimFront.cpp')
     tweaked_qual_mat <- gapQualityTweaker_cpp(reads_and_qual$reads, qual_mat)

     aligned_with_qual <-
     ShortReadQ(sread = DNAStringSet(tweaked_qual_mat$reads),
                quality = BStringSet(tweaked_qual_mat$quals),
                id = BStringSet(names(aligned_seqs)))

# for normal for loop
#     if (is.null(all_bins_aligned_with_qual )){
#       all_bins_aligned_with_qual <- aligned_with_qual
#     } else {
#       all_bins_aligned_with_qual <- shortReadQ_forced_append(all_bins_aligned_with_qual,
#                                                              aligned_with_qual)
#     }
    aligned_with_qual
  }
  all_bins_aligned_with_qual <- tmp_x[[1]]
  for (i in 2:length(tmp_x)){
    all_bins_aligned_with_qual <- shortReadQ_forced_append(all_bins_aligned_with_qual,
                                                           tmp_x[[i]])
  }
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
