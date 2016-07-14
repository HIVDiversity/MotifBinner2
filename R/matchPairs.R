#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

matchPairs <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  stopifnot(length(op_args$data_source) == 2)
  stopifnot(sort(names(op_args$data_source)) == c("fwd", "rev"))
  stopifnot(op_args$cache)

  data_source_indx_fwd <- grep(op_args$data_source['fwd'], names(all_results))
  stopifnot(length(data_source_indx_fwd) == 1)
  data_source_indx_rev <- grep(op_args$data_source['rev'], names(all_results))
  stopifnot(length(data_source_indx_rev) == 1)
  
  seq_dat_fwd <- all_results[[data_source_indx_fwd]]$seq_dat
  seq_dat_rev <- all_results[[data_source_indx_rev]]$seq_dat

  tmp <- data.frame(x=as.character(seq_dat_fwd@id),stringsAsFactors=F)
  fwd_names <- separate(data = tmp, col = x, into = c('raw_name', 'pid_fwd'), sep = ".PID:")
  tmp <- data.frame(x=as.character(seq_dat_rev@id),stringsAsFactors=F)
  rev_names <- separate(data = tmp, col = x, into = c('raw_name', 'pid_rev'), sep = ".PID:")

  merged_names <- merge(fwd_names, rev_names, by = 'raw_name')
  merged_names$new_name <- paste(merged_names$raw_name,
                                 "_PID:",
                                 merged_names$pid_fwd,
                                 '_',
                                 merged_names$pid_rev,
                                 sep = '')
  fwd_kept <- seq_dat_fwd[match(paste(merged_names$raw_name, '_PID:', merged_names$pid_fwd, sep = ''), 
                                as.character(seq_dat_fwd@id))]
  fwd_trim <- seq_dat_fwd[!(as.character(seq_dat_fwd@id) %in% paste(merged_names$raw_name, '_PID:', merged_names$pid_fwd, sep = ''))]

  rev_kept <- seq_dat_rev[match(paste(merged_names$raw_name, '_PID:', merged_names$pid_rev, sep = ''), 
                                as.character(seq_dat_rev@id))]
  rev_trim <- seq_dat_rev[!(as.character(seq_dat_rev@id) %in% paste(merged_names$raw_name, '_PID:', merged_names$pid_rev, sep = ''))]
  fwd_kept@id <- BStringSet(

  per_read_metrics <- data.frame('read_exists' = rep(1, length(fwd_kept)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  kept_dat <- list(fwd = fwd_kept,
                   rev = rev_kept)
  class(result) <- 'matchPairs'
  result$seq_dat <- kept_dat
  result$input_dat <- list(fwd = seq_dat_fwd,
                           rev = seq_dat_rev)
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

saveToDisk.matchPairs <- function(result, config, seq_dat)
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

computeMetrics.matchPairs <- function(result, config, seq_dat)
{
  return(result)
}

print.matchPairs <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: matchPairs')
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
