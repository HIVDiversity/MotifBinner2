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

  fwd_raw_names <- data.frame(x=as.character(seq_dat_fwd@id),stringsAsFactors=F)
  fwd_names <- separate(data = fwd_raw_names, col = x, into = c('raw_name', 'pid_fwd'), sep = ".PID:")
  rev_raw_names <- data.frame(x=as.character(seq_dat_rev@id),stringsAsFactors=F)
  rev_names <- separate(data = rev_raw_names, col = x, into = c('raw_name', 'pid_rev'), sep = ".PID:")

  if (config$header_format == "SRA"){
    fwd_names$raw_name <- gsub(" length=[0-9]*", "", gsub("\\.1 ", "_", fwd_names$raw_name))
    rev_names$raw_name <- gsub(" length=[0-9]*", "", gsub("\\.2 ", "_", rev_names$raw_name))
    fwd_raw_names <- gsub(" length=[0-9]*_", "_", gsub("\\.1 ", "_", fwd_raw_names$x))
    rev_raw_names <- gsub(" length=[0-9]*_", "_", gsub("\\.2 ", "_", rev_raw_names$x))
  } else {
    fwd_raw_names <- fwd_raw_names$x
    rev_raw_names <- rev_raw_names$x
  }

  merged_names <- merge(fwd_names, rev_names, by = 'raw_name')
  merged_names$new_name <- paste(merged_names$raw_name,
                                 "_PID:",
                                 merged_names$pid_fwd,
                                 '_',
                                 merged_names$pid_rev,
                                 sep = '')
  
  
  fwd_kept <- seq_dat_fwd[match(paste(merged_names$raw_name, '_PID:', merged_names$pid_fwd, sep = ''), 
                                fwd_raw_names)]

  fwd_trim <- seq_dat_fwd[!(fwd_raw_names %in% paste(merged_names$raw_name, '_PID:', merged_names$pid_fwd, sep = ''))]

  rev_kept <- seq_dat_rev[match(paste(merged_names$raw_name, '_PID:', merged_names$pid_rev, sep = ''), 
                                rev_raw_names)]
  rev_trim <- seq_dat_rev[!(rev_raw_names %in% paste(merged_names$raw_name, '_PID:', merged_names$pid_rev, sep = ''))]
  fwd_kept@id <- BStringSet(merged_names$new_name)
  rev_kept@id <- BStringSet(merged_names$new_name)

  per_read_metrics <- rbind(
                      data.frame('has_pair' = rep(1, length(fwd_kept))),
                      data.frame('has_pair' = rep(0, length(fwd_trim) + length(rev_trim))))
  trim_steps <- list(step1 = list(name = 'has_pair',
                                  threshold = 1,
                                  comparator = `>=`,
                                  breaks = c(Inf, 1, 0, -Inf)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  kept_dat <- list(fwd = fwd_kept,
                   rev = rev_kept)
  trim_dat <- list(fwd = fwd_trim,
                   rev = rev_trim)
  class(result) <- 'matchPairs'
  result$seq_dat <- kept_dat
  result$trim_dat <- trim_dat
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
  kept <- result$seq_dat
  trimmed <- result$trim_dat

  if (length(kept[['fwd']]) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_fwd_kept_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(kept[['fwd']], tmp_name, compress=F)
  }
  if (length(kept[['rev']]) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_rev_kept_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(kept[['rev']], tmp_name, compress=F)
  }
  if (length(trimmed[['fwd']]) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_fwd_trimmed_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(trimmed[['fwd']], tmp_name, compress=F)
  }
  if (length(trimmed[['rev']]) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_rev_trimmed_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(trimmed[['rev']], tmp_name, compress=F)
  }
  return(result)
}

genSummary_matchPairs <- function(result)
{
  summary_tab <-
    rbind(
  genSummary_comb(kept = result$seq_dat$fwd,
                  trimmed = BiocGenerics::append(
                              BiocGenerics::append(
                                result$seq_dat$rev,
                                result$trim_dat$fwd),
                              result$trim_dat$rev),
                  op = class(result),
                  parameter = 'fwd_with_rev')
  ,
  genSummary_comb(kept = result$seq_dat$rev,
                  trimmed = BiocGenerics::append(result$trim_dat$fwd,
                                   result$trim_dat$rev),
                  op = class(result),
                  parameter = 'rev_with_fwd *')
  ,
  genSummary_comb(kept = result$trim_dat$fwd,
                  trimmed = result$trim_dat$rev,
                  op = class(result),
                  parameter = 'fwd_only')
  , 
  genSummary_comb(kept = result$trim_dat$rev,
                  trimmed = result$trim_dat$rev[0],
                  op = class(result),
                  parameter = 'rev_only')
  )
  return(summary_tab)
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
