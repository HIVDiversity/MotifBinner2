#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

binSizeCheck <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  seq_dat <- all_results[[data_source_indx]]$seq_dat
  
  per_read_metrics <- data.frame(seq_name = as.character(seq_dat@id),
                                 stringsAsFactors = F)
  per_read_metrics$pid <- gsub('^.*PID:([ACGT]*)_([ACGT]*).*$', '\\1\\2', per_read_metrics$seq_name)

  n_per_pid <- data.frame(table(per_read_metrics$pid))
  n_per_pid$Var1 <- as.character(n_per_pid$Var1)
  names(n_per_pid) <- c('pid', 'seq_left')

  per_read_metrics <- merge(per_read_metrics, n_per_pid,
                            all.x=TRUE, by.x = 'pid',
                            by.y = 'pid')
  seq_dat <- seq_dat[order(as.character(seq_dat@id))]
  per_read_metrics <- per_read_metrics[order(per_read_metrics$seq_name),]

  min_bin_size <- op_args$min_bin_size
  if (is.null(min_bin_size))
  {
    min_bin_size <- 3
  }
  trim_steps <- list(step1 = list(name = 'seq_left',
                                  threshold = min_bin_size,
                                  comparator = `>=`,
                                  breaks = c(Inf, min_bin_size+c(1,0,-1), -Inf)
                                  )
                    )

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  kept_dat <- getKept(result, seq_dat=seq_dat)
  class(result) <- 'binSizeCheck'
  if (op_args$cache){
    result$seq_dat <- kept_dat
  }
  result$input_dat <- seq_dat
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

saveToDisk.binSizeCheck <- function(result, config, seq_dat)
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

computeMetrics.binSizeCheck <- function(result, config, seq_dat)
{
  return(result)
}

print.binSizeCheck <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: binSizeCheck')
  cat('\n-------------------')
  cat('\nKept Sequences:\n')
  print(result$summary[,c('parameter', 'k_seqs', 'k_mean_length', 'k_mean_qual')])
  cat('\n-------------------')
  cat('\nTrimmed Sequences:\n')
  print(result$summary[,c('parameter', 't_seqs', 't_mean_length', 't_mean_qual')])
  invisible(result)
}
