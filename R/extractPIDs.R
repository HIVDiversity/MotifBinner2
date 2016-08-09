#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

extractPIDs <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  pid_in_which_fragment <- op_args$pid_in_which_fragment
  pid_gaps_allowed <- op_args$pid_gaps_allowed
  pattern_to_chop_from_names <- op_args$pattern_to_chop_from_names
  
  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)

  seq_dat <- all_results[[data_source_indx]]$seq_dat
  pid_holding_metrics <- all_results[[data_source_indx]]$metrics$per_read_metrics
  pid_holding_metrics <- pid_holding_metrics[pid_holding_metrics$id %in% as.character(seq_dat@id),]

  if (is.null(pid_in_which_fragment)){
    per_read_metrics <- pid_holding_metrics[,'id',drop=F]
    per_read_metrics$total_pid_gaps <- 0
    per_read_metrics$read <- ''
    trim_breaks <- c(-Inf, 0, Inf)
  } else {
    per_read_metrics <- pid_holding_metrics[,c('id', 
        paste('prefix_fragment_', pid_in_which_fragment, sep = ''), 
        paste('read_fragment_', pid_in_which_fragment, sep = ''), 
        paste('read_qual_fragment_', pid_in_which_fragment, sep = ''), 
        paste('read_gaps_fragment_', pid_in_which_fragment, sep = ''),
        paste('prefix_gaps_fragment_', pid_in_which_fragment, sep = '')
      )]
    names(per_read_metrics) <- gsub('_fragment_[0-9]$', '', names(per_read_metrics))
    per_read_metrics$total_pid_gaps <- per_read_metrics$read_gaps + per_read_metrics$prefix_gaps
    trim_breaks <- c(-Inf, 0:3, Inf)
  }
  seq_ids <- as.character(seq_dat@id)
  per_read_metrics <- per_read_metrics[match(seq_ids, per_read_metrics$id),]
  seq_ids <- gsub(pattern_to_chop_from_names, '', seq_ids)
  per_read_metrics$id <- seq_ids
  seq_ids <- paste(seq_ids, '_PID:', per_read_metrics$read, sep = '')
  seq_dat@id <- BStringSet(seq_ids)
  trim_steps <- list(step1 = list(name = 'total_pid_gaps',
                                  threshold = pid_gaps_allowed,
                                  comparator = `<=`,
                                  breaks = trim_breaks
                                  )
                    )

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  kept_dat <- getKept(result, seq_dat=seq_dat)
  class(result) <- 'extractPIDs'
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

saveToDisk.extractPIDs <- function(result, config, seq_dat)
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

computeMetrics.extractPIDs <- function(result, config, seq_dat)
{
  return(result)
}

print.extractPIDs <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: extractPIDs')
  cat('\n-------------------')
  cat('\nKept Sequences:\n')
  print(result$summary[,c('parameter', 'k_seqs', 'k_mean_length', 'k_mean_qual')])
  cat('\n-------------------')
  cat('\nTrimmed Sequences:\n')
  print(result$summary[,c('parameter', 't_seqs', 't_mean_length', 't_mean_qual')])
  invisible(result)
}

