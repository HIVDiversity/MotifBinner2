#' Removes all instances of a specified character from a set of sequences
#'
#' Note that the current consensus building operation assigns a quality of zero
#' to both gaps and N's. These are the two primary types of characters to
#' remove and hence it does not make sense to report on the quality of the
#' letters that were removed.
#'
#' @inheritParams applyOperation
#' @export

removeChars <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  seq_dat <- all_results[[data_source_indx]]$seq_dat

  no_gaps <- removeChars_cpp(as.character(seq_dat@sread),
                         as.character(seq_dat@quality@quality),
                         '-')

  no_gaps_dat <- ShortReadQ(DNAStringSet(no_gaps$sread),
                            BStringSet(no_gaps$qual),
                            seq_dat@id)

  per_read_metrics <- data.frame(seq_length = width(no_gaps_dat@sread),
                                 n_gaps = no_gaps$n_chars)

  trim_steps <- list(step1 = list(name = 'n_gaps',
                                  threshold = Inf,
                                  comparator = `<=`,
                                  breaks = c(-Inf, 1, 10, 100, Inf)
                                  )
                    )

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  kept_dat <- no_gaps_dat
  class(result) <- 'removeChars'
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

saveToDisk.removeChars <- function(result, config, seq_dat)
{
  tmp_name <- file.path(result$config$op_dir, 
    paste(config$base_for_names, '_kept_', result$config$op_args$name, '.fastq', sep = ''))
  writeFastq(result$seq_dat, tmp_name, compress=F)
  return(result)
}

computeMetrics.removeChars <- function(result, config, seq_dat)
{
  return(result)
}

print.removeChars <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: removeChars')
  cat('\n-------------------')
  cat('\nKept Sequences:\n')
  print(result$summary[,c('parameter', 'k_seqs', 'k_mean_length', 'k_mean_qual')])
  cat('\n-------------------')
  cat('\nTrimmed Sequences:\n')
  print(result$summary[,c('parameter', 't_seqs', 't_mean_length', 't_mean_qual')])
  invisible(result)
}
