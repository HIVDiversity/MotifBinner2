#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

qualTrim <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  seq_dat <- all_results[[data_source_indx]]$seq_dat

  qual_mat <- as(FastqQuality(quality(quality(seq_dat))), 'matrix')

  avg_qual_threshold <- op_args$avg_qual
  bad_base_threshold <- op_args$bad_base_threshold
  max_bad_bases <- op_args$max_bad_bases
  if (is.null(avg_qual_threshold)) { avg_qual_threshold <- 20 }
  if (is.null(bad_base_threshold)) { bad_base_threshold <- 10 }
  if (is.null(max_bad_bases)) { max_bad_bases <- .05 }
  
  count_under <- function(x, threshold=10){sum(x<threshold, na.rm=T)}
  
  per_read_metrics <- data.frame(seq_len = width(seq_dat@sread),
                                 avg_qual = apply(qual_mat, 1, mean, na.rm=T),
                                 bad_bases = apply(qual_mat, 1, count_under, op_args$bad_base_threshold)
                                 )
  per_read_metrics$perc_bad <- per_read_metrics$bad_bases / per_read_metrics$seq_len

  trim_steps <- list(step1 = list(name = 'avg_qual',
                                  threshold = avg_qual_threshold,
                                  comparator = `>=`,
                                  breaks = c(Inf, avg_qual_threshold+c(5,1,0,-5), -Inf)
                                  ),
                     step2 = list(name = 'perc_bad',
                                  threshold = max_bad_bases,
                                  comparator = `<=`,
                                  breaks = c(-Inf, max_bad_bases+c(-0.02, -0.01, 0, 0.01), Inf)
                                  )
                    )

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  kept_dat <- getKept(result, seq_dat=seq_dat)
  class(result) <- 'qualTrim'
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

saveToDisk.qualTrim <- function(result, config, seq_dat)
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

computeMetrics.qualTrim <- function(result, config, seq_dat)
{
  return(result)
}

print.qualTrim <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: qualTrim')
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
