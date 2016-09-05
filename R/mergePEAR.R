#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

mergePEAR <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  stopifnot(length(op_args$data_source) == 1)
  stopifnot(op_args$cache)

  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  
  seq_dat_fwd <- all_results[[data_source_indx]]$seq_dat$fwd
  seq_dat_rev <- all_results[[data_source_indx]]$seq_dat$rev

  fwd_file_name <- file.path(op_dir, paste('fwd', '.fasta', sep = ''))
  writeFastq(seq_dat_fwd,
             fwd_file_name,
             compress=FALSE)
  rev_file_name <- file.path(op_dir, paste('rev', '.fasta', sep = ''))
  writeFastq(seq_dat_rev,
             rev_file_name,
             compress=FALSE)

  merged_file_name <- file.path(op_dir, paste('merged', '.fastq', sep = ''))
  system(paste('pear', 
               ' -f ', fwd_file_name, 
               ' -r ', rev_file_name,
               ' -j ', config$ncpu,
               ' -o ', merged_file_name, sep = ''))
  merged_dat <- readFastq(paste(merged_file_name, '.assembled.fastq', sep =''))
  discarded_dat <- shortReadQ_forced_append(list(
      'fwd' = readFastq(paste(merged_file_name, '.unassembled.forward.fastq', sep ='')),
      'rev' = readFastq(paste(merged_file_name, '.unassembled.reverse.fastq', sep =''))))
  file.remove(fwd_file_name)
  file.remove(rev_file_name)
  file.rename(paste(merged_file_name, '.discarded.fastq', sep = ''),
              gsub('.fastq$', '_trimmed_bad_merges.fastq', merged_file_name))
  file.rename(paste(merged_file_name, '.unassembled.forward.fastq', sep = ''),
              gsub('.fastq$', '_trimmed_fwd.fastq', merged_file_name))
  file.rename(paste(merged_file_name, '.unassembled.reverse.fastq', sep = ''),
              gsub('.fastq$', '_trimmed_rev.fastq', merged_file_name))

  per_read_metrics <- rbind(
                      data.frame('merged' = rep(1, length(merged_dat))),
                      data.frame('merged' = rep(0, length(discarded_dat))))
  trim_steps <- list(step1 = list(name = 'merged',
                                  threshold = 1,
                                  comparator = `>=`,
                                  breaks = c(Inf, 1, -Inf)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'mergePEAR'
  result$seq_dat <- merged_dat
  result$input_dat <- shortReadQ_forced_append(list(
                      'merged' = merged_dat,
                      'discarded' = discarded_dat))
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

saveToDisk.mergePEAR <- function(result, config, seq_dat = NULL)
{
  # output already produced by PEAR operation
  return(result)
}

computeMetrics.mergePEAR <- function(result, config, seq_dat)
{
  return(result)
}

print.mergePEAR <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: mergePEAR')
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
