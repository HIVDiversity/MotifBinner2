#' Trims the ends of the sequences so that they all start at the same place
#' @inheritParams applyOperation
#' @export

trimAffixes <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  seq_dat <- all_results[[data_source_indx]]$seq_dat

  primer_seq <- op_args$primer_seq
  primer_lens <- op_args$primer_lens

  tmp_result <- trimAffixes_internal(seq_dat, primer_seq, primer_lens)
  trim_seq_dat <- tmp_result$seq_dat
  per_read_metrics <- tmp_result$per_read_metrics

  min_score <- op_args$min_score
  if (min_score < 0) stop('Min score must be > 0')

  trim_steps <- list(step1 = list(name = 'read_front_gaps',
                                  threshold = op_args$front_gaps_allowed,
                                  breaks = c(-Inf, 0, 1, 2, 5, Inf)
                                  ),
                     step2 = list(name = 'score',
                                  threshold = min_score,
                                  comparator = `>=`,
                                  breaks = rev(c(-Inf, 0, round(0.5*nchar(primer_seq),0),
                                             round(0.75*nchar(primer_seq),0),
                                             nchar(primer_seq)-3,
                                             nchar(primer_seq)-2,
                                             nchar(primer_seq)-1,
                                             nchar(primer_seq)-0,
                                             Inf))
                                  )
                    )
  result <- list(trim_steps = trim_steps,                                                           
                 metrics = list(per_read_metrics = per_read_metrics))                               

  kept_dat <- getKept(result, seq_dat=trim_seq_dat)                                                      

  if (op_args$cache){
    result$seq_dat <- kept_dat
  }
  result$input_dat <- trim_seq_dat
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)

  class(result) <- 'trimAffixes'
  return(result)
}

trimAffixes_internal <- function(seq_dat, primer_seq, primer_lens)
{
  trimmed <- trimFront_cpp(as.character(seq_dat@sread),
                          as.character(seq_dat@quality@quality),
                          primer_seq, primer_lens)
  per_read_metrics <- data.frame(
    score = trimmed$score,
    read_front_gaps = trimmed$read_front_gaps,
    prefix_front_gaps = trimmed$prefix_front_gaps
    )
  for (i in 1:length(primer_lens)){
    tmp <- data.frame(
      read_fragment_ = trimmed$read_fragments[,i],
      read_qual_fragment_ = trimmed$read_qual_fragments[,i],
      prefix_fragment_ = trimmed$prefix_fragments[,i],
      read_gaps_fragment_ = trimmed$read_fragment_gaps[,i],
      read_bases_fragment_ = trimmed$read_fragment_bases[,i],
      prefix_gaps_fragment_ = trimmed$prefix_fragment_gaps[,i],
      prefix_bases_fragment_ = trimmed$prefix_fragment_bases[,i]
      )
    names(tmp) <- paste(names(tmp), i, sep = '')
    per_read_metrics <- cbind(per_read_metrics, tmp)
  }
  trimmed_seq_dat <- ShortReadQ(sread = DNAStringSet(trimmed$rest_of_read),
                                quality = BStringSet(trimmed$rest_of_read_qual),
                                id = seq_dat@id)
  return(list(seq_dat = trimmed_seq_dat,
              per_read_metrics = per_read_metrics))
}

saveToDisk.trimAffixes <- function(result, config, seq_dat)
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

computeMetrics.trimAffixes <- function(result, config, seq_dat)
{
  return(result)
}

print.trimAffixes <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: trimAffixes')
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
#library(Rcpp)
#library(ShortRead)
#
#Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
#
##seq_dat <- readFastq('/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq')
#seq_dat <- readFastq('test_dat.fastq')
#seq_dat <- seq_dat[1:4]
#
#sourceCpp('trimAffixes.cpp')
#prefix <- "TATGGGAYSAAAGYCTMAARCCATGTG"

#trimmed <- trimAffixes_internal(seq_dat, prefix)


