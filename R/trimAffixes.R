#' Trims the ends of the sequences so that they all start at the same place
#' @inheritParams applyOperation
#' @export

trimAffixes <- function(all_results, config)
{
  op_dir <- file.path(config$output_dir, config$prefix_for_names,
                      paste('n', sprintf("%03d", length(all_results)+1), '_trimAffixes', sep = ''))
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)

  kept <- list()
  trimmed <- list()
  metrics <- list()

  fwd_primer <- config$fwd_primer
  rev_primer <- config$rev_primer

  fwd_seq <- all_results[[length(all_results)]]$kept[['fwd_reads']]
  fwd_result <- trimAffixes_internal(fwd_seq, fwd_primer)

  rev_seq <- all_results[[length(all_results)]]$kept[['rev_reads']]
  rev_result <- trimAffixes_internal(rev_seq, rev_primer)

  result <- list(fwd_result = fwd_result,
                 rev_result = rev_result)
#  result <- list(kept = kept,
#                 trimmed = trimmed,
#                 metrics = metrics,
#                 step_num = length(all_results)+1,
#                 op_dir = op_dir)
  class(result) <- 'trimAffixes'
  return(result)
}

trimAffixes_internal <- function(seq_dat, prefix, min_score = 0.7, front_gaps_allowed = 0)
{
  if (is.null(min_score)) {min_score <- -Inf}
  if (min_score < 1 & min_score > 0){min_score <- -nchar(prefix)*(1-min_score)}

  trimmed <- trimEnds_cpp(as.character(seq_dat@sread),
                          as.character(seq_dat@id),
                          as.character(seq_dat@quality@quality),
                          prefix)
  trimmed <- list(seq_dat = ShortReadQ(sread = DNAStringSet(trimmed$sread),
                                       quality = BStringSet(trimmed$qual),
                                       id = BStringSet(trimmed$id)),
                  trim_stats = data.frame(score = trimmed$score,
                                          bases_trimmed = trimmed$trim_spot,
                                          gaps_at_front_of_read = trimmed$gaps_at_front_of_read))

  kept_list <- trimmed$trim_stats$score > min_score &
               trimmed$trim_stats$gaps_at_front_of_read <= front_gaps_allowed
  kept <- list(seq_dat = trimmed$seq_dat[kept_list],
               trim_stats = trimmed$trim_stats[kept_list,])
  trimmed <- list(seq_dat = trimmed$seq_dat[!kept_list],
                  trim_stats = trimmed$trim_stats[!kept_list,])
  return(list(kept = kept,
              trimmed = trimmed))
}



saveToDisk.trimAffixes <- function(result, config)
{
  for (data_set_name in names(result$kept)){
    seq_dat <- result$kept[[data_set_name]]
    if (length(seq_dat) > 0)
    {
      writeFastq(seq_dat, file.path(result$op_dir, 
        paste(config$prefix_for_names, '_kept_', data_set_name, '.fastq', sep = '')), compress=F)
    }
  }
  for (data_set_name in names(result$trimmed)){
    seq_dat <- result$trimmed[[data_set_name]]
    if (length(seq_dat) > 0)
    {
      writeFastq(seq_dat, file.path(result$op_dir, 
        paste(config$prefix_for_names, '_trimmed_', data_set_name, '.fastq', sep = '')), compress=F)
    }
  }
  return(result)
}

genSummary.trimAffixes <- function(result, config)
{
  summary_tab <- rbind(
    genSummary_internal(operation = 'trimAffixes',
                        parameters = 'fwd_reads',
                        kept_seq_dat = result$kept$fwd_reads,
                        trimmed_seq_dat = result$trimmed$fwd_reads),
    genSummary_internal(operation = 'trimAffixes',
                        parameters = 'rev_reads',
                        kept_seq_dat = result$kept$rev_reads,
                        trimmed_seq_dat = result$trimmed$rev_reads))
  result$summary <- summary_tab
  write.csv(summary_tab, file.path(result$op_dir, 'trimAffixes_summary.csv'), row.names=FALSE)
  return(result)
}

computeMetrics.trimAffixes <- function(result, config)
{
  return(result)
}

print.trimAffixes <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: trimAffixes')
  cat('\n-------------------')
  print(names(result))
#  cat('\nKept Sequences:\n')
#  print(result$summary[,c('parameters', 'seqs_kept', 'mean_length_kept', 'mean_qual_kept')])
#  cat('\n-------------------')
#  cat('\nTrimmed Sequences:\n')
#  print(result$summary[,c('parameters', 'seqs_trimmed', 'mean_length_trimmed', 'mean_qual_trimmed')])
  return(result)
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


