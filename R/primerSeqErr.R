#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

primerSeqErr <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  stopifnot(length(op_args$data_source) == 1 | length(op_args$data_source) == 2)
  stopifnot(all(names(op_args$data_source) %in% c("fwd", "rev"))) # must have fwd and/or rev primers
  stopifnot(!op_args$cache) #makes no sense to cache this...

  data_source_name <- names(op_args$data_source)[2]
  for (data_source_name in names(op_args$data_source)){
    data_source_indx <- grep(op_args$data_source[data_source_name], names(all_results))
    stopifnot(length(data_source_indx) == 1)

    dat <- all_results[[data_source_indx]]$metrics$per_read_metrics
    n_fragments <- length(grep('read_fragment_', names(dat)))
    reads <- NULL
    primers <- NULL
    quals <- NULL
    for (i in 1:n_fragments){
      reads <- paste(reads, dat[,paste('read_fragment_', i, sep = '')], sep = '')
      primers <- paste(primers, dat[,paste('prefix_fragment_', i, sep = '')], sep = '')
      quals <- paste(quals, dat[,paste('read_qual_fragment_', i, sep = '')], sep = '')
    }
    stopifnot(all(sapply(reads, nchar) == sapply(primers, nchar)))
    stopifnot(all(sapply(reads, nchar) == sapply(quals, nchar)))

    tallies_list <- tallyPrimerSeqErrors_cpp(reads, primers, quals)
    tallies_df <- data.frame(from = character(0),
                    to = character(0),
                    qual = character(0),
                    count = numeric(0),
                    stringsAsFactors = FALSE)
    for (i in names(tallies_list)){
      ci <- intToUtf8(i)
      if (ci == '1') {ci <- 'ins'}
      if (ci == '2') {ci <- 'del'}
      for (j in names(tallies_list[[i]])){
        cj <- intToUtf8(j)
        if (cj == '1') {cj <- 'ins'}
        if (cj == '2') {cj <- 'del'}
        for (k in names(tallies_list[[i]][[j]])){
          ck <- intToUtf8(k)
          tallies_df <- rbind(tallies_df, 
                              data.frame(from = ci,
                                         to = cj,
                                         qual = ck,
                                         count = tallies_list[[i]][[j]][k],
                                         stringsAsFactors = FALSE))      
        }
      }
    }
    tallies_df$data_source <- data_source_name
    
#    write.csv(data.frame(reads = reads, primers = primers, quals = quals), '/fridge/data/primer_seq_err.csv',
#                         row.names=F)


  }






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
  class(result) <- 'primerSeqErr'
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

saveToDisk.primerSeqErr <- function(result, config, seq_dat)
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

genSummary_primerSeqErr <- function(result)
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

computeMetrics.primerSeqErr <- function(result, config, seq_dat)
{
  return(result)
}

print.primerSeqErr <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: primerSeqErr')
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
