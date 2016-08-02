#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

binSeqErr <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  stopifnot(length(op_args$data_source) == 2)
  stopifnot(all(names(op_args$data_source) %in% c("bin_msa", "cons"))) # must have fwd and/or rev primers
  stopifnot(!op_args$cache) #makes no sense to cache this...

  msa_indx <- grep(op_args$data_source['bin_msa'], names(all_results))
  msa_dat <- all_results[[msa_indx]]$seq_dat

  cons_indx <- grep(op_args$data_source['cons'], names(all_results))
  cons_dat <- all_results[[cons_indx]]$seq_dat

  all_tallies <- NULL
  bin_name <- as.character(cons_dat@id)[1]

  counter <- 0
  for (bin_name in as.character(cons_dat@id)){
    counter <- counter+1
    print(counter)
    cons_seq <- cons_dat[as.character(cons_dat@id) == bin_name]
    bin_seq <- msa_dat[grep(bin_name, as.character(msa_dat@id))]
    
    cons_seq <- rep(as.character(cons_seq@sread), length(bin_seq))

    tallies_list <- tallyPrimerSeqErrors_cpp(as.character(bin_seq@sread), 
                                             cons_seq, 
                                             as.character(bin_seq@quality@quality))
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
    tallies_df$data_source <- bin_name
    all_tallies <- rbind(all_tallies, tallies_df)
  }

  if (length(unique(all_tallies$data_source))){
    x <- all_tallies %>%
      group_by(from, to, qual) %>%
      summarize(count = sum(count))
    x <- as.data.frame(x)
    x$data_source <- 'all'
    all_tallies <- rbind(all_tallies, as.data.frame(x))
  }

  per_read_metrics <- data.frame('read_exists' = rep(1, length(seq_dat)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics,
                                primer_sequencing_stats = all_tallies))
  class(result) <- 'binSeqErr'
  if (op_args$cache){
    stop('Do not cache data in steps that do not alter the datasets')
  }
  result$input_dat <- seq_dat
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

saveToDisk.binSeqErr <- function(result, config, seq_dat)
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

genSummary_binSeqErr <- function(result)
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

computeMetrics.binSeqErr <- function(result, config, seq_dat)
{

  #silly little private function
  let_to_int <- function(x){
    if (x == 'A') return(1)
    if (x == 'C') return(2)
    if (x == 'G') return(3)
    if (x == 'T') return(4)
  }

  error_parameters <- list()
  for (data_source_name in unique(result$metrics$primer_sequencing_stats$data_source))
  {
    dat <- subset(result$metrics$primer_sequencing_stats,
                  data_source == data_source_name)
    dat$match <- dat$from == dat$to
    dat$match_count <- dat$match * dat$count
    
    denominator_for_indels <- sum(dat$count) - sum(subset(dat, from == 'del', select = 'count', drop=T))
    del_rate <- sum(subset(dat, from == 'del', select = 'count', drop=T))/denominator_for_indels
    ins_rate <- sum(subset(dat, from == 'ins', select = 'count', drop=T))/denominator_for_indels
    #1 - (n_deletions / denominator_for_indels)
    #1 - (n_insertions / denominator_for_indels)
    
    subs_dat <- subset(dat, 
                       from %in% c('A', 'C', 'G', 'T') &
                         to %in% c('A', 'C', 'G', 'T'))
    subs_rate <- 1-(sum(subs_dat$match_count) / sum(subs_dat$count))
    
    qual_subs <- subs_dat %>%
           group_by(match, qual) %>%
           summarize(numerator = sum(count))
    qual_subs <- data.frame(qual_subs)
    qual_totals <- subs_dat %>%
                     group_by(qual) %>%
                     summarize(denominator = sum(count))
    qual_totals <- data.frame(qual_totals)
    qual_subs <- merge(qual_subs, qual_totals, all.x = TRUE)
    
    subs_by_qual <- 
    subset(qual_subs, 
           match == TRUE, 
           select = c('qual', 'denominator', 'numerator'))
    names(subs_by_qual) <- c('qual', 'total', 'matches')
    mismatches_by_qual <- 
    subset(qual_subs, 
           match == FALSE, 
           select = c('qual', 'numerator'))
    names(mismatches_by_qual) <- c('qual', 'mismatches')
    subs_by_qual <- merge(subs_by_qual, mismatches_by_qual, all.x = TRUE)
    subs_by_qual$matches[is.na(subs_by_qual$matches)] <- 0
    subs_by_qual$mismatches[is.na(subs_by_qual$mismatches)] <- 0
    subs_by_qual <- subset(subs_by_qual, total > 1000)
    
    subs_by_qual$num_qual <- sapply(subs_by_qual$qual, utf8ToInt) - 33
    subs_by_qual$rate <- subs_by_qual$mismatches / subs_by_qual$total
    subs_by_qual$SE <- sqrt(subs_by_qual$rate*(1-subs_by_qual$rate)/subs_by_qual$total)
    subs_by_qual$E <- qnorm(.975)*subs_by_qual$SE
    subs_by_qual$LB <- subs_by_qual$rate - subs_by_qual$E
    subs_by_qual$UB <- subs_by_qual$rate + subs_by_qual$E
    
    
    
    from_let <- 'A'
    to_let <- 'A'
    to_let <- 'T'
    rate_mat <- matrix(rep(0, 4*4),ncol=4,nrow=4)
    for (from_let in c('A', 'C', 'G', 'T')){
      for (to_let in c('A', 'C', 'G', 'T')){
        count <- sum(subset(dat, from == from_let & to == to_let, select='count', drop=T))
        total <- sum(subset(dat, from == from_let, select='count', drop=T))
        rate <- count / total
        rate_mat[let_to_int(from_let), let_to_int(to_let)] <- rate
      }
    }
    rate_mat <- data.frame(round(rate_mat, 5))
    names(rate_mat) <- c('A', 'C', 'G', 'T')
    row.names(rate_mat) <- c('A', 'C', 'G', 'T')


    error_parameters[[data_source_name]] <- list(
      subs_rate = subs_rate,
      ins_rate = ins_rate,
      del_rate = del_rate,
      rate_mat = rate_mat,
      subs_by_qual = subs_by_qual)
  }
  result$metrics$error_parameters <- error_parameters
  return(result)
}

print.binSeqErr <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: binSeqErr')
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
