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

  all_tallies <- NULL
  seq_dat <- NULL
  for (data_source_name in names(op_args$data_source)){
    data_source_indx <- grep(op_args$data_source[data_source_name], names(all_results))
    stopifnot(length(data_source_indx) == 1)

    dat <- all_results[[data_source_indx]]$metrics$per_read_metrics
    seq_dat <- shortReadQ_forced_append(list(seq_dat, all_results[[data_source_indx]]$seq_dat))
    n_fragments <- length(grep('read_fragment_', names(dat)))
    reads <- NULL
    primers <- NULL
    quals <- NULL
    for (i in 1:n_fragments){
      reads <- paste(reads, dat[,paste('read_fragment_', i, sep = '')], sep = '')
      primers <- paste(primers, dat[,paste('prefix_fragment_', i, sep = '')], sep = '')
      quals <- paste(quals, dat[,paste('read_qual_fragment_', i, sep = '')], sep = '')
    }
#    stopifnot(all(sapply(reads, nchar) == sapply(primers, nchar)))
#    stopifnot(all(sapply(reads, nchar) == sapply(quals, nchar)))

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
  class(result) <- 'primerSeqErr'
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

saveToDisk.primerSeqErr <- function(result, config, seq_dat)
{
  return(result)
}

computeMetrics.primerSeqErr <- function(result, config, seq_dat)
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
