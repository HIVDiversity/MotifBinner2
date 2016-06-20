#' Generates a summary for a operation or the full process
#' @inheritParams saveToDisk
#' @export

genSummary <- function(result, config)
{
  op_args <- result$config$op_args

  # load data
  if (!op_args$cache)
  {
    stop('implement getKept now')
  } else 
  {
    seq_dat <- result$seq_dat
  }

  for (trim_step in result$trim_steps)
  {
    trim_dat <- result$metrics$per_read_metrics[,trim_step$name,drop=T]
    stopifnot(length(trim_dat) == length(seq_dat))
    if (length(trim_step$breaks) == 1)
    {
      stopifnot(all(trim_dat == trim_step$threshold))
      kept <- seq_dat
      trimmed <- seq_dat[0]
      parameter <- paste(trim_step$name, ' = ', trim_step$threshold, sep = '')
      summary_tab <- genSummary_comb(kept = kept,
                                     trimmed = trimmed,
                                     op = class(result),
                                     parameter = parameter)
    } else {
      stop('implement trimming with multiple categories now')
    }
  }
  write.csv(
    summary_tab, 
    file.path(
      result$config$op_dir, 
      paste(result$config$op_full_name, '.csv', sep = '')),
    row.names = FALSE
    )
  result$summary <- summary_tab
  return(result)
}

genSummary_comb <- function(kept, trimmed, op, parameter)
{
  cbind(
    data.frame(op = op,
               parameter = parameter),
    genSummary_one(kept, 'k'),
    genSummary_one(trimmed, 't')
    )
}

genSummary_one <- function(seq_dat, names_pref)
{
  if (length(seq_dat) == 0){
    ret_val <-
    data.frame(seqs = 0,
               q25_length = NA,
               mean_length = NA,
               q75_length = NA,
               q25_qual = NA,
               mean_qual = NA,
               q75_qual = NA,
               prop_gaps = NA,
               prop_non_ACGT = NA,
               prop_AT = NA)
  } else {
    read_widths <- width(seq_dat)
    seq_qual <- seq_dat@quality
    per_read_quality <- letterFrequency(seq_qual@quality , names(encoding(seq_qual)))
    per_read_quality <- (per_read_quality %*% matrix(encoding(seq_qual), ncol=1))[,1]
    per_read_quality <- per_read_quality / read_widths
    conMat <- consensusMatrix(seq_dat@sread)
    tot_lets <- sum(conMat)
    prop_gaps <- round(sum(conMat['-',])/tot_lets,5)
    prop_non_ACGT <- round(sum(conMat[5:18,])/tot_lets,5)
    prop_AT <- round(sum(conMat[c('A','T'),])/tot_lets,5)

    ret_val <-
    data.frame(seqs = length(seq_dat),
               q25_length = round(quantile(read_widths, 0.25, names=FALSE),2),
               mean_length = round(mean(read_widths, na.rm=TRUE),2),
               q75_length = round(quantile(read_widths, 0.75, names=FALSE),2),
               q25_qual = round(quantile(per_read_quality, 0.25, names = FALSE),2),
               mean_qual = round(mean(per_read_quality, na.rm=T),2),
               q75_qual = round(quantile(per_read_quality, 0.75, names = FALSE),2),
               prop_gaps = round(prop_gaps,4),
               prop_non_ACGT = round(prop_non_ACGT,4),
               prop_AT = round(prop_AT,4))
  }
  names(ret_val) <- paste(names_pref, names(ret_val), sep = '_')
  return(ret_val)
}


#' formats a summary table for markdown
#' @param summary_tab A data.frame as produced by genSummary

kable_summary <- function(summary_tab)
{
  summary_tab$seqs_in <- summary_tab$k_seqs + summary_tab$t_seqs
  cat('\n\nThe number of sequences kept and discarded:\n\n')
  print(
  kable(summary_tab[,c('parameter', 'seqs_in', 'k_seqs', 't_seqs')])
  )
  cat('\n\nStatistics of the kept sequences\n\n')
  kept_tab <- summary_tab[,grep('k_', names(summary_tab))]
  names(kept_tab) <- gsub('k_', '', names(kept_tab))
  print(kable(kept_tab))
  cat('\n\nStatistics of the trimmed sequences\n\n')
  trimmed_tab <- summary_tab[,grep('t_', names(summary_tab))]
  names(trimmed_tab) <- gsub('t_', '', names(trimmed_tab))
  print(kable(trimmed_tab))
}
