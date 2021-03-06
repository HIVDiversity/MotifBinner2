# How does a summary table work?
# It summarizes the data that was kept was the parameters are relaxed.

excl_inf_max <- function(x)
{
  x <- x[x!= Inf]
  x <- x[x!= -Inf]
  x <- x[x!= 0]
  if (length(x) == 0){
    1
  } else {
    max(x)
  }
}

excl_inf_min <- function(x)
{
  x <- x[x!= Inf]
  x <- x[x!= -Inf]
  x <- x[x!= 0]
  if (length(x) == 0){
    1
  } else {
    min(x)
  }
}


#' Generates summary table when there is one trim step with one value
#' @inheritParams saveToDisk
#' @export

genSummary_case1 <- function(result, config, seq_dat)
{
  trim_step <- result$trim_steps[[1]]
  trim_dat <- result$metrics$per_read_metrics[,trim_step$name,drop=T]
  stopifnot(all(trim_dat == trim_step$threshold))
  kept <- seq_dat
  trimmed <- seq_dat[0]
  parameter <- paste(trim_step$name, ' = ', trim_step$threshold, sep = '')
  summary_tab <- genSummary_comb(kept = kept,
                                 trimmed = trimmed,
                                 op = class(result),
                                 parameter = parameter)
  return(summary_tab)
}

#' Generates summary table when there is one trim step with many values
#' @inheritParams saveToDisk
#' @export

genSummary_case3 <- function(result, config, seq_dat)
{
  summary_tab <- NULL
  trim_step <- result$trim_steps[[1]]
  trim_dat <- result$metrics$per_read_metrics[,trim_step$name,drop=T]
  if ('comparator' %in% names(trim_step))
  {
    comparator = trim_step$comparator
  } else {
    comparator = `<=`
  }
  for (break_indx in 2:length(trim_step$breaks))
  {
    kept_vec <- comparator(trim_dat, trim_step$breaks[break_indx]) &
                !comparator(trim_dat, trim_step$breaks[break_indx-1])
    trim_vec <- !comparator(trim_dat, trim_step$breaks[break_indx])
    kept_seqs <- seq_dat[kept_vec]
    curr_seq_dat <- seq_dat[trim_vec]
    if (excl_inf_max(c(abs(trim_step$breaks[break_indx-1]), abs(trim_step$breaks[break_indx]))) > 10000 |
        excl_inf_min(c(abs(trim_step$breaks[break_indx-1]), abs(trim_step$breaks[break_indx]))) < 0.01){
      parameter <- paste(trim_step$name, ' (', 
                         format(trim_step$breaks[break_indx-1], scientific = T, digits = 4), 
                         ',', 
                         format(trim_step$breaks[break_indx], scientific = T, digits = 4), 
                         ']', sep = '')
    } else {
      parameter <- paste(trim_step$name, ' (', 
                         round(trim_step$breaks[break_indx-1], 3), 
                         ',', 
                         round(trim_step$breaks[break_indx], 3), 
                         ']', sep = '')
    }
    if (comparator(trim_step$threshold, trim_step$breaks[break_indx]) &
        !(comparator(trim_step$threshold, trim_step$breaks[break_indx-1])))
    {
      parameter <- paste(parameter, ' *', sep = '')
    }
  
    summary_tab <- rbind(summary_tab,
      genSummary_comb(kept = kept_seqs,
                      trimmed = curr_seq_dat,
                      op = class(result),
                      parameter = parameter))
  }
  return(summary_tab)
}

#' Generates summary table when there is many trim steps each with many values
#' @inheritParams saveToDisk
#' @export

genSummary_case4 <- function(result, config, seq_dat)
{
  summary_tab <- NULL
  trim_criteria <- matrix(TRUE, nrow = length(seq_dat), ncol = length(result$trim_steps))
  i <- 0
  summary_tab <- NULL
  for (trim_step in result$trim_steps)
  {
    i <- i+1
    trim_dat <- result$metrics$per_read_metrics[,trim_step$name,drop=T]
    stopifnot(length(trim_dat) == length(seq_dat))
    if ('comparator' %in% names(trim_step))
    {
      comparator = trim_step$comparator
    } else {
      comparator = `<=`
    }
    for (break_indx in 2:length(trim_step$breaks))
    {
      kept_vec <- comparator(trim_dat, trim_step$breaks[break_indx]) &
                  !comparator(trim_dat, trim_step$breaks[break_indx-1])
      trim_vec <- !comparator(trim_dat, trim_step$breaks[break_indx])
      comp_kept_vec <- apply(trim_criteria, 1, all) & kept_vec
      comp_trim_vec <- apply(trim_criteria, 1, all) & trim_vec
      kept_seqs <- seq_dat[comp_kept_vec]
      curr_seq_dat <- seq_dat[comp_trim_vec]
      if (excl_inf_max(c(abs(trim_step$breaks[break_indx-1]), abs(trim_step$breaks[break_indx]))) > 10000 |
          excl_inf_min(c(abs(trim_step$breaks[break_indx-1]), abs(trim_step$breaks[break_indx]))) < 0.01){
        parameter <- paste(trim_step$name, ' (', 
                           format(trim_step$breaks[break_indx-1], scientific = T, digits = 4), 
                           ',', 
                           format(trim_step$breaks[break_indx], scientific = T, digits = 4), 
                           ']', sep = '')
      } else {
        parameter <- paste(trim_step$name, ' (', 
                           round(trim_step$breaks[break_indx-1], 3), 
                           ',', 
                           round(trim_step$breaks[break_indx], 3), 
                           ']', sep = '')
      }
      if (comparator(trim_step$threshold, trim_step$breaks[break_indx]) &
          !(comparator(trim_step$threshold, trim_step$breaks[break_indx-1])))
      {
        parameter <- paste(parameter, ' *', sep = '')
      }
      summary_tab <- rbind(summary_tab,
        genSummary_comb(kept = kept_seqs,
                        trimmed = curr_seq_dat,
                        op = class(result),
                        parameter = parameter))
    }
    trim_criteria[,i] <- comparator(trim_dat, trim_step$threshold)
  }
  return(summary_tab)
}

#' Generates a summary for a operation or the full process
#' @inheritParams saveToDisk
#' @export

genSummary <- function(result, config, seq_dat)
{
  op_args <- result$config$op_args
## 4 cases: (roughly order of complexity)
## #1: 1 trim_step + single value: everything kept -> hacky summary_tab
## #2: >1 trim_step + any step with single value: invalid
## #3: 1 trim_step + multiple values: normal process
## #4: >1 trim_steps + multiple values: complex:
#### iterate over trim steps and breaks.
#### complex building up of selection criteria to obtain correct datasets

## Specialize summaries for weird step - should maybe rather use OOP instead of
## ifs here...
  if (class(result) == 'matchPairs')
  {
    summary_tab <- genSummary_matchPairs(result)
  } else if (class(result) == 'processBadPIDs') {
    summary_tab <- genSummary_processBadPIDs(result)
  } else if (class(result) == 'dataTracing') {
    summary_tab <- NULL
  } else {
## So if a specialized genSummary is not needed:
  ## CASE 1: 1 trim_step + single_value
    if (length(result$trim_steps) == 1 & length(result$trim_steps[[1]]$breaks) == 1)
    {
      summary_tab <- genSummary_case1(result, config, seq_dat)
    }
  
  ## CASE 2: >1 trim_steps + any step with single value
    if (length(result$trim_steps) > 1)
    {
      for (i in 1:length(result$trim_steps)){
        if (length(result$trim_steps[[i]]$breaks) == 1){
          stop('dont use a single value trimming step in a trimming procedure with multiple trim steps - it does not make sense and it is untested')
        }
      }
    }
  
  ## CASE 3: 1 trim_step + multiple values
    if (length(result$trim_steps) == 1 & length(result$trim_steps[[1]]$breaks) > 1)
    {
      summary_tab <- genSummary_case3(result, config, seq_dat)
    }
  
  ## CASE 4: >1 trim_step + multiple values
    if (length(result$trim_steps) > 1){
      summary_tab <- genSummary_case4(result, config, seq_dat)
    }
  } ## else - this block runs the generic version of genSummary

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
  cat('\n\nTable: The number of sequences kept and discarded.This table shows the effect of increasingly less strict trimming based one or more parameters. The first column shows the name of the parameter and the range of values to which the row corresponds. Each row is applying the same filter with an expanded parameter range to those reads that were still excluded by the criteria specified in the previous row. Hence you will note that the seqs_in (number of sequences evaluated) column value always corresponds with the t_seqs (number of trimmed sequences) of the previous row. This table is useful to see if a small tweak to the trimming parameters will include or exclude a large number of reads.\n\n')
  #print( kable(summary_tab[,c('parameter', 'seqs_in', 'k_seqs', 't_seqs')]))
  print( format_table(summary_tab[,c('parameter', 'seqs_in', 'k_seqs', 't_seqs')], 
                      align = 'l'))
  cat('\n\nTable: Statistics of the kept sequences. The rows correspond to those of the previous table.\n\n')
  kept_tab <- summary_tab[,grep('k_', names(summary_tab))]
  names(kept_tab) <- gsub('k_', '', names(kept_tab))
#  print(kable(kept_tab, format='markdown'))
  print(format_table(kept_tab))
  cat('\n\nTable: Statistics of the trimmed sequences. The rows correspond to those of the previous two tables.\n\n')
  trimmed_tab <- summary_tab[,grep('t_', names(summary_tab))]
  names(trimmed_tab) <- gsub('t_', '', names(trimmed_tab))
#  print(kable(trimmed_tab, format='markdown'))
  print(format_table(trimmed_tab))
}
