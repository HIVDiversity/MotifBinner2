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

  data_table_rbindlist <- function(...){
    data.table::rbindlist(..., use.names = FALSE, fill = FALSE, idcol=NULL)
  }

#  counter <- 0
  all_tallies <- foreach(bin_name = as.character(cons_dat@id)) %dopar% {
#  for (bin_name in as.character(cons_dat@id)){
#    counter <- counter+1
#    print(counter)
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
    #all_tallies <- rbind(all_tallies, tallies_df)
    tallies_df
  }
  all_tallies <- data.frame(data.table::rbindlist(all_tallies))

  per_read_metrics <- data.frame('read_exists' = rep(1, length(seq_dat)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics,
                                bin_sequencing_stats = all_tallies))
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
  return(result)
}

computeMetrics.binSeqErr <- function(result, config, seq_dat)
{
  no_weirdness <- subset(result$metrics$bin_sequencing_stats, 
                         from %in% c('A', 'C', 'G', 'T', 'del', 'ins') &
                         to %in% c('A', 'C', 'G', 'T', 'del', 'ins')
                         )
  no_weirdness$match <- no_weirdness$from == no_weirdness$to
  no_weirdness$match[no_weirdness$from %in% c('del', 'ins')] <- FALSE
  no_weirdness$match_count <- no_weirdness$match * no_weirdness$count

  data_source_name <- unique(result$metrics$bin_sequencing_stats$data_source)[1]

  error_parameters <- list()
  i <- 0
  for (data_source_name in unique(result$metrics$bin_sequencing_stats$data_source))
  {
    i <- i+1
    x <- no_weirdness %>%
           filter(data_source == data_source_name) %>%
           group_by(from, to) %>%
           summarize(countx = sum(count),
                     match_countx = sum(match_count))
    x <- data.frame(x)
    n_dels <- subset(x, from == 'del', select = 'countx', drop = T)
    if (length(n_dels) == 0) {n_dels <- 0}
    n_ins <- subset(x, from == 'ins', select = 'countx', drop = T)
    if (length(n_ins) == 0) {n_ins <- 0}

    del_rate <- n_dels / (sum(x$countx) - n_ins)
    ins_rate <- n_ins / (sum(x$countx) - n_ins)
    sub_rate <- sum(x$match_countx) / (sum(x$countx) - n_dels - n_ins)

    error_parameters[[i]] <- data.frame(bin = data_source_name,
                                        del_rate = del_rate,
                                        ins_rate = ins_rate,
                                        sub_rate = sub_rate,
                                        stringsAsFactors = FALSE)
  }
  error_parameters <- data.frame(data.table::rbindlist(error_parameters))

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
