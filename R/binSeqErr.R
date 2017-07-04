#' Computes the sequencing parameters for each bin
#' @inheritParams applyOperation
#' @export

binSeqErr <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  stopifnot(all(names(op_args$data_source) %in% c("bin_msa_fwd", "bin_msa_rev", "bin_msa_merged",
                                                  "cons_fwd", "cons_rev", "cons_merged",
                                                  "primer_err"))) # must have fwd and/or rev primers
  stopifnot(!op_args$cache) #makes no sense to cache this...

  ref_err_indx <- grep(op_args$data_source[['primer_err']], names(all_results))
  ref_err_dat <- all_results[[ref_err_indx]]$metrics$error_parameters$all

  # Case 1 - merged
  if ("bin_msa_merged" %in% names(op_args$data_source))
  {
    stopifnot("cons_merged" %in% names(op_args$data_source))
    msa_indx <- grep(op_args$data_source[['bin_msa_merged']], names(all_results))
    cons_indx <- grep(op_args$data_source[['cons_merged']], names(all_results))

    msa_dat <- all_results[[msa_indx]]$seq_dat
    cons_dat <- all_results[[cons_indx]]$seq_dat

    all_tallies <- binSeqErr_internal(msa_dat = msa_dat, cons_dat = cons_dat)
  }

  # Case 2 - non-merged
  if ("bin_msa_fwd" %in% names(op_args$data_source))
  {
    stopifnot(all(c("bin_msa_fwd", "bin_msa_rev", "cons_fwd", "cons_rev", "primer_err") %in% names(op_args$data_source)))
    
    msa_indx <- grep(op_args$data_source[['bin_msa_fwd']], names(all_results))
    cons_indx <- grep(op_args$data_source[['cons_fwd']], names(all_results))

    msa_dat_fwd <- all_results[[msa_indx]]$seq_dat
    cons_dat_fwd <- all_results[[cons_indx]]$seq_dat

    all_tallies_fwd <- binSeqErr_internal(msa_dat = msa_dat_fwd, cons_dat = cons_dat_fwd)
    
    msa_indx <- grep(op_args$data_source[['bin_msa_rev']], names(all_results))
    cons_indx <- grep(op_args$data_source[['cons_rev']], names(all_results))

    msa_dat_rev <- all_results[[msa_indx]]$seq_dat
    cons_dat_rev <- all_results[[cons_indx]]$seq_dat

    all_tallies_rev <- binSeqErr_internal(msa_dat = msa_dat_rev, cons_dat = cons_dat_rev)

    all_tallies <- merge(all_tallies_fwd, all_tallies_rev, 
                         all.x = TRUE, all.y = TRUE,
                         by = c("from", "to", "qual", "data_source"))
    all_tallies$count.x[is.na(all_tallies$count.x)] <- 0
    all_tallies$count.y[is.na(all_tallies$count.y)] <- 0
    all_tallies$count <- all_tallies$count.x + all_tallies$count.y
    all_tallies$count.x <- NULL
    all_tallies$count.y <- NULL
    msa_dat <- shortReadQ_forced_append(list('fwd' = msa_dat_fwd,
                                             'rev' = msa_dat_rev))
  }

  per_read_metrics <- data.frame('read_exists' = rep(1, length(msa_dat)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics,
                                bin_sequencing_stats = all_tallies,
                                ref_err_dat = ref_err_dat))
  class(result) <- 'binSeqErr'
  if (op_args$cache){
    stop('Do not cache data in steps that do not alter the datasets')
  }
  result$input_dat <- msa_dat
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

#' Group consensus sequences and their read together
#'
#' Constructs a list with a single element for each bin. This element is a list
#' that contains the consensus sequence as well as the aligned reads of the
#' bins.
#' @inheritParams binSeqErr_internal

group_cons_bins_internal <- function(cons_dat, msa_dat)
{
  cons_pids <- gsub('_[^ACGTacgt]*$', '', as.character(cons_dat@id)) # strip off bin size counts
  cons_pids <- gsub('^.*_', '', cons_pids) # strip off base for names prefix
  cons_dat <- cons_dat[order(cons_pids)]
  cons_pids <- sort(cons_pids)
  msa_pids <- gsub('_[^ACGTacgt]*$', '', gsub('^.*PID:_', '', as.character(msa_dat@id)))
  msa_dat <- msa_dat[order(msa_pids)]
  msa_pids <- sort(msa_pids)
  stopifnot(all(msa_pids %in% cons_pids))

  # construct a list of all the consensus sequences and their associated bins
  dat_list <- list()
  cons_pids_i <- 1
  current_msa_indxs <- NULL
  for (msa_pids_i in 1:length(msa_pids)){
    if (msa_pids[msa_pids_i] != cons_pids[cons_pids_i]){
      # moved to next pid
      # save the current pid's sequences
      dat_list[[cons_pids[cons_pids_i]]] <- list(cons_dat_bin = cons_dat[cons_pids_i],
                                                msa_dat_bin = msa_dat[current_msa_indxs])
      # move/reset indexes to next pids
      cons_pids_i <- cons_pids_i + 1
      stopifnot(msa_pids[msa_pids_i] == cons_pids[cons_pids_i])
      current_msa_indxs <- NULL
    } else {
      current_msa_indxs <- c(current_msa_indxs, msa_pids_i)
    }
  }
  # remember to save the last element in the dat_list!!!
  dat_list[[cons_pids[cons_pids_i]]] <- list(cons_dat_bin = cons_dat[cons_pids_i],
                                            msa_dat_bin = msa_dat[current_msa_indxs])
  return(dat_list)
}

#' Reformats the tallies produced by tallyPrimerSeqError_cpp
#'
#' @param tallies_list The tallies list produced by tallyPrimerSeqErrors_cpp
#' @return A data.frame with the tallies of the sequencing parameters

reformat_tallies <- function(tallies_list)
{
  tallies_list_of_dfs <- list()
  inner_loop_counter <- 1
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
        tallies_list_of_dfs[[inner_loop_counter]] <- data.frame(from = ci,
          to = cj, qual = ck, count = tallies_list[[i]][[j]][k], 
          stringsAsFactors = FALSE)
        inner_loop_counter <- inner_loop_counter + 1
      }
    }
  }
  return(data.frame(data.table::rbindlist(tallies_list_of_dfs)))
}

#' Computes the sequencing error rate for each bin
#'
#' Under the assumption that there was zero PCR recombination and no PID
#' collisions, the only source of heterogeneity in the bins are sequencing
#' errors. Thus a comparison between each read in a bin and the bin's consensus
#' sequence can be taken as an indication of the sequencing error rates.
#'
#' A full characterization of the errors is returned that tallies each match
#' and mismatch and the quality at which they happened.
#'
#' @param cons_dat The consensus sequences. Their names must start with the PID
#'   followed by an underscore and then anything.
#' @param msa_dat The reads with sequences of each bin aligned. The sequence
#'   name must contain PID:_ followed by the PID somewhere in the name.

binSeqErr_internal <- function(cons_dat, msa_dat)
{
  dat_list <- group_cons_bins_internal(cons_dat = cons_dat, msa_dat = msa_dat)

  all_tallies <- NULL
  all_tallies <- foreach(bin_name = names(dat_list)) %dopar% {
    bin_seq <- dat_list[[bin_name]]$msa_dat_bin
    cons_seq <- dat_list[[bin_name]]$cons_dat_bin
    cons_seq <- rep(as.character(cons_seq@sread), length(bin_seq))

    tallies_list <- tallyPrimerSeqErrors_cpp(as.character(bin_seq@sread), 
                                             cons_seq, 
                                             as.character(bin_seq@quality@quality))

    tallies_df <- reformat_tallies(tallies_list)
    tallies_df$data_source <- bin_name

    tallies_df
  }
  return(data.frame(data.table::rbindlist(all_tallies)))
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
    sub_rate <- 1-(sum(x$match_countx) / (sum(x$countx) - n_dels - n_ins))

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
