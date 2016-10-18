#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

dataTracing <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  i <- 1
  chain_stats <- NULL
  for (i in 1:length(op_args$data_source)){
    chain <- strsplit(names(op_args$data_source)[i], '\\.')[[1]][1]
    link_n <- as.numeric(strsplit(names(op_args$data_source)[i], '\\.')[[1]][2])
    link_indx <- grep(op_args$data_source[i], names(all_results))
    link <- all_results[[link_indx]]
    link_name <- link$config$op_full_name
    link_output <- sum(genKeptVector(link))

    chain_stats <- rbind(chain_stats,
      data.frame(chain = chain,
                 link_n = link_n,
                 link_name = link_name,
                 link_output = link_output,
                 stringsAsFactors = FALSE))
    rm(link)
  }

  result <- list(metrics = list(chain_stats = chain_stats))
  class(result) <- 'dataTracing'
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

saveToDisk.dataTracing <- function(result, config, seq_dat)
{
  return(result)
}

computeMetrics.dataTracing <- function(result, config, seq_dat)
{
  return(result)
}

print.dataTracing <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: dataTracing')
  cat('\n-------------------\n')
  print(result$metrics$chain_stats)
  invisible(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
