#' Prepares config list for inclusion in final report
#'
#' Converts the configuration used to call MotifBinner to a data.frame that can
#' be included in the final report.
#'
#' The configuration used to specify which opertation should be applied and
#' what settings to use is contained in a nested list of lists. These lists are
#' typically built using the buildConfig function. By including these settings
#' and the version of MotifBinner in the final report an easy way to document
#' the exact procedures that were used to produce a dataset.
#'
#' The main action performed by this prepConfig operation is just to convert
#' the nested list of lists that holds the config to a data.frame and to
#' provide a convenient way, using the standard approach of MotifBinner, to add
#' this information into the final report.
#'
#' @inheritParams applyOperation
#' @export

prepConfig <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)

  flat_op_list <- data.frame(op_num = character(0),
                             setting_1 = character(0),
                             setting_2 = character(0),
                             value = character(0),
                             stringsAsFactors = FALSE)

  input_files <- data.frame(op_num = character(0),
                            value = character(0),
                            stringsAsFactors = FALSE)

  other_settings <- data.frame(opt_name = character(0),
                               value = character(0),
                               stringsAsFactors = FALSE)

  in_loadData <- FALSE
  for (op_num in names(config$operation_list)){
    for (i in names(config$operation_list[[op_num]])){
      if (length(config$operation_list[[op_num]][[i]]) > 1){
        print('multi')
        for (j in names(config$operation_list[[op_num]][[i]])){
          the_value <- config$operation_list[[op_num]][[i]][[j]]
          flat_op_list <- rbind(flat_op_list,
            data.frame(op_num = op_num,
                       setting_1 = i,
                       setting_2 = j,
                       value = ifelse(is.null(the_value), "NULL", the_value),
                       stringsAsFactors = FALSE))
        }
      } else {
        the_value <- config$operation_list[[op_num]][[i]]
        if (!is.null(the_value)){
          if (the_value == "loadData"){
            in_loadData <- TRUE
          }
        }
        if (in_loadData & i == "data_source"){
          flat_op_list <- rbind(flat_op_list,
            data.frame(op_num = op_num,
                       setting_1 = i,
                       setting_2 = "",
                       value = "See note below",
                       stringsAsFactors = FALSE))
          input_files <- rbind(input_files,
            data.frame(op_num = op_num,
                       value = ifelse(is.null(the_value), "NULL", the_value),
                       stringsAsFactors = FALSE))

        } else {
          flat_op_list <- rbind(flat_op_list,
            data.frame(op_num = op_num,
                       setting_1 = i,
                       setting_2 = "",
                       value = ifelse(is.null(the_value), "NULL", the_value),
                       stringsAsFactors = FALSE))
        }
      }
      
    }
    in_loadData <- FALSE
  }

  for (opt_name in names(config)){
    if (opt_name != "operation_list")
    {
      other_settings <- rbind(other_settings,
        data.frame(opt_name = opt_name,
                   value = config[[opt_name]],
                   stringsAsFactors = FALSE))
    }
  }

  names(flat_op_list) <- c("Op. Number  ", "Setting  ", "Sub-Setting  ", "Value  ")

  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(flat_op_list = flat_op_list,
                                input_files = input_files,
                                other_settings = other_settings
                                ))
  class(result) <- 'prepConfig'
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

computeMetrics.prepConfig <- function(result, config, seq_dat)
{
  return(result)
}

saveToDisk.prepConfig <- function(result, config, seq_dat)
{
  return(result)
}

print.prepConfig <- function(result, config)
{
  cat('\n-------------------\n')
  cat(  'Operation: prepConfig\n')
  cat(  '-------------------\n')
#  cat('\nLoaded Sequences: Bases:\n')
#  print(result$seq_dat@sread)
#  cat('\nLoaded Sequences: Qualities:\n')
#  print(result$seq_dat@quality@quality)
#  cat('\nSummary:\n')
  print(result$summary[,c('parameter', 'k_seqs', 'k_mean_length', 'k_mean_qual')])
  invisible(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
