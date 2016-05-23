
#' Performs any final processing operations
#' @inheritParams applyOperation
#' @export

allResults <- function(all_results, config)
{
  return(all_results)
}

saveToDisk.allResults <- function(result, config)
{
  return(result)
}

genSummary.allResults <- function(result, config)
{
  summary_tab <- NULL
  for (operation in names(result))
  {
    summary_tab <- rbind(summary_tab, result[[operation]]$summary)
  }
  result$summary <- summary_tab
  write.csv(summary_tab, file.path(config$output_dir, config$prefix_for_names, 
                                   paste(config$prefix_for_names, '_summary.csv', sep = '')), 
            row.names=FALSE)
  return(result)
}

print.allResults <- function(result, config)
{
  cat('\nallResults\n')
  print(result$summary)
  return(result)
}

