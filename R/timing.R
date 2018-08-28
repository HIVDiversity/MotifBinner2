#' Internal function that reforms a list storing exe times
#' @param x a list of execution times

sliceTime <- function(x)
{
  z <- summary(x)['user']
  names(z) <- NULL
  z
}

sliceTime_allResults <- function(x)
{
  all_operations <- names(x)
  all_operations <- all_operations[grep('^n[0-9]{3}_', all_operations)]
  all_times <- NULL
  for (op in all_operations)
  {
    z <- sapply(x[[op]]$timing, sliceTime)
    z['Total'] <- sum(z)
    z <- round(z, 2)
    z['Operation'] <- op
    z <- t(data.frame(z))
    row.names(z) <- NULL
    all_times <- rbind(all_times, z)
  }
  all_times
}

#' Produces the timing tables in the reports
#' @param result A result object with a timing element
#' @export

timingTable <- function(result)
{
  if (class(result) == 'allResults')
  {
    timingTable_allResults(result)
  } else {
    timingTable_other(result)
  }
}

timingTable_allResults <- function(result)
{
  cat(paste('\nTable: Timing for all the operations', '.\n', sep = ''))
  #print(kable(sliceTime_allResults(result)))
  print(format_table(sliceTime_allResults(result), 
                     align = 'l'))
}

timingTable_other <- function(result)
{
  cat(paste('\nTable: Timing for the steps in ', class(result), '.\n', sep = ''))
#  print(
#    kable(data.frame(step = names(sapply(result$timing, sliceTime)),
#               time = sapply(result$timing, sliceTime),
#               row.names = NULL))
#  )
  print(
    format_table(data.frame(step = names(sapply(result$timing, sliceTime)),
                   time = sapply(result$timing, sliceTime),
                   row.names = NULL),
                 align = 'l')
  )
}
