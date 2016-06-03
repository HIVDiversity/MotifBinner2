#' This function creates a scoring matrix for aligning IUPAC sequences

makeIUPACScoreMat <- function()
{
  lets <- list(
  U = "T",
  A = "A",
  C = "C",
  M = c("A", "C"),
  G = "G",
  R = c("A", "G"),
  S = c("G", "C"),
  V = c("A", "C", "G"),
  T = "T",
  W = c("A", "T"),
  Y = c("C", "T"),
  H = c("A", "C", "T"),
  K = c("G", "T"),
  D = c("A", "G", "T"),
  B = c("C", "G", "T"),
  N = c("A", "C", "G", "T")
  )
  cat("    ")
  for (i in names(lets)){
    cat("   ")
    cat(i)
  }
  cat("\n")
  for (i in names(lets)){
    cat("   ")
    cat(i)
    for (j in names(lets)){
      ilets <- c(lets[[i]], i)
      jlets <- c(lets[[j]], j)
      any_match <- any(ilets %in% jlets)
      cat(ifelse(any_match, "  0,", " -1,"))
    }
    cat("\n")
  }

}
