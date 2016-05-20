#' Generates a summary for a operation or the full process
#' @inheritParams saveToDisk
#' @export

genSummary <- function(result, config)
{
  UseMethod('genSummary')
}

#' Computes the summary tables of an operation
#' @export

genSummary_internal <- function(operation, parameters, good_seq_dat, bad_seq_dat)
{
  cbind(data.frame(operation = operation,
                   parameters = parameters),
        genSummary_internal_one(good_seq_dat),
        genSummary_internal_one(bad_seq_dat))
}

genSummary_internal_one <- function(seq_dat)
{
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

  data.frame(seqs = length(seq_dat),
             q25_length = quantile(read_widths, 0.25, names=FALSE),
             mean_length = mean(read_widths, na.rm=TRUE),
             q75_length = quantile(read_widths, 0.75, names=FALSE),
             q25_qual = quantile(per_read_quality, 0.25, names = FALSE),
             mean_qual = mean(per_read_quality, na.rm=T),
             q75_qual = quantile(per_read_quality, 0.75, names = FALSE),
             prop_gaps = prop_gaps,
             prop_non_ACGT = prop_non_ACGT,
             prop_AT = prop_AT)
}
