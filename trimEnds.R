library(Rcpp)
library(ShortRead)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

#seq_dat <- readFastq('/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq')
seq_dat <- readFastq('test_dat.fastq')
seq_dat <- seq_dat[1:4]

sourceCpp('trimEnds.cpp')
prefix <- "TATGGGAYSAAAGYCTMAARCCATGTG"

trimEnds_internal <- function(seq_dat, prefix, min_score = 0.7, front_gaps_allowed = 0)
{
  if (is.null(min_score)) {min_score <- -Inf}
  if (min_score < 1 & min_score > 0){min_score <- -nchar(prefix)*(1-min_score)}

  trimmed <- trimEnds_cpp(as.character(seq_dat@sread),
                          as.character(seq_dat@id),
                          as.character(seq_dat@quality@quality),
                          prefix)
  trimmed <- list(seq_dat = ShortReadQ(sread = DNAStringSet(trimmed$sread),
                                       quality = BStringSet(trimmed$qual),
                                       id = BStringSet(trimmed$id)),
                  trim_stats = data.frame(score = trimmed$score,
                                          bases_trimmed = trimmed$trim_spot,
                                          first_nongap = trimmed$first_nongap))

  kept_list <- trimmed$trim_stats$score > min_score &
               trimmed$trim_stats$first_nongap <= front_gaps_allowed
  kept <- list(seq_dat = trimmed$seq_dat[kept_list],
               trim_stats = trimmed$trim_stats[kept_list,])
  trimmed <- list(seq_dat = trimmed$seq_dat[!kept_list],
                  trim_stats = trimmed$trim_stats[!kept_list,])
  return(list(kept = kept,
              trimmed = trimmed))
}

trimmed <- trimEnds_internal(seq_dat, prefix)


