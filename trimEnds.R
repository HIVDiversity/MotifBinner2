library(Rcpp)
library(ShortRead)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

#seq_dat <- readFastq('/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq')
seq_dat <- readFastq('test_dat.fastq')
seq_dat <- seq_dat[1:4]

sourceCpp('trimEnds.cpp')
prefix <- "TATGGGAYSAAAGYCTMAARCCATGTG"

trimmed <- 
trimEnds_cpp(as.character(seq_dat@sread),
             as.character(seq_dat@id),
             as.character(seq_dat@quality@quality),
             prefix)

trimmed <- list(seq_dat = ShortReadQ(sread = DNAStringSet(trimmed$sread),
                                     quality = BStringSet(trimmed$qual),
                                     id = BStringSet(trimmed$id)),
                trim_stats = data.frame(score = trimmed$score,
                                        bases_trimmed = trimmed$trim_spot,
                                        first_nongap = trimmed$first_nongap))
print(trimmed$seq_dat@sread)
print(trimmed$trim_stats)


