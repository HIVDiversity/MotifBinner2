library(Rcpp)
library(ShortRead)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

seq_dat <- readFastq('test_dat.fastq')
seq_dat <- seq_dat[1:4]

sourceCpp('trimEnds.cpp')
prefix <- "TATGGGAYSAAAGYCTMAARCCATGTG"
prefix <- "ATATATATATATATATATATATATATA"

trimEnds_cpp(as.character(seq_dat@sread),
             as.character(seq_dat@id),
             as.character(seq_dat@quality@quality),
             prefix)
print(seq_dat)


