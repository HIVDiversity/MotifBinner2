library(Rcpp)
library(ShortRead)

seq_dat <- readFastq('test_dat.fastq')

prefix <- 'TATGGGAYSAAAGYCTMAARCCATGTG'

print(seq_dat@sread)

cat('\n\n\n\n\n\n\n\n\n'); sourceCpp('trimFront.cpp')
trimmed <- trimFront_cpp(as.character(seq_dat@sread),
                        as.character(seq_dat@quality@quality),
                        prefix, 13, 7, 7,
                        10)

trimmed$sread

