library(Rcpp)
library(ShortRead)

seq_dat <- readFastq('test_dat.fastq')

sourceCpp('trimEnds.cpp')

print(seq_dat)
print_cpp(as.character(seq_dat@sread),
          as.character(seq_dat@id),
          as.character(seq_dat@quality@quality))
