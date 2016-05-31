library(Rcpp)
library(ShortRead)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

seq_dat <- readFastq('test_dat.fastq')

sourceCpp('trimEnds.cpp')

print(seq_dat)
print_cpp(as.character(seq_dat@sread),
          as.character(seq_dat@id),
          as.character(seq_dat@quality@quality))
