MotifBinner
===========

MotifBinner processes high-throughput sequencing data of an RNA virus
population that was sequenced using the Primer ID approach as described in Jabara et
al, 2011. A random sequence tag is included in the initial primer during the
conversion from RNA to cDNA so that each input template is tagged with a unique
primer ID. After amplification with PCR and sequencing, each sequence with the
same PID should be from the same input template.

This data is cleaned by scanning each sequence for a primer ID. A prefix and a
suffix that surrounds the primer ID must be supplied.  If the prefix and the
suffix is found and the distance between them matches the primer ID length
(supplied by user), then the letters between the prefix and suffix is taken as
the primer ID for that sequence. Fuzzy matching is supported. All sequences
with the same primer ID are then grouped together in bins so that each bin
consists of sequences which had the same primer ID.

Bins that are made up of sequences whose primer IDs probably have sequencing
errors in them are then discarded using the consensus cutoff approach (see the
vignette on consensus cutoffs and Zhou et al, 2015 for more details). Each bin
is then inspected and the most outlying sequences are removed so that the bins
can be aligned without difficulty. The alignments are used to generate
consensus sequences using a majority rule. The consensus sequences together
with a report on the binning of the dataset is reported.

## Installation Instructions for Ubuntu 14.04

Make sure you have a recent version of R:
http://stackoverflow.com/questions/10476713/how-to-upgrade-r-in-ubuntu. Follow
these instructions to set up the correct repositiory for apt.

Make sure that both r-base and r-base-dev is installed
```{sh}
sudo apt-get install r-base r-base-dev
```

Next, install devtools' depedancies and muscle using apt-get:
```{sh}
sudo apt-get install libssl-dev libxml2-dev libcurl4-gnutls-dev muscle
```

A very recent (at the time this README was written) version of pandoc is required:
Download a binary package from
https://github.com/jgm/pandoc/releases/1.15 and then install it with dpkg:
```{sh}
wget https://github.com/jgm/pandoc/releases/download/1.15/pandoc-1.15-1-amd64.deb
sudo dpkg -i pandoc-1.15-1-amd64.deb
```

Next install all the dependencies from Bioconductor. This is done from within
an R session:
```{r}
sudo R
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")
biocLite("ShortRead")
```

Then, from within R, install devtools:
```{r}
install.packages('devtools', repo = 'http://cran.rstudio.com/')
```

At last install MotifBinner:
From a local file:

```{r}
library(devtools)
install_local('/path/to/file/MotifBinner_x.y.z.tar.gz')
```

Please note that you must use install_local from devtools - install.packages
will not work. Change /path/to/file to the path to the installation file on
your computer and x.y.z to match the installation file you have.

Or using the bit_bucket repo:
```{r}
library(devtools)
install_bitbucket('hivdiversity/MotifBinner2', 
  auth_user = 'username', password = 'password')
```

Lastly, MotifBinner includes a script that can be run from the commandline. You
need to put this script somewhere convenient ('/usr/bin' for example)
```{r}
file.symlink(from = file.path(find.package('MotifBinner2'), 'MotifBinner2.R'),
             to = '/usr/bin')
```

## Usage

### Within R

```{r}
library(MotifBinner2)
help('process_file')
```

This will display the help for the main function in MotifBinner.

### From the command line

```{sh}
MotifBinner2 -h
```

or (depending on your installation):

```{sh}
MotifBinner2.R -h
```

This will display help for all the options and an example call to MotifBinner.

## Bibliography

* Cassandra B Jabara, Corbin D Jones, Jeffrey Roach, Jeffrey A Anderson, and
Ronald Swanstrom.  Accurate sampling and deep sequencing of the HIV- 1 protease
gene using a Primer ID. Proceedings of the National Academy of Sciences of the
United States of America, 108(50):20166–71, December 2011.

* Shuntai Zhou, Corbin Jones, Piotr Mieczkowski and Ronald Swanstrom. Primer ID
validates template sampling depth and greatly reduced the error rate of
next-generation sequencing of HIV1 genomic RNA populations. J Virol. 2015 Aug
15;89(16):8540-55. doi: 10.1128/JVI.00522-15. Epub 2015 Jun 3.

