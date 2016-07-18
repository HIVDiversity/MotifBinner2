# Tools to simulate the consensus cutoff of Zhou et al 2015

#' Simulate a parent and its main offspring size
#'
#' @param n bin size (1000)
#' @param pl PID length (9)
#' @param phi Sequencing read error rate (1/50)
#' @export

sim_one_parent_main_off <- function(n = 1000, pl = 9, phi = (1/50)){
  pid <- sample(c('A', 'C', 'G', 'T'), pl, replace = TRUE)
  alpha <- rbinom(1, n*pl, phi)
  err_indx <- sample(1:(n*pl), alpha, replace = FALSE)
  read_pids <- rep(pid, n)

  wrong_lets <- sapply(read_pids[err_indx], mutate_letter)
  read_pids[err_indx] <- wrong_lets
  read_pids <- t(matrix(read_pids, nrow = pl))
  read_pids <- apply(read_pids, 1, paste0, collapse = '')
  count_pids <- sort(table(read_pids), TRUE)
  attr(count_pids, 'dimnames') <- NULL
  return(list(parent = count_pids[1], 
              main_off = count_pids[2]))
}

#' Mutates a nucleotide to another
#'
#' @param x The letter to mutate
#' @export

mutate_letter <- function(x){
  alphabet <- c('A', 'C', 'G', 'T')
  stopifnot(x %in% alphabet)
  alphabet <- alphabet[alphabet != x]
  return(sample(alphabet, 1))
}

#' Simulate range of parents and their main offspring sizes
#' 
#' @param n a vector of bin_sizes (seq(100, 10000, 100))
#' @param pl PID length (9)
#' @param phi Sequencing read error rate (1/50)
#' @param repeats The number of repeats to run for each bin size
#' @param vebose If TRUE print loop progress
#' @param ncpu Number of cores to use for simulation
#' @export

sim_all_parents_main_off <- function(ns = seq(100, 10000, 100), pl = 9, 
                                     phi = (1/50), repeats = 1000, 
                                     verbose = FALSE, ncpu = 8){
  results_big <- NULL
  registerDoMC(cores=ncpu) ##
  timing_dat <- NULL
  run_times <- NULL
  pre_bind_times <- NULL
  bin_sizes <- NULL
  plotted <- FALSE
  for (n in ns){
    start_time <- Sys.time()
    results_small <- NULL
    results_small <- foreach(i=1:repeats, .combine = rbind) %dopar% { ##
      x <- sim_one_parent_main_off(n, pl, phi)
      parent = x$parent
      main_off = x$main_off
      data.frame(n = n, i = i, parent = parent, main_off = main_off)
    }
    pre_bind_time <- Sys.time() - start_time
    results_big <- rbind(results_big, results_small)
    run_time <- Sys.time() - start_time
    timing_dat <- rbind(timing_dat, 
                        data.frame(bin_size = n,
                                   pre_bind_time = pre_bind_time,
                                   run_time = run_time))
    if (nrow(timing_dat) > 10){
      mod <- lm(run_time ~ bin_size, data = timing_dat)
      final_pred <- mod$coefficients[1] + mod$coefficients[2]*max(ns)
      names(final_pred) <- NULL
      comp_time <- ((final_pred - run_time)/2 + run_time)*sum(ns > n)
      if (verbose) {
        if(!plotted){
          plot(x = min(ns), y = run_time, 
               xlim=c(min(ns), max(ns)), ylim = c(0, final_pred*2),
               col = 'white')
          points(pre_bind_time ~ bin_size, data = timing_dat, col = 'red', cex = 2)
          points(run_time ~ bin_size, data = timing_dat, col = 'black', cex = 2)
          plotted <- TRUE
        } else {
          points(pre_bind_time ~ bin_size, data = timing_dat[nrow(timing_dat),], 
                 col = 'red', cex = 2)
          points(run_time ~ bin_size, data = timing_dat[nrow(timing_dat),], 
                 col = 'black', cex = 2)
          abline(mod, col = rgb(n/max(ns),0,0,0.1))
        }
        print(c(n, run_time, final_pred, comp_time)); print(Sys.time())
      }
    }
  }
  return(list(results = results_big,
              timing_dat = timing_dat))
}

#' Wrapper that performs simulation and process the results
#'
#' The function calls the sim_all_parents_main_off function and saves the
#' result to a folder together with a file that contains a description of the
#' simulation parameters.
#' 
#' @param output_dir The output directory where the simulation result must be
#' placed
#' @param file_name The name of the output file in which the results will be
#' stored. Must end in .csv.
#' @param size_start The smallest bin size to simulate
#' @param size_stop The largest bin size to simulate
#' @param size_step The size of the steps to take to reach the largest bin size
#' @param pl PID length (9)
#' @param phi Sequencing read error rate (1/50)
#' @param repeats The number of repeats to run for each bin size
#' @param vebose If TRUE print loop progress
#' @param ncpu Number of cores to use for simulation
#' @export

sim_wrapper <- function(output_dir, file_name, size_start = 100, 
                        size_stop = 10000, size_step = 100, pl = 9, 
                        phi = (1/50), repeats = 1000, verbose = FALSE, 
                        ncpu = 8){
  stopifnot(file.exists(output_dir))
  stopifnot(grepl('.csv$', file_name))
  ns <- seq(size_start, size_stop, size_step)
  x <- sim_all_parents_main_off(ns = ns, pl = pl, phi = phi, 
                                repeats = repeats, verbose = verbose, 
                                ncpu = ncpu)
  write.csv(x$results, file.path(output_dir, file_name), row.names=FALSE)
  index_file <- data.frame(dir = output_dir,
                           file_name = file_name,
                           size_start = size_start, 
                           size_stop = size_stop, 
                           size_step = size_step, 
                           pl = pl, 
                           phi = phi, 
                           repeats = repeats)
  write.csv(index_file, file.path(output_dir, gsub('.csv', '_index.csv', file_name)),
            row.names = FALSE)
  return(x$timing_dat)
}

#' Builds the consensus cutoff models given a dataset index file
#'
#' @param index_file The file that describes the available datasets
#' @param verbose Prints progress and fitting informations and plots models to
#' active graphics device.
#' @param models_file If not NULL, then the models data.frame will be written
#' to this file.
#' @export

build_cutoff_models <- function(index_file, verbose = FALSE,
                                models_file = NULL){
  datasets <- read.csv(index_file, stringsAsFactors = FALSE)
  models <- NULL
  for (i in 1:nrow(datasets)){
    if(verbose){print(datasets[i,])}
    dataset <- 
      datasets %>% 
      slice(i) %$%
      read.csv(file.path(dir, file_name), stringsAsFactors = F)
    dataset$main_off[is.na(dataset$main_off)] <- 0
    mv_stats <-
      dataset %$% 
      moving_statistic(main_off, parent, statistic = quantile, bw = 40, 
                       n_steps = 500, probs = 0.99) %>%
      subset(., n > quantile(.$n, 0.025))
    row.names(mv_stats) <- NULL
    nlc <- nls.control(maxiter = 1000,
                       minFactor = 1/(2^15))
    mod_lm <- lm(mv_stats$stat ~ mv_stats$x)
    mod_nls <- nls2(stat ~ a + b*x + d*(x^e), 
                data = mv_stats, 
                start = data.frame(a=c(-10,10), 
                                   b=c(0.001, 0.01), 
                                   d=c(0.2, 0.7), 
                                   e=c(0.2, 0.6)),
                control = nlc)
    stopifnot(sum(residuals(mod_nls)^2) < sum(residuals(mod_lm)^2))

    if (verbose){
      plot(mv_stats$stat ~ mv_stats$x)
      abline(mod_lm)
      lines(mv_stats$x, predict(mod_nls), lwd = 3, col = 'red')
    }
    models <- rbind(models,
      data.frame(input_file = datasets[i,'file_name'],
      pl = datasets[i,'pl'],
      phi = datasets[i,'phi'],
      n_obs = nrow(mv_stats),
      rss = sum(residuals(mod_nls)^2),
      a = summary(mod_nls)$coefficients[1,1],
      b = summary(mod_nls)$coefficients[2,1],
      d = summary(mod_nls)$coefficients[3,1],
      e = summary(mod_nls)$coefficients[4,1]))
  }
  print(models)
  if (!is.null(models_file)){
    if (file.exists(models_file)){
      old_models <- read.csv(models_file, stringsAsFactors = FALSE)
      models <- rbind(models, old_models)
    }
    write.csv(models, models_file)
  }
  return(models)
}

#' Plots the results of a consensus_cutoff simulation
#' @param dat The data simulated by sim_all_parents_main_off
#' @export

plot_consensus_cutoff_sim <- function(dat){
  mv_stats <- moving_statistic(dat$main_off, dat$parent, statistic = quantile,
                               bw = 50, n_steps = 100, probs = 0.99)
  p <- ggplot(dat, aes(x = parent, y = main_off)) + 
    geom_point(alpha = 0.02) +
    geom_smooth(se = TRUE)+ 
    geom_line(data = mv_stats, aes(x = x, y = stat))
  return(p)
}

#' Computes a moving statistic over a sliding window.
#' @export

moving_statistic <- function(y, x, statistic, bw = 5, n_steps = 30, ...){
  lb <- min(x) + bw/2
  ub <- max(x) - bw/2
  results <- NULL
  for (i in 1:n_steps){
    c_mid <- lb + (i-1)*((ub - lb)/(n_steps-1))
    c_lb <- c_mid - bw/2
    c_ub <- c_mid + bw/2
    c_dat <- y[x >= c_lb & x < c_ub]
    c_n <- length(c_dat)
    c_stat <- statistic(c_dat, ...)
    results <- rbind(results,
                     data.frame(x = c_mid,
                                i = i,
                                n = c_n,
                                stat = c_stat,
                                mean = mean(c_dat),
                                sd = sd(c_dat)))
  }
  return(results)
}

#' Makes a series of numbers monotonically increasing or decreasing
#'
#' Monotonicity is achieved by substituting a number with the number preceeding
#' it if it violates the monotonicity criteria.
#' @param x The vector of numbers to make monotonic
#' @param direction The direction in which the numbers should be made monotonic
#' @export

make_monotonic <- function(x, direction = 'increasing'){
  for (i in 2:length(x)){
    if (direction == 'increasing'){
      if (x[i-1] > x[i]){
        x[i] <- x[i-1]
      }
    } else if (direction == 'decreasing'){
      if (x[i-1] < x[i]){
        x[i] <- x[i-1]
      }
    } else {
      stop('invalid value for direction - specify either increasing or decreasing')
    }
  }
  return(x)
}

#' Returns the consensus cutoff to use based on a given models file
#'
#' This function will find the model that matches your input parameters and use
#' it to predict the appropriate consensus cutoff for your dataset
#' @param pl PID length (9)
#' @param main_offspring_size The size of the largest bin in the dataset
#' produced by the binner. NB! This must be the size of the bin BEFORE outlier
#' removal
#' @param models A data.frame containing the parameters for various models. As
#' produced by build_cutoff_models.
#' @param phi Sequencing read error rate (1/50)
#' @export

get_consensus_cutoff <- function(pl, main_offspring_size, models, phi=1/50){
  model <- models[models$pl == pl & round(models$phi,5) == round(phi,5),]
  stopifnot(nrow(model)==1)
  x <- main_offspring_size
  model$a + model$b * x + model$d * (x^model$e)
}

#' Applies the consensus cutoff to an existing consensuses file
#'
#' Reads in a consensuses file and a pid_seq_name file. The largest bin is
#' computed from the data in the pid_seq_name file. Using the information about
#' the largest bin and the supplied sequencing error rate and PID length the
#' consensus cutoff is computed. From the PIDs in the pid_seq_name file, a list
#' of all PIDs associated with bins larger than the consensus cutoff is made.
#' The sequences in the consensuses file is compared to that list and the ones
#' that are based on small bins are removed from the file. The consensuses are
#' then written out to a file with the same name as the input file, but the
#' suffix cc added to it.
#' @param pid_seq_name_file The name of the file that contains the list of
#' sequence_names and the pids found in them
#' @param consensuses_file The name of the file with the consensus sequences.
#' @param sequencing_error_rate The sequencing error rate to use when looking
#' up the consensus cutoff model (only required if 'consensus_cutoff' is not
#' specified)
#' @param motif_length The length of the pid.
#' @param consensus_cutoff The consensus_cutoff number to use. If this is
#' specified, that number will be used to choose which sequences to use.
#' Alternatively if it is set to 'model' (the default) then a consensus cutoff 
#' will be computed based on the size of the largest bin and the sequencing_error_rate.
#' @export

apply_consensus_cutoff <- function(pid_seq_name_file, consensuses_file, sequencing_error_rate = 1/50,
                                   motif_length = 8, consensus_cutoff = 'model'){
  stopifnot(grepl('.fasta$', consensuses_file))
  pid_counts <- read.csv(pid_seq_name_file, stringsAsFactors = F)
  pid_counts <- data.frame(pid = names(table(pid_counts$motif)),
                           count = as.numeric(table(pid_counts$motif)),
                           stringsAsFactors = F)
  row.names(pid_counts) <- NULL
  if (consensus_cutoff == 'model'){
    largest_bin <- max(pid_counts$count)
    cc <- get_consensus_cutoff(motif_length, largest_bin, 
                               get_cutoff_models(), sequencing_error_rate)
  } else {
    stopifnot(as.character(consensus_cutoff) != as.character(as.numeric(consensus_cutoff)))
    cc <- as.numeric(consensus_cutoff)
  }
  allowed_pids <- subset(pid_counts, count >= cc)$pid

  consensus_dat <- read_sequence_file(consensuses_file)
  x <- gregexpr(paste0('[ACGT]{', motif_length, '}'), names(consensus_dat))
  PIDs <- NULL
  for (i in seq_along(x)){
    if(length(x[[1]]) != 1) {stop('Non or multiple PIDS found in fasta header')}
    start_pos <- x[[i]]
    stop_pos <- start_pos + attr(x[[i]], 'match.length')
    PIDs <- c(PIDs, substr(names(consensus_dat)[i], start_pos, stop_pos))
  }
  consensus_dat <- consensus_dat[PIDs %in% allowed_pids]
  writeXStringSet(consensus_dat, gsub('.fasta$', '_cc.fasta', consensuses_file))
  return(NULL)
}

#' Returns a data.frame of the stored models
#' 
#' A number of consensus cutoff models have been pre-built. This function
#' returns a data.frame with there parameters and the parameters that describes
#' the datasets that they are applicable to
#'
#' @export

get_cutoff_models <- function(){
  built_models <- data.frame(
       input_file = c("cc_100_15k_2k_10_1in50_100.csv",
                      "cc_100_15k_2k_11_1in50_100.csv",
                      "cc_100_15k_2k_12_1in50_100.csv",
                      "cc_100_15k_2k_7_1in100_100.csv",
                      "cc_100_15k_2k_7_1in50_100.csv",
                      "cc_100_15k_2k_8_1in100_100.csv",
                      "cc_100_15k_2k_8_1in50_100.csv",
                      "cc_100_15k_2k_9_1in50_100.csv",
                      "cc_100_15k_2k_10_1in100_100.csv",
                      "cc_100_15k_2k_11_1in100_100.csv",
                      "cc_100_15k_2k_11_1in75_100.csv",
                      "cc_100_15k_2k_12_1in100_100.csv",
                      "cc_100_15k_2k_12_1in75_100.csv",
                      "cc_100_15k_2k_13_1in100_100.csv",
                      "cc_100_15k_2k_13_1in50_100.csv",
                      "cc_100_15k_2k_13_1in75_100.csv",
                      "cc_100_15k_2k_14_1in100_100.csv",
                      "cc_100_15k_2k_14_1in50_100.csv",
                      "cc_100_15k_2k_14_1in75_100.csv",
                      "cc_100_15k_2k_15_1in100_100.csv",
                      "cc_100_15k_2k_15_1in50_100.csv",
                      "cc_100_15k_2k_5_1in100_100.csv",
                      "cc_100_15k_2k_6_1in100_100.csv",
                      "cc_100_15k_2k_9_1in100_100.csv",
                      "cc_100_15k_2k_10_1in75_100.csv",
                      "cc_100_15k_2k_15_1in75_100.csv",
                      "cc_100_15k_2k_5_1in50_100.csv",
                      "cc_100_15k_2k_5_1in75_100.csv",
                      "cc_100_15k_2k_6_1in50_100.csv",
                      "cc_100_15k_2k_6_1in75_100.csv",
                      "cc_100_15k_2k_7_1in75_100.csv",
                      "cc_100_15k_2k_8_1in75_100.csv",
                      "cc_100_15k_2k_9_1in75_100.csv"), 
       pl = c(10L, 11L, 12L, 7L, 7L, 8L, 8L, 9L, 10L, 11L, 11L, 12L, 12L, 13L,
              13L, 13L, 14L, 14L, 14L, 15L, 15L, 5L, 6L, 9L, 10L, 15L, 5L, 5L,
              6L, 6L, 7L, 8L, 9L), phi = c(0.02, 0.02, 0.02, 0.01, 0.02, 0.01,
              0.02, 0.02, 0.01, 0.01, 0.0133333333333333, 0.01,
              0.0133333333333333, 0.01, 0.02, 0.0133333333333333, 0.01, 0.02,
              0.0133333333333333, 0.01, 0.02, 0.01, 0.01, 0.01,
              0.0133333333333333, 0.0133333333333333, 0.02, 0.0133333333333333,
              0.02, 0.0133333333333333, 0.0133333333333333, 0.0133333333333333,
              0.0133333333333333), n_obs = c(486L, 487L, 487L, 487L, 487L, 487L,
              487L, 487L, 487L, 487L, 487L, 486L, 487L, 487L, 487L, 487L, 487L,
              487L, 487L, 487L, 487L, 487L, 487L, 487L, 487L, 487L, 487L, 487L,
              487L, 487L, 487L, 487L, 487L), rss = c(76.7152460970319,
              72.4145898678498, 72.5493084772589, 62.8679342999809,
              84.6054837175178, 68.4106926808367, 73.461317178893,
              73.8850009173147, 75.423654885015, 65.4430696634311,
              70.0638707675187, 61.4513622846392, 67.2727559612141,
              66.5002852075071, 70.6822348350969, 69.0458624013487,
              63.3891546006475, 69.5744360848198, 64.3930026414534,
              66.9319603486196, 77.475264909517, 70.8350327032835,
              75.4981394396915, 63.2384878114452, 64.4332862166375,
              68.7133890118472, 89.9834771975547, 79.0187447869306,
              86.7597211883153, 72.0815385539171, 82.6048631223133,
              70.9185524053787, 73.8410355744976), 
       a = c(0.686119175619573, 1.73150008968541, 1.6689814839558,
             1.39594844864604, 1.45724072113527, 1.13454133726276,
             0.654702205998597, 1.12751745781316, 0.572811288958486,
             0.616882582872697, 1.7900438872287, 1.1651297961472,
             1.84692875918164, 1.53366335254838, 1.68600606030997,
             0.71783322195996, 1.74922045523059, 1.37475258900425,
             0.770151879592033, 0.846515935775951, 1.54608301668143,
             1.43679112714864, 1.24138776950173, 1.67010212704258,
             0.717615253638883, 1.12510030643743, 2.02975468096652,
             0.761376001220897, 0.548861821809462, 0.749191300277057,
             0.814536308473333, 0.898742974732295, 0.883807790810071), 
       b = c(0.00708553041975079, 0.00674873980904373, 0.00673529584125957,
             0.00326845452028359, 0.00671523334790717, 0.00337809913819421,
             0.00705426400899282, 0.00684542471161186, 0.00350993343108102,
             0.0035263086672735, 0.00420037752769893, 0.0033698725232357,
             0.00435429978569428, 0.00343665930009632, 0.00661576143592767,
             0.00469812567319346, 0.00324235211431031, 0.00678727254448242,
             0.00472896833548091, 0.00350544397040905, 0.00669614252982428,
             0.0033215200352989, 0.00334206266594306, 0.00317111404779269,
             0.00463690836830452, 0.00462275226303281, 0.00651872747189751,
             0.00463234898288566, 0.00697230409372073, 0.00453026231690821,
             0.00459438375221039, 0.00452613510108387, 0.00463655694555551), 
       d = c(0.464618960352229, 0.260093246987408, 0.274471190461056,
             0.173866197712127, 0.251031671567515, 0.235849502696832,
             0.442831621814862, 0.328338519846985, 0.342629779714319,
             0.351631479194327, 0.162669184630174, 0.241725190500931,
             0.181028132427466, 0.219937737941938, 0.251138701327526,
             0.389425471536826, 0.163746688302077, 0.316281615813166,
             0.390987576993547, 0.325236073913584, 0.292189798277939,
             0.163331452489704, 0.193502318971726, 0.146376628371237,
             0.353293420989731, 0.317843428301467, 0.173025075691225,
             0.308885179700544, 0.407644226988411, 0.285019687702588,
             0.309752611242941, 0.284747934682996, 0.325530583506462), 
       e = c(0.438362806688629, 0.511611268591322, 0.507076569335002,
             0.517161863411317, 0.51335439984887, 0.480575976265248,
             0.442707489861248, 0.482826712938625, 0.438238622135624,
             0.435061433328426, 0.550925438969628, 0.482410949818365,
             0.533203356096027, 0.488197732334187, 0.521591146603833,
             0.44031044264272, 0.530100508479616, 0.492311304626231,
             0.4388998967175, 0.446987373338232, 0.503847566433139,
             0.517266017954479, 0.500619813707586, 0.541949281455506,
             0.451177663259607, 0.465481972619024, 0.555448366631774,
             0.45805224555244, 0.452342212706244, 0.473753346235387,
             0.463159983851921, 0.476484347184577, 0.458117105960893))
  return(built_models)
}

#' Plots a consensus cutoff model
#'
#' The following model is assumed: y = a + bx + dx^e.
#'
#' @param a The a parameter for the model
#' @param b The b parameter for the model
#' @param d The d parameter for the model
#' @param e The e parameter for the model
#' @export

plot_consensus_cutoff_model <- function(a, b, d, e){
  f <- function(x, a, b, d, e){a + x*b + d*(x^e)}
  dat <- data.frame(x = 10:2000)
  dat$y <- f(10:2000, a, b, d, e)
  ggplot(dat, aes(x = x, y = y)) + geom_line()
}
