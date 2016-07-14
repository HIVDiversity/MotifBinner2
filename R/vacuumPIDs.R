explorePIDRelationships <- function(){
  # Probability that that bin1 is offspring of bin2:
  # pid_len = 8:14
  # b1_size = 1:10
  # pid_dist = 1:3
  # b2_size = 20:1000
  # seq_err = 0.01/3 <- must be for specific base!
  # no_seq_err = 0.99

  # conclusions:
  # Consider all bins within 3 pid dist of a >cc bin for offspringiness
  # Allow user to set:
  # minimum bin size (mbs) to consider - default to max(3, trunc(cuttoff/2))
  # consensus cutoff (value or model)
  # remember to use the true bin size and not the observed bin size
  # if a bin > consensus cutoff (cc) then consider it own bin
  # if a bin < cc and bin > mbs and p(offspring) > ? merge into bigger bin
  # if a bin < cc and bin > mbs and p(offspring) < ? - make it its own bin

  seq_err <- 0.01/3
  no_seq_err <- 0.99
  pid_lens <- 8:14
  pid_dists <- 1:3

  results <- list()
  for (b1_size in 1:10){
    print (b1_size)
    for (b2_size_eff in seq(25,10000, 25)){
      tmp_df <- data.frame(b1_size = b1_size,
                           b2_size_eff = b2_size_eff,
                           pid_dist = 0,
                           pid_len = 0,
                           off_spring_prob = rep(0, length(pid_lens)*length(pid_dists))
                           )
      i <- 0
      for (pid_dist in pid_dists){
        for (pid_len in pid_lens){
          success_prob <- (seq_err^pid_dist)*(no_seq_err^(pid_len-pid_dist))
          off_spring_prob <- dbinom(b1_size, b2_size_eff, success_prob)
          i <- i + 1
          tmp_df[i,c('pid_dist', 'pid_len', 'off_spring_prob')] <-
            c(pid_dist, pid_len, off_spring_prob)
        }
      }
      results[[length(results)+1]] <- tmp_df
    }
  }
  result_df = data.table::rbindlist(results)
  result_df$b1_size <- as.factor(result_df$b1_size)

  ggplot(result_df, aes(x = b2_size_eff, y = off_spring_prob,
                        col = b1_size)) +
    geom_line(aes(group = b1_size)) +
    facet_grid(pid_dist ~ pid_len)

  ggplot(subset(result_df, pid_dist == 1 & pid_len == 9 & b2_size_eff < 1000), 
         aes(x = b2_size_eff, y = off_spring_prob,
             col = b1_size)) +
    geom_line(aes(group = b1_size)) +
    facet_grid(pid_dist ~ pid_len)
}
