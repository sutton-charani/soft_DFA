# Functions library for Detrended Distributional Analysis (DDA)
library(nonlinearTseries); library(ggplot2); library(jpeg); library(ggpubr); library(ggdark)
library(seewave); library(dtw); library(MFDFA); library(entropy)

url <- "https://raw.githubusercontent.com/sutton-charani/uncertain_DFA/main/code/www/tabeau.jpg"
download.file(url, destfile = "tableau.jpg")
tableau <- readJPEG("tableau.jpg")
###############
tab_plot <- function(x, y, title=NULL, txt_size=20){
  df <- data.frame(x, y)
  ggplot(df, aes(x, y)) +
    background_image(tableau) +
    geom_line() +
    dark_theme_light() +
    ggtitle(title) +
    theme(text = element_text(size = txt_size),
          plot.title = element_text(hjust = 0.5))
}
###############
sim_fractal_signal <- function(sig=sig, pattern=pattern, n_convol, fractal_factor=0.5){
  n <- length(pattern)
  final_sig <- sig
  for (i_convol in 1:n_convol){
    interpolated_signal <- #final_sig +#* 1/fractal_factor + 
      approx(x=1:length(final_sig), 
             y=final_sig, 
             n=length(final_sig) + (length(final_sig)-1)*(n-2))$y
    convoluted_signal <- interpolated_signal
    ibreak=1
    for (ibreak in 1:eval(length(final_sig)-1)){
      segt <- convoluted_signal[seq(from=(ibreak-1)*(length(pattern)-1)+1, 
                                    to=(ibreak-1)*(length(pattern)-1)+4)] 
      segt <- segt + pattern * fractal_factor
      convoluted_signal[seq(from=(ibreak-1)*(length(pattern)-1)+1, 
                            to=(ibreak-1)*(length(pattern)-1)+4)] <- segt 
    }
    final_sig <- convoluted_signal
  }
  return(final_sig)
}
###############
find_main_pattern <- function(sig, pattern_size){
  # Find the window (starting from the origin) that seems to be the most repeated in the whole signal.
  # 
  # criteria: simarity in terms of dtw and score = 
  # n <- length(sig)
  # scores <- c()
  # for (i in seq(2, n/2)){
  #   wind1 <- sig[1 : i]
  #   score <- 0
  #   for (j in seq(i+1, n-1)){
  #     wind2 <- sig[j : n]
  #     score <- score + dtw(wind1, wind2)$distance
  #   }
  #   scores <- c(scores, score)
  # }
  # main_pattern <- sig[seq(1, which.min(scores))]
  kl_dist <- c()
  for (istart in 1 : eval(length(sig) - pattern_size)){
    wind <- sig[istart : eval(istart + pattern_size-1)]
    kl_dist <- c(kl_dist, kl.dist(rep(1/pattern_size, pattern_size), wind)$D)
  }
  main_pattern <- sig[which.min(kl_dist) : eval(which.min(kl_dist) + pattern_size - 1)]
  
  return(main_pattern)
}
###############
dda_complexity <- function(sig, pattern_size, verbose=F){
  if (verbose){
    cat("main paterrn search...")
  }
  main_pattern <- find_main_pattern(sig, pattern_size)
  if (verbose){
    cat("ok")
  }
  dtw_dist <- c()
  if (verbose){
    cat(length(seq(1, floor(length(sig)/pattern_size))), "window to try...")
  }
  for (n in seq(1, floor(length(sig)/pattern_size))){
    if (verbose){
      cat(n, " ")
    }
    wind <- sig[seq(1, n * pattern_size)]
    dtw_dist <- c(dtw_dist, dtw(main_pattern, wind)$normalizedDistance)
  }
  return(sum(dtw_dist))
}
###############

###############

###############









