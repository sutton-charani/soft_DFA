# Functions library for Detrended Distributional Analysis (DDA)
library(nonlinearTseries); library(ggplot2); library(jpeg); library(ggpubr); library(ggdark)
library(seewave); library(dtw)

url <- "https://raw.githubusercontent.com/sutton-charani/uncertain_DFA/main/code/www/tabeau.jpg"
download.file(url, destfile = "tableau.jpg")
tableau <- readJPEG("tableau.jpg")
###############
tab_plot <- function(x, y, title=NULL){
  df <- data.frame(x, y)
  ggplot(df, aes(x, y)) +
    background_image(tableau) +
    geom_line() +
    dark_theme_light() +
    ggtitle(title) +
    theme(text = element_text(size = 20),
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
  kl_dist <- c()
  for (istart in 1 : eval(length(sig) - pattern_size)){
    wind <- sig[istart : eval(istart + pattern_size-1)]
    kl_dist <- c(kl_dist, kl.dist(rep(1/pattern_size, pattern_size), wind)$D)
  }
  main_pattern <- sig[which.min(kl_dist) : eval(which.min(kl_dist) + pattern_size - 1)]
  
  return(main_pattern)
}
###############
dda_complexity <- function(sig, pattern_size){
  main_pattern <- find_main_pattern(sig, pattern_size)
  dtw_dist <- c()
  for (n in 1 : eval(floor(length(sig)/pattern_size))){
    wind <- sig[1 : eval(n * pattern_size)]
    dtw_dist <- c(dtw_dist, dtw(main_pattern, wind)$distance)
  }
  return(sum(dtw_dist))
}
###############

###############

###############









