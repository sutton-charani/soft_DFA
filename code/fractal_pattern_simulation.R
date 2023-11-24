
# Variables and screen cleaning
graphics.off(); cat("\014"); rm(list=ls()) ;options(warn=-1)

library(nonlinearTseries); library(ggplot2); library(ggpubr)
library(ggdark); library(ggthemes)

set.seed(123)
pattern_values <- c(8, 7)
pattern <- c(0, pattern_values, 0)
n <- length(pattern)
n_values <- length(pattern_values)

p_pattern <- ggplot(data.frame(x=1:length(pattern), 
                               y=pattern), aes(x,y)) + 
  geom_line() +
  xlab("") + ylab("") +
  ggtitle("initial pattern") +
  dark_theme_light() +
  theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y=element_blank())
p_pattern

fractal_factor <- 0.2
n_convol <- 8

final_sig <- pattern
alphas <- c()
plots <- list()
plots[[1]] <- p_pattern
i_convol=1
for (i_convol in 1:n_convol){
  interpolated_signal <- #final_sig +#* 1/fractal_factor + 
    approx(x=1:length(final_sig), 
           y=final_sig, 
           n=length(final_sig) + (length(final_sig)-1)*(n-2))$y#fractal_factor*length(final_sig)
  convoluted_signal <- interpolated_signal
  ibreak=1
  for (ibreak in 1:eval(length(final_sig)-1)){
    segt <- convoluted_signal[seq(from=(ibreak-1)*(length(pattern)-1)+1, 
                                  to=(ibreak-1)*(length(pattern)-1)+4)] 
    segt <- segt + pattern * fractal_factor
    convoluted_signal[seq(from=(ibreak-1)*(length(pattern)-1)+1, 
                          to=(ibreak-1)*(length(pattern)-1)+4)] <- segt 
  }
  #new_sig <- new_sig * (max(pattern)/max(new_sig))
  final_sig <- convoluted_signal
  alpha <- as.numeric(estimate(dfa(final_sig, do.plot=F)))
  cat(alpha, " ")
  plots[[i_convol+1]] <- ggplot(data.frame(x=1:length(final_sig), y=final_sig), aes(x,y)) +
    geom_line() +
    xlab("") + ylab("") +
    ggtitle(paste(i_convol, "convolutions -> alpha =", round(alpha, 3)))+
    dark_theme_light() +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.text.y=element_blank())
  alphas <- c(alphas, alpha)
  print(ggplot(data.frame(x=1:length(alphas), 
                    y=alphas), aes(x,y)) + 
    geom_line() +
    xlab("nb convolutions") + ylab("alpha") +
    ggtitle(paste("Complexity evolution for", i_convol, "nested convolutions")) +
      dark_theme_light())
}

p_evol_alpha <- ggplot(data.frame(x=1:length(alphas), 
                                  y=alphas), aes(x,y)) + 
  geom_line() +
  xlab("nb convolutions") + ylab("alpha") +
  ggtitle(paste("Complexity evolution for", i_convol, "nested convolutions")) +
  dark_theme_light()

ggarrange(p_initial_sig, p_initial, plots[[2]], 
          plots[[3]], plots[[4]], plots[[5]], 
          plots[[6]], plots[[7]], p_evol_alpha, 
          nrow=3, ncol=3)

ggarrange(p_initial, plots[[2]], plots[[3]], plots[[4]], plots[[5]], 
          plots[[6]], plots[[7]], plots[[8]], plots[[9]], plots[[10]], 
          plots[[11]], plots[[12]], plots[[13]], plots[[14]], p_evol_alpha,
          nrow=3, ncol=5, main=paste("nbreaks =", nbreaks, 
                                     ", fractal_factor=", fractal_factor))

ggarrange(plots[[1]], plots[[2]], 
          plots[[3]], plots[[4]],
          nrow=2, ncol=2)

ggarrange(plots[[2]], plots[[3]], plots[[4]],
          plots[[5]], plots[[6]], plots[[7]],
          plots[[8]], plots[[9]], plots[[10]],
          nrow=3, ncol=3)

ggarrange(plots[[1]], plots[[2]], plots[[3]], 
          plots[[4]], plots[[5]], plots[[6]], 
          plots[[7]], plots[[8]], p_evol_alpha, 
          nrow=3, ncol=3)

ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
          plots[[5]], plots[[6]], plots[[7]], plots[[8]], 
          plots[[9]], plots[[10]], plots[[11]], NULL,
          nrow=3, ncol=4)











