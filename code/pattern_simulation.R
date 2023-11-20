
# Variables and screen cleaning
graphics.off(); cat("\014"); rm(list=ls()) ;options(warn=-1)

library(nonlinearTseries); library(ggplot2); library(ggpubr)

set.seed(123)
sig <- rnorm(100)
nbreaks <- 5
h <- hist(sig, breaks=nbreaks, plot=F)
pattern <- h$density

p_initial_sig <- ggplot(data.frame(x=1:length(sig), 
                                   y=sig), aes(x,y)) + 
  geom_line() +
  xlab("") + ylab("") +
  ggtitle("initial signal") +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), 
        axis.text.x=element_blank())
p_initial <- ggplot(data.frame(x=h$mids, 
                               y=pattern), aes(x,y)) + 
  geom_line() +
  xlab("") + ylab("") +
  ggtitle("initial pattern") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y=element_blank())
p_initial

fractal_factor <- 5
n <- length(pattern)

final_sig <- pattern
alphas <- c()
sigs <- list()
sigs[[1]] <- pattern
plots <- list()
plots[[1]] <- p_initial
i_window=4
n_convol <- 10
for (i_convol in 1:n_convol){
  new_sig <- final_sig * 1/fractal_factor + 
    approx(x=1:length(final_sig), 
           y=final_sig, 
           n=fractal_factor*length(final_sig))$y
  new_sig <- new_sig * (max(pattern)/max(new_sig))
  # sigs[[i_convol+1]] <- new_sig
  final_sig <- new_sig
  alpha <- as.numeric(estimate(dfa(final_sig, do.plot=F)))
  cat(alpha, " ")
  plots[[i_convol+1]] <- ggplot(data.frame(x=1:length(final_sig), y=final_sig), aes(x,y)) +
    geom_line() +
    ggtitle(paste(i_convol, "convolutions -> alpha =", round(alpha, 3)))+
    theme_bw()
  alphas <- c(alphas, alpha)
}

ggarrange(plots[[1]], plots[[2]], 
          plots[[3]], plots[[4]],
          nrow=2, ncol=2)

ggarrange(plots[[1]], plots[[2]], plots[[3]], 
          plots[[4]], plots[[5]], plots[[6]], 
          plots[[7]], plots[[8]], plots[[9]], 
          nrow=3, ncol=3)

ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
          plots[[5]], plots[[6]], plots[[7]], plots[[8]], 
          plots[[9]], plots[[10]], plots[[11]], NULL,
          nrow=3, ncol=4)











