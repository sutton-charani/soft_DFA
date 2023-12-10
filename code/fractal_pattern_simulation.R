
# Variables and screen cleaning
graphics.off(); cat("\014"); rm(list=ls()) ;options(warn=-1)


source("/Users/nsc/Google Drive/Recherche/Projets/Escarres/code/my_lib_dda.R")
#devtools::source_url("https://raw.githubusercontent.com/sutton-charani/uncertain_DFA/main/code/my_lib_dda.R")

pattern <- sample(1:10, size=7)# c(2, 8, 7, 1)

p_pattern <- tab_plot(1:length(pattern), pattern, title="Input pattern")

n_convol <- 1
conv_sig <- sim_fractal_signal(sig=pattern, pattern=pattern, n_convol)
dfa_alpha <- as.numeric(estimate(dfa(conv_sig, do.plot=F)))
mdfa_alpha <- mean(MFDFA(conv_sig, scale=10:100, q=-10:10)$spec$hq)
dda_alpha <- dda_complexity(conv_sig, length(pattern))
p_conv_sig <- tab_plot(1:length(conv_sig), conv_sig, 
                       title=paste0(n_convol, "-Convoluated signal \n-> (DFA-MDFA-DDA) alpha = (", 
                                    round(dfa_alpha, 2), ", ", 
                                    #round(mdfa_alpha, 2), ", ", 
                                    round(dda_alpha), ")"),
                       txt_size=15)
ggarrange(p_pattern, p_conv_sig)

n_convol <- 5
scale <- 10:100
q <- -10:10
sig <- pattern

vars <- c()
entropies <- c()
dfa_complexities <- c()
mdfa_complexities <- c()
dda_complexities <- c()
for (i_convol in 1 : n_convol){
  cat(i_convol, " ")
  
  sig <- sim_fractal_signal(sig, pattern, 1)
  
  vars <- c(vars, var(sig))
  entropies <- c(entropies, entropy(sig))
  dfa_alpha <- as.numeric(estimate(dfa(sig, do.plot=F)))
  # if (length(sig) < 200){
  #   mdfa_alpha <- NA
  # } else {
  #   mdfa_alpha <- mean(MFDFA(sig, scale=scale, q=q)$spec$hq)
  # }
  dda_alpha <- dda_complexity(sig, length(pattern))
  dfa_complexities <- c(dfa_complexities, dfa_alpha)
  #mdfa_complexities <- c(mdfa_complexities, mdfa_alpha)
  dda_complexities <- c(dda_complexities, dda_alpha)
}
p_conv_sig <- tab_plot(1:length(sig), sig, 
                       title=paste0(n_convol, "-Convoluated signal \n-> (DFA-MDFA-DDA) alpha = (", 
                                    round(dfa_alpha, 2), ", ", 
                                    #round(mdfa_alpha, 2), ", ", 
                                    round(dda_alpha), 
                                    ")"),
                       txt_size=10)
p_var_complexity <- tab_plot(1:length(vars), vars, title="Variance evolution", txt_size=10)
p_entropy_complexity <- tab_plot(1:length(entropies), entropies, title="Entropy evolution", txt_size=10)
p_dfa_complexity <- tab_plot(1:length(dfa_complexities), dfa_complexities, title="DFA complexity evolution", txt_size=10)
#p_mdfa_complexity <- tab_plot(1:length(mdfa_complexities), mdfa_complexities, title="MDFA complexity evolution", txt_size=10)
p_dda_complexity <- tab_plot(1:length(dda_complexities), dda_complexities, title="DDA complexity evolution", txt_size=10)
ggarrange(p_pattern, p_conv_sig, p_var_complexity, 
          p_entropy_complexity, p_dfa_complexity, 
          #p_mdfa_complexity, 
          p_dda_complexity, 
          nrow=2, ncol=3)















