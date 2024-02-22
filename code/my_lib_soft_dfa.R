
# devtools::source_url(paste0("https://raw.githubusercontent.com/sutton-charani/",
#                             "possibilistic_linear_regression/",
#                             "main/code/my_lib_possibilistic_linear_regression.R"))
############################################################################
dfa_complexity <- function(x, window_growing_factor=2, n_wind=10, doplot=F, size=25){
  
  N <- length(x)
  
  # Detrending
  X <- cumsum(x - mean(x))
  
  # Windowing
  terms <- 10
  while (tail(terms, 1) < N){
    new_term <- ceiling(window_growing_factor * tail(terms, 1))
    terms <- c(terms, new_term)
  }
  terms <- terms[-length(terms)]
  
  fluctuations <- c()
  term=terms[1]
  for (term in terms){
    fluctuations_1term <- c()
    for (iwindow in 1 : n_wind){
      start <- sample(N-term, 1)
      stop <- start + term - 1
      window <- X[start:stop]
      reg_model <- lm(y ~ x, data.frame(x=1:term, y=window))
      fluctuations_1term <- c(fluctuations_1term, sqrt(mean(reg_model$residuals^2)))
    }
    fluctuations <- rbind(fluctuations, data.frame(term=rep(term, n_wind), fluc=fluctuations_1term))
  }
  terms <- unique(fluctuations$term)
  fluctuations_av <- data.frame()
  for (term in terms){
    fluctuations1term=fluctuations[fluctuations$term == term, ]
    fluctuations_av <- rbind(fluctuations_av, data.frame(term=term, fluc=mean(fluctuations1term$fluc)))
    }
  log_reg_model <- lm(fluc ~ term, log(fluctuations))
  alpha <- log_reg_model$coefficients[['term']]
  
  if (doplot){
    p <- ggplot(log(fluctuations), aes(term, fluc)) + 
      geom_point()+
      geom_abline(intercept=log_reg_model$coefficients[['(Intercept)']], slope=alpha, col="red", linewidth=1) +
      xlab("log(term)") + ylab("log(fluctuation)") +
      theme_bw() + 
      theme(text = element_text(size = size)) 
  }
  
  return(list(alpha=alpha, plot=p))
}
###########################################################################
soft_dfa <- function(x, window_growing_factor=2, n_wind=10, doplot=F){
  
  N <- length(x)
  
  # Detrending
  X <- cumsum(x - mean(x))
  
  # Windowing
  terms <- 10
  while (tail(terms, 1) < N){
    new_term <- ceiling(window_growing_factor * tail(terms, 1))
    terms <- c(terms, new_term)
  }
  terms <- terms[-length(terms)]
  
  fluctuations <- c()
  term=terms[1]
  for (term in terms){
    fluctuations_1term <- c()
    for (iwindow in 1 : n_wind){
      start <- sample(N-term, 1)
      stop <- start + term - 1
      window <- X[start:stop]
      reg_model <- lm(y ~ x, data.frame(x=1:term, y=window))
      fluctuations_1term <- c(fluctuations_1term, sqrt(mean(reg_model$residuals^2)))
    }
    fluctuations <- rbind(fluctuations, data.frame(term=rep(term, n_wind), fluc=fluctuations_1term))
  }
  
  log_df <- log(fluctuations)
  soft_reg <- possibilistic_linear_regression(x=log_df$term, y=log_df$fluc, do_plot=T)
  alpha_evid <- soft_reg$slope_possibility
  
  if (doplot){
    p <- soft_reg$plot + xlab("log(term)") + ylab("log(fluctuation)")
  } else {
    p <- NULL
  }
  
  return(list(alpha_evid, plot=p))
}
