
devtools::source_url(paste0("https://raw.githubusercontent.com/sutton-charani/",
                            "possibilistic_linear_regression/",
                            "main/code/my_lib_possibilistic_linear_regression.R"),
                     sha1="da5ae5bb22fae50aa2f83c691e594117701dbd94")

dfa_complexity <- function(x, window_growing_factor=2, n_wind=10, doplot=F){
  
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
      #fluctuations_1term <- c(fluctuations_1term, sqrt(mean((window - mean(window))^2)))
    }
    fluctuations <- rbind(fluctuations, data.frame(term=rep(term, n_wind), fluc=fluctuations_1term))
  }
  
  log_reg_model <- lm(fluc ~ term, log(fluctuations))
  alpha <- log_reg_model$coefficients[['term']]
  
  if (doplot){
    tab_plot(x=log(fluctuations$term), y=log(fluctuations$fluc), type="points", 
             x_lab="log(term)", y_lab="log(fluctuation)") +
      geom_smooth(method=lm, se=FALSE, col='red')
  }
  
  return(alpha)
}
############################################################################
soft_dfa_complexity <- function(x, window_growing_factor=2, n_wind=10, doplot=F){
  
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
      #fluctuations_1term <- c(fluctuations_1term, sqrt(mean((window - mean(window))^2)))
    }
    fluctuations <- rbind(fluctuations, data.frame(term=rep(term, n_wind), fluc=fluctuations_1term))
  }
  
  if (!is.null(fluctuations)){
    log_reg_model <- lm(fluc ~ term, log(fluctuations))
    alpha <- log_reg_model$coefficients[['term']]
    soft_log_reg_model <- imprecise_regression(x=log(fluctuations$term), y=log(fluctuations$fluc))
    alpha_min <- soft_log_reg_model$slope_min
    alpha_max <- soft_log_reg_model$slope_max
    if (length(soft_log_reg_model$slope_min) == 0){
      alpha <- NA
      alpha_min <- NA
      alpha_max <- NA
    }
  } else {
    alpha <- NA
    alpha_min <- NA
    alpha_max <- NA
  }
  
  if (doplot){
    tab_plot(x=log(fluctuations$term), y=log(fluctuations$fluc), type="points", 
             x_lab="log(term)", y_lab="log(fluctuation)") +
      geom_smooth(method=lm, se=FALSE, col='red')
  }
  
  return(list(alpha=alpha, alpha_min=alpha_min, alpha_max=alpha_max))
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
      #fluctuations_1term <- c(fluctuations_1term, sqrt(mean((window - mean(window))^2)))
    }
    fluctuations <- rbind(fluctuations, data.frame(term=rep(term, n_wind), fluc=fluctuations_1term))
  }
  
  #log_reg_model <- lm(fluc ~ term, log(fluctuations))
  #alpha <- log_reg_model$coefficients[['term']]
  log_df <- log(fluctuations)
  alpha_possib <- possibilistic_linear_regression(x=log_df$term, y=log_df$fluc)$slope_possibility
  
  if (doplot){
    tab_plot(x=log(fluctuations$term), y=log(fluctuations$fluc), type="points", 
             x_lab="log(term)", y_lab="log(fluctuation)") +
      geom_smooth(method=lm, se=FALSE, col='red')
  }
  
  return(alpha_possib)
}
