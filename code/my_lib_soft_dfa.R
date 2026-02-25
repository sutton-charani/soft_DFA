
library(ggplot2); library(reshape2)
################################
geom_evid_band <- function(intercept_conf, intercept_min, intercept_max, slope_conf, slope_min, slope_max, dataframe, 
                           band_slope=T, band_diags=T, band_border_type="dashed", band_border_width=1, alpha=0.25, 
                           band_diag_col="deepskyblue4", band_border_col="purple"){
  band_color <- "black"
  x <- dataframe$x
  if (band_slope & band_diags){
    list(geom_abline(intercept=intercept_conf, slope=slope_conf, color='red'),
         # diag - slope min
         geom_segment(aes(x=min(x), y=intercept_max + 2*slope_conf*min(x),
                          xend=max(x), yend=intercept_min + slope_conf*(max(x)+min(x))), # slope_conf*min(x) + intercept_max + slope_conf*x)
                      col=band_diag_col),
         # diag - slope max
         geom_segment(aes(x=min(x), y=intercept_min + 2*slope_conf*min(x),
                          xend=max(x), yend=intercept_max + slope_conf*(max(x)+min(x))),
                      col=band_diag_col),
         # left
         geom_segment(aes(x=min(x), y=intercept_min + 2*slope_conf*min(x),
                          xend=min(x), yend=intercept_max + 2*slope_conf*min(x)),
                      linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
         # right
         geom_segment(aes(x=max(x), y=intercept_min + slope_conf*(max(x)+min(x)),
                          xend=max(x), yend=intercept_max + slope_conf*(max(x)+min(x))),
                      linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
         # top
         geom_segment(aes(x=min(x), y=intercept_max + 2*slope_conf*min(x),
                          xend=max(x), yend=intercept_max + slope_conf*(max(x)+min(x))),
                      linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
         # bottom
         geom_segment(aes(x=min(x), y=intercept_min + 2*slope_conf*min(x),
                          xend=max(x), yend=intercept_min + slope_conf*(max(x)+min(x))),
                      linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
         # band
         geom_ribbon(aes(ymin=slope_conf*min(x) + intercept_min + slope_conf*x, 
                         ymax=slope_conf*min(x) + intercept_max + slope_conf*x),
                     fill=band_color, alpha=alpha))
  } else {
    if (band_slope & !band_diags){
      list(geom_abline(intercept=intercept_conf, slope=slope_conf, color='red'),
           # left
           geom_segment(aes(x=min(x), y=intercept_min + 2*slope_conf*min(x),
                            xend=min(x), yend=intercept_max + 2*slope_conf*min(x)),
                        linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
           # right
           geom_segment(aes(x=max(x), y=intercept_min + slope_conf*(max(x)+min(x)),
                            xend=max(x), yend=intercept_max + slope_conf*(max(x)+min(x))),
                        linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
           # top
           geom_segment(aes(x=min(x), y=intercept_max + 2*slope_conf*min(x),
                            xend=max(x), yend=intercept_max + slope_conf*(max(x)+min(x))),
                        linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
           # bottom
           geom_segment(aes(x=min(x), y=intercept_min + 2*slope_conf*min(x),
                            xend=max(x), yend=intercept_min + slope_conf*(max(x)+min(x))),
                        linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
           # band
           geom_ribbon(aes(ymin=slope_conf*min(x) + intercept_min + slope_conf*x, 
                           ymax=slope_conf*min(x) + intercept_max + slope_conf*x),
                       fill=band_color, alpha=alpha))
    } else {
      if (!band_slope & band_diags){
        list(# diag - slope min
          geom_segment(aes(x=min(x), y=intercept_max + 2*slope_conf*min(x),
                           xend=max(x), yend=intercept_min + slope_conf*(max(x)+min(x))), # slope_conf*min(x) + intercept_max + slope_conf*x)
                       col=band_diag_col),
          # diag - slope max
          geom_segment(aes(x=min(x), y=intercept_min + 2*slope_conf*min(x),
                           xend=max(x), yend=intercept_max + slope_conf*(max(x)+min(x))),
                       col=band_diag_col),
          # left
          geom_segment(aes(x=min(x), y=intercept_min + 2*slope_conf*min(x),
                           xend=min(x), yend=intercept_max + 2*slope_conf*min(x)),
                       linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
          # right
          geom_segment(aes(x=max(x), y=intercept_min + slope_conf*(max(x)+min(x)),
                           xend=max(x), yend=intercept_max + slope_conf*(max(x)+min(x))),
                       linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
          # top
          geom_segment(aes(x=min(x), y=intercept_max + 2*slope_conf*min(x),
                           xend=max(x), yend=intercept_max + slope_conf*(max(x)+min(x))),
                       linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
          # bottom
          geom_segment(aes(x=min(x), y=intercept_min + 2*slope_conf*min(x),
                           xend=max(x), yend=intercept_min + slope_conf*(max(x)+min(x))),
                       linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
          # band
          geom_ribbon(aes(ymin=slope_conf*min(x) + intercept_min + slope_conf*x, 
                          ymax=slope_conf*min(x) + intercept_max + slope_conf*x),
                      fill=band_color, alpha=alpha))
      } else { # case (!slope & !diags)
        list(
          # left
          geom_segment(aes(x=min(x), y=intercept_min + 2*slope_conf*min(x),
                           xend=min(x), yend=intercept_max + 2*slope_conf*min(x)),
                       linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
          # right
          geom_segment(aes(x=max(x), y=intercept_min + slope_conf*(max(x)+min(x)),
                           xend=max(x), yend=intercept_max + slope_conf*(max(x)+min(x))),
                       linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
          # top
          geom_segment(aes(x=min(x), y=intercept_max + 2*slope_conf*min(x),
                           xend=max(x), yend=intercept_max + slope_conf*(max(x)+min(x))),
                       linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
          # bottom
          geom_segment(aes(x=min(x), y=intercept_min + 2*slope_conf*min(x),
                           xend=max(x), yend=intercept_min + slope_conf*(max(x)+min(x))),
                       linetype=band_border_type, linewidth=band_border_width, col=band_border_col),
          # band
          geom_ribbon(aes(ymin=slope_conf*min(x) + intercept_min + slope_conf*x, 
                          ymax=slope_conf*min(x) + intercept_max + slope_conf*x),
                      fill=band_color, alpha=alpha)
        )
      }
    }
  }
} 
################################
empirical_conf_int <- function(x, y, confidence=0.95, do_plot=F, band_slope=T, band_diags=T, points_size=10,
                               band_border_type="dashed", size=25, band_border_width=1, band_diag_col="deepskyblue4", 
                               band_border_col="purple"){
  dataframe <- data.frame(x=x, y=y)
  reg_model <- lm(y~x, dataframe)
  slope <- reg_model$coefficients[['x']]
  intercept <- reg_model$coefficients[['(Intercept)']]
  
  err_top <- max(y - predict(reg_model, dataframe))
  err_bottom <- max(predict(reg_model, dataframe) - y)
  
  dataframe$d <- abs(dataframe$y - slope*dataframe$x - intercept)
  dataframe <- dataframe[order(dataframe$d, decreasing=F), ]
  irow_max <- floor(confidence * nrow(dataframe))
  df_conf <- dataframe[1 : irow_max, ]
  reg_model_conf <- lm(y~x, df_conf)
  slope_conf <- reg_model_conf$coefficients[['x']]
  intercept_conf <- reg_model_conf$coefficients[['(Intercept)']]
  
  intercept_min <- intercept_conf - max(df_conf$d)
  intercept_max <- intercept_conf + max(df_conf$d)
  slope_min <- (slope_conf*(max(df_conf$x) - min(df_conf$x)) + intercept_min - intercept_max) / 
    (max(df_conf$x) - min(df_conf$x))
  slope_max <- (slope_conf*(min(df_conf$x) - max(df_conf$x)) + intercept_min - intercept_max) / 
    (min(df_conf$x) - max(df_conf$x))
  
  if (is.na(slope_conf)){ # this can happen when x is constant
    slope_conf <- mean(y)
    slope_min <- slope_conf
    slope_max <- slope_conf
  }
  
  if (do_plot){
    p <- ggplot(dataframe, aes(x, y)) + 
      #geom_point(size=points_size) +
      theme_bw() +
      theme(text=element_text(size=size), plot.title=element_text(hjust=0.5)) +
      geom_evid_band(intercept_conf, intercept_min, intercept_max, slope_conf, slope_min, slope_max, dataframe, 
                     band_slope=band_slope, band_diags=band_diags, band_border_type=band_border_type, 
                     band_border_width=band_border_width, band_diag_col=band_diag_col, band_border_col=band_border_col)
  } else {
    p <- NULL
  }
  
  
  result <- list(intercept=intercept_conf, intercept_min=intercept_min,intercept_max=intercept_max, 
                 slope=slope_conf, slope_min=slope_min,slope_max=slope_max, plot=p)
  return(result)
}
################################
soft_linear_regression <- function(x, y, initial_confidences=c(0.1, 0.5, 0.95), do_plot=F, size=25, show_slope=F, 
                                   show_diags=F, alpha_fac=1, band_border_type=NA, band_border_width=1, 
                                   band_diag_col="deepskyblue4", band_border_col="purple", alpha=25, points_size=10){
  
  precise_reg_model <- lm(y~x, data.frame(x=x, y=y))
  precise_slope <- precise_reg_model$coefficients[['x']]
  precise_intercept <- precise_reg_model$coefficients[['(Intercept)']]
  
  intercept_interval <- data.frame()
  slope_interval <- data.frame()
  slope_confs <- c()
  intercept_confs <- c()
  for (confidence in initial_confidences){
    soft_regression_intervals <- empirical_conf_int(x, y, confidence=confidence, points_size=points_size, do_plot=do_plot)
    intercept_interval <- rbind(intercept_interval, unlist(soft_regression_intervals[c('intercept_min', 'intercept_max')]))
    slope_interval <- rbind(slope_interval, unlist(soft_regression_intervals[c('slope_min', 'slope_max')]))
    slope_confs <- c(slope_confs, soft_regression_intervals$slope)
    intercept_confs <- c(intercept_confs, soft_regression_intervals$intercept)
  }
  names(slope_interval) <- c('slope_min', 'slope_max')
  
  #slope_pignistic_expectation <- sum(initial_confidences/sum(initial_confidences) * colMeans(slope_interval))
  slope_intervals_means <- rowMeans(slope_interval)
  slope_pignistic_expectation <- sum(slope_intervals_means * initial_confidences) / sum(initial_confidences)
  slope_min_content <- sum(initial_confidences/sum(initial_confidences) * slope_interval$slope_min)
  slope_max_contour <- sum(initial_confidences/sum(initial_confidences) * slope_interval$slope_max)
  
  intercept_interval <- rbind(intercept_interval, c(-Inf, +Inf))
  slope_interval <- rbind(slope_interval, c(-Inf, +Inf))
  mass_omega <- 1/(length(x) - 1) # Shafer uncertainty model
  confidences <- initial_confidences/(sum(initial_confidences) + mass_omega)
  names(intercept_interval) <- c("intercept_min", "intercept_max")
  names(slope_interval) <- c("slope_min", "slope_max")
  intercept_interval$mass <- c(confidences, mass_omega)
  slope_interval$mass <- c(confidences, mass_omega)
  
  if (do_plot){
    dataframe <- data.frame(x=x, y=y)
    reg_model <- lm(y~x, dataframe)
    slope <- reg_model$coefficients[['x']]
    intercept <- reg_model$coefficients[['(Intercept)']]
    
    p <- ggplot(dataframe, aes(x, y)) + 
      theme_bw() + xlab("x") + ylab("y") +
      theme(text=element_text(size=size),
            plot.title=element_text(hjust=0.5))
    for (i in seq(from=length(initial_confidences), to=1, by=-1)){
      command_line <- paste0("p <- p + geom_evid_band(intercept_confs[", i,
                             "], intercept_interval$intercept_min[", i,
                             "], intercept_interval$intercept_max[", i,
                             "], slope_confs[", i,
                             "], slope_interval$slope_min[", i,
                             "], slope_interval$slope_max[", i,
                             "], dataframe, band_slope=F, band_diags=show_diags, band_border_type=band_border_type, alpha=alpha_fac * (1 - e1071::sigmoid(1-confidences[", i,
                             "])), band_border_width=band_border_width, band_diag_col=band_diag_col, band_border_col=band_border_col)")
      # alpha=alpha_fac * (1-confidences[", i, "]) -> alpha=alpha_fac * (1 - e1071::sigmoid(confidences[", i, ?
      eval(parse(text=command_line))
    }
    #p <- p + geom_point(size=points_size)
    if (show_slope){
      p <- p + geom_abline(intercept=intercept, slope=slope, col="red", linewidth=1)
    }
  } else {
    p <- NULL
  }
  
  result <- list(precise_slope=precise_slope, precise_intercept=precise_intercept, 
                 slope_possibility=slope_interval, intercept_possibility=intercept_interval, slope_confs=slope_confs, 
                 slope_EbetP=slope_pignistic_expectation, slope_min_content=slope_min_content, slope_max_contour=slope_max_contour,
                 plot=p)
  return(result)
}
################################
compute_window_fluctuations <- function(window, fluc_type="dfa"){
  
  if (fluc_type == "dfa"){
    if (is.null(dim(window))){ # if time series
      reg_model <- lm(y ~ x, data.frame(x=1:length(window), y=window))
      fluc <- sqrt(mean((predict(reg_model)-window)^2))
    } else { # image
      melted_win <- reshape2::melt(window)
      n <- nrow(melted_win)
      reg_model <- lm(value ~ ., melted_win)
      melted_predicted_win <- predict(reg_model, melted_win)
      fluc <- sqrt(mean((melted_win$value - melted_predicted_win)^2))
    }
  } else if (fluc_type == "hurst"){
    R <- max(window) - min(window)
    S <- sd(window)
    fluc <- R/S
  } else {
    fluc <- var(window)
  }
  
  return(fluc)
}
################################
compute_fluctuations <- function(x, n_wind=10, n_sizes=20, min_win_size=10, max_win_size=length(x)/2, 
                                 #x_limit_sizes=c(10, ncol(x)/2), y_limit_sizes=c(10, nrow(x)/2), 
                                 min_win_side=10, max_win_side=round(nrow(x)/2), n_wind_per_surface=100,
                                 windowing="uniform", fluc_type="dfa", verbose=F){
  # n_wind=10; n_sizes=20; min_win_size=10; max_win_size=length(x)/2; fluc_type="dfa"; windowing="uniform"
  
  if (is.null(dim(x))){ # if time series
    
    # Profile computation
    X <- cumsum(x - mean(x))
    n <- length(X)
    
    # Windowing
    if (windowing == "complete"){
      max_win_size <- n
      window_sizes <- seq(min_win_size, max_win_size)
    } else {
      window_sizes <- round(seq(min_win_size, max_win_size, length.out=n_sizes))
    }
    
    fluctuations <- c()
    window_size=window_sizes[1]
    for (window_size in window_sizes){
      fluctuations_1window_size <- c()
      iwindow=1
      if (windowing == "uniform"){
        starts <- round(seq(1, n - window_size, length.out=n_wind))
      }  else if (windowing == "complete"){
        starts <- seq(1, n - window_size + 1)
        n_wind <- length(starts)
      }
      iwindow=1
      for (iwindow in 1 : n_wind){
        if (windowing == "random"){
          start <- sample(n - window_size + 1, 1)
        } else { # (windowing == "uniform"){
          start <- starts[iwindow]
        }
        stop <- start + window_size - 1
        window <- X[start:stop]
        fluc <- compute_window_fluctuations(window)
        fluctuations_1window_size  <- c(fluctuations_1window_size, fluc)
      }
      fluctuations <- rbind(fluctuations, data.frame(window_size=rep(window_size, n_wind), fluc=fluctuations_1window_size))
    }
  } else { # image
    
    im <- x
    
    if (length(dim(im)) > 2){
      im <- rgb_2gray(im)
    }
    
    x_dim <- ncol(im)
    y_dim <- nrow(im)
    if (verbose) cat('y_dim =', y_dim)
    melted_im <- reshape2::melt(im)
    N <- nrow(melted_im)
    
    # Profile computation
    if (verbose) cat('\nprofile computation')
    mu <- mean(im)
    IM <- matrix(NA, nrow=y_dim, ncol=x_dim)
    if (verbose)cat('\ni = ')
    for (i in 1:y_dim){
      if (verbose)cat(i, ' ')
      for (j in 1:x_dim){
        IM[i, j] <- sum(im[1:i, 1:j] - mu)
      }
    }
    
    # Adatation of windows limit sizes to the image dimension
    win_sides <- round(seq(from=min_win_side, to=max_win_side, length.out=n_sizes))
    
    fluctuations <- c()
    # x_win_size=x_win_sizes[1]; y_win_size=y_win_sizes[1]
    win_side=win_sides[1]
    if (verbose) cat("\n", length(win_sides), "window sides to process: ")
    for (win_side in win_sides){
      if (verbose) cat(win_side, ' ')
      
      # Windows centers computation
      if (is.integer(sqrt(n_wind_per_surface))){
        x_centers <- round(seq(from=win_side/2, to=x_dim-win_side/2-1, length.out=sqrt(n_wind_per_surface)))
        y_centers <- round(seq(from=win_side/2, to=y_dim-win_side/2-1, length.out=sqrt(n_wind_per_surface)))
      } else {
        x_centers <- round(seq(from=win_side/2, to=x_dim-win_side/2-1, length.out=floor(sqrt(n_wind_per_surface))))
        y_centers <- round(seq(from=win_side/2, to=y_dim-win_side/2-1, length.out=floor(sqrt(n_wind_per_surface))))
      }
      fluctuations_1win_dim <- c()
      x_center=x_centers[1]; y_center=y_centers[1]#x_center=tail(x_centers, 1); y_center=tail(y_centers, 1)
      for (x_center in x_centers){
        for (y_center in y_centers){
          window <- IM[seq(from=round(y_center - win_side/2) + 1, to=round(y_center + win_side/2)),
                       seq(from=round(x_center - win_side/2) + 1, to=round(x_center + win_side/2))]
          if (sd(window) > 0){
            fluc <- compute_window_fluctuations(window, fluc_type="dfa")
            fluctuations_1win_dim <- c(fluctuations_1win_dim, fluc)
          }
        }
      }
      if (length(fluctuations_1win_dim) > 0){
        fluctuations <- rbind(fluctuations,
                              data.frame(win_side=win_side, fluc=fluctuations_1win_dim))
      }
    }
  }
  
  return(fluctuations)
}
################################
soft_dfa <- function(x, n_wind=10, n_sizes=20, min_win_size=10, max_win_size=300, 
                     initial_confidences=c(0.1, 0.5, 0.95), doplot=F, alpha_fac=1){
  
  fluctuations <- compute_fluctuations(x, n_wind=n_wind, n_sizes=n_sizes,
                                       min_win_size=min_win_size, max_win_size=max_win_size)
  log_df <- log10(fluctuations)
  soft_reg <- soft_linear_regression(x=log_df$term, y=log_df$fluc, 
                                     initial_confidences=initial_confidences, 
                                     do_plot=doplot, alpha_fac=alpha_fac)
  alpha_evid <- soft_reg$slope_possibility
  alpha_evid_EbetP <- soft_reg$slope_EbetP
  alpha_max_contour <- soft_reg$slope_max_contour
  alpha_min_content <- soft_reg$slope_min_content
  
  if (doplot){
    p <- soft_reg$plot + xlab("log(term)") + ylab("log(fluctuation)")
  } else {
    p <- NULL
  }
  
  return(list(fluctuations=fluctuations, alpha_evid=alpha_evid, 
              alpha_EbetP=alpha_evid_EbetP, alpha_min_content=alpha_min_content, alpha_max_contour=alpha_max_contour, plot=p))
}
################################
fluctions_cloud <- function(ts, min_win_size=10, max_win_size=length(ts)-2, n_sizes=1000, n_wind=100,
                            transparency=0.2, points_size=0.1, slope=F, loc_slopes=F, do.plot=F){
  # min_win_size=10; max_win_size=length(ts)-2; n_sizes=1000; n_wind=100; transparency=0.2; points_size=0.1; slope=F; loc_slopes=F; do.plot=F
  fluctuations <- compute_fluctuations(ts, n_sizes=n_sizes, n_wind=n_wind, 
                                       min_win_size=min_win_size, max_win_size=max_win_size)

  log_fluc <- log10(fluctuations)
  names(log_fluc) <- c('log_n', 'log_fluc')
  
  lin_mod <- lm(log_fluc ~ log_n, data=log_fluc)  
  
  p <- ggplot(log_fluc, aes(x=log_n, y=log_fluc)) + 
    geom_point(size=points_size, alpha=transparency, color="white") + 
    theme_void() +
    theme(plot.background=element_rect(fill="#04304B"),
          plot.title=element_text(size=40, colour="white", hjust=0.5),
          plot.margin=unit(c(-1, -1, -1, -1), "cm"))
  if (slope){
    p <- p + geom_abline(intercept=lin_mod$coefficients[['(Intercept)']],
                         slope=lin_mod$coefficients[['log_n']], color='red')
  }
  if (loc_slopes){
    p <- p + stat_summary(geom="point", fun.y="mean", col="green", size=points_size*5) +
      stat_summary(geom="line", fun.y="mean", col="green", size=points_size*3) +
      stat_summary(geom="point", fun.y="median", col="blue", size=points_size*5) +
      stat_summary(geom="line", fun.y="median", col="blue", size=points_size*3)
  }
  if (do.plot){
    p
  }
  return(list(log_fluc=log_fluc, p=p))
}
################################
compute_cloud_fluctuations <- function(x, min_win_size=10, fluc_type="dfa"){
  
  N <- length(x)
  
  # Detrending
  X <- cumsum(x - mean(x))
  
  max_win_size <- N
  
  # Windowing
  window_sizes <- seq(min_win_size, max_win_size)
  #log_window_sizes <- log(window_sizes)
  # terms <- exp(log_window_sizes)
  terms <- window_sizes
  
  fluctuations <- c()
  term=terms[1]
  for (term in terms){
    fluctuations_1term <- c()
    n_wind <- N - term/2
    start=1
    for (start in seq(1, N - term + 1)){
      stop <- start + term - 1
      window <- X[start:stop]
      fluc <- compute_window_fluctuations(window, fluc_type)
      fluctuations_1term <- c(fluctuations_1term, fluc)
    }
    fluctuations <- rbind(fluctuations, 
                          data.frame(term=rep(term, length(fluctuations_1term)), fluc=fluctuations_1term))
  }
  return(fluctuations)
}
################################
plot_cloud_fluctuations <- function(log_fluc, transparency=NULL, points_size=NULL){
  
  n <- nrow(log_fluc)
  if (is.null(points_size)){
    if (n < 30){
      points_size <- 20
    } else if(n < 100){
      points_size <- 8.5
    } else if(n < 200){
      points_size <- 3
    } else{
      points_size <- 0.05
    }
  }
  
  if (is.null(points_size)){
    if (n < 30){
      transparency <- 0.2
    } else {
      transparency <- 0.1
    } 
  }
  
  p <- ggplot(log_fluc, aes(x=log_n, y=log_fluc)) + 
    geom_point(size=points_size, alpha=transparency, color="white") + 
    theme_void() +
    theme(plot.background=element_rect(fill="#04304B"),
          plot.title=element_text(size=40, colour="white", hjust=0.5),
          plot.margin=unit(c(-1, -1, -1, -1), "cm"))
  return(p)
}
################################
compute_fluctuations_im <- function(im, n_sizes=30, n_wind_per_size=100, verbose=F,
                                    x_limit_sizes=c(10, ncol(im)/2), y_limit_sizes=c(10, nrow(im)/2)){ 
  
  if (length(dim(im)) > 2){
    im <- rgb_2gray(im)
  }
  
  x_dim <- ncol(im)
  y_dim <- nrow(im)
  melted_im <- reshape2::melt(im)
  N <- nrow(melted_im)
  
  # Detrending
  melted_im$value <- melted_im$value - mean(melted_im$value)
  
  # Adatation of windows limit sizes to the image dimension
  x_limit_sizes[2] <- min(x_limit_sizes[2], round(x_dim/2))
  y_limit_sizes[2] <- min(y_limit_sizes[2], round(y_dim/2))
  
  # Windowing
  if (is.integer(sqrt(n_sizes))){
    x_win_sizes <- round(seq(from=x_limit_sizes[1], to=x_limit_sizes[2], length.out=sqrt(n_sizes)))
    y_win_sizes <- round(seq(from=y_limit_sizes[1], to=y_limit_sizes[2], length.out=sqrt(n_sizes)))
  } else {
    x_win_sizes <- round(seq(from=x_limit_sizes[1], to=x_limit_sizes[2], length.out=ceiling(sqrt(n_sizes))))
    y_win_sizes <- round(seq(from=y_limit_sizes[1], to=y_limit_sizes[2], length.out=round(sqrt(n_sizes))))
  }
  
  fluctuations <- c()
  x_win_size=x_win_sizes[1]; y_win_size=y_win_sizes[1]
  if (verbose) paste(length(x_win_sizes) * length(y_win_sizes), "window sizes to process")
  for (x_win_size in x_win_sizes){
    if (verbose) cat('\n')
    for (y_win_size in y_win_sizes){
      if (verbose) cat(paste0(", (", x_win_size, ', ', y_win_size, ')'))
      
      # Windows centers computation
      if (is.integer(sqrt(n_wind_per_size))){
        x_centers <- round(seq(from=x_win_size/2, to=x_dim-x_win_size/2-1, length.out=sqrt(n_wind_per_size)))
        y_centers <- round(seq(from=y_win_size/2, to=y_dim-y_win_size/2-1, length.out=sqrt(n_wind_per_size)))
      } else {
        x_centers <- round(seq(from=x_win_size/2, to=x_dim-x_win_size/2-1, length.out=floor(sqrt(n_wind_per_size))))
        y_centers <- round(seq(from=y_win_size/2, to=y_dim-y_win_size/2-1, length.out=floor(sqrt(n_wind_per_size))))
      }
      
      fluctuations_1win_dim <- c()
      x_center=x_centers[1]; y_center=y_centers[1]#x_center=tail(x_centers, 1); y_center=tail(y_centers, 1)
      for (x_center in x_centers){
        for (y_center in y_centers){
          window <- im[seq(from=round(y_center - y_win_size/2) + 1, to=round(y_center + y_win_size/2)),
                       seq(from=round(x_center - x_win_size/2) + 1, to=round(x_center + x_win_size/2))]
          if (sd(window) > 0){
            melted_win <- reshape2::melt(window)
            n <- nrow(melted_win)
            reg_model <- lm(value ~ ., melted_win)
            melted_predicted_win <- predict(reg_model, melted_win)
            predicted_win <- matrix(melted_predicted_win, nrow=nrow(window), ncol=ncol(window), byrow=T)
            fluctuations_1win_dim <- c(fluctuations_1win_dim, sqrt(mean((window - predicted_win)^2)))
          }
        }
      }
      if (length(fluctuations_1win_dim) > 0){
        fluctuations <- rbind(fluctuations,
                              data.frame(fluc=fluctuations_1win_dim, x_size=x_win_size, y_size=y_win_size))
      }
    }
  }
  
  return(fluctuations)
}
################################
spatial_dfa <- function(im, n_sizes=20, n_wind_per_size=30, verbose=F,
                        x_limit_sizes=c(10, ncol(im)/2), y_limit_sizes=c(10, nrow(im)/2)){
  
  # Fluctuations computation
  flucs <- compute_fluctuations_im(im, n_sizes=n_sizes, n_wind_per_size=n_wind_per_size, verbose=F,
                                   x_limit_sizes=x_limit_sizes, y_limit_sizes=y_limit_sizes)
  flucs$win_size <- flucs$x_size * flucs$y_size
  # Average fluctuations computation (for each window size)
  av_flucs <- c()
  for (x_size in unique(flucs$x_size)){
    for (y_size in unique(flucs$y_size)){
      flucs_1win_side <- flucs[flucs$x_size == x_size & flucs$y_size == y_size, ]
      av_flucs <- rbind(av_flucs, colMeans(flucs_1win_side))
    }
  }
  
  # Alpha coefficient computation
  logflucs <- log(flucs)
  lin_mod <- lm(fluc ~ win_size, logflucs)
  alpha <- lin_mod$coefficients['win_size']
  return(list(alpha=alpha, av_flucs=av_flucs))
}


