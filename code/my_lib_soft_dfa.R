
library(ggplot2); library(reshape2); library(av); library(grid); library(gridExtra)
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
    if (verbose){
      cat(length(window_sizes), 'scales to compute: ')
    }
    for (iwindow_size in 1:length(window_sizes)){
      window_size <- window_sizes[iwindow_size]
      if (verbose){
        cat(iwindow_size, ' ')
      }
      fluctuations_1window_size <- c()
      iwindow=1
      if (windowing == "uniform"){
        starts <- round(seq(1, n - window_size, length.out=n_wind))
      } else if (windowing == "complete"){
        starts <- seq(1, n - window_size + 1)
        n_wind <- length(starts)
      } else if(windowing == "random"){
        starts <- sample(n - window_size + 1, n_wind)
      }
      stops <- starts + window_size - 1
      iwindow=1
      for (iwindow in 1 : n_wind){
        # if (windowing == "random"){
        #   start <- sample(n - window_size + 1, 1)
        # } else { # (windowing == "uniform"){
        #   start <- starts[iwindow]
        # }
        start <- starts[iwindow]
        stop <- stops[iwindow]
        # stop <- start + window_size - 1
        window <- X[start:stop]
        fluc <- compute_window_fluctuations(window)
        fluctuations_1window_size  <- c(fluctuations_1window_size, fluc)
      }
      fluctuations <- rbind(fluctuations, data.frame(window_size=rep(window_size, n_wind), fluc=fluctuations_1window_size, start=starts, stop=stops))
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
    if (verbose) cat("\n", length(win_sides), "scales to process: ")
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
                              data.frame(win_side=win_side, fluc=fluctuations_1win_dim, window_index=win_sides))
      }
    }
  }
  return(fluctuations)
}
################################
soft_dfa <- function(x, n_wind=10, n_sizes=20, min_win_size=10, max_win_size=300, 
                     initial_confidences=c(0.1, 0.5, 0.95), doplot=F, alpha_fac=1,
                     short_long_term_log_threshold=1){
  
  fluctuations <- compute_fluctuations(x, n_wind=n_wind, n_sizes=n_sizes,
                                       min_win_size=min_win_size, max_win_size=max_win_size)
  log_df <- log10(fluctuations)
  log_df_short_term <- log_df[log_df$window_size <= short_long_term_log_threshold, ]
  log_df_long_term <- log_df[log_df$window_size > short_long_term_log_threshold, ]
  
  soft_reg <- soft_linear_regression(x=log_df$window_size, y=log_df$fluc, 
                                     initial_confidences=initial_confidences, 
                                     do_plot=doplot, alpha_fac=alpha_fac)
  soft_reg_short_term <- soft_linear_regression(x=log_df_short_term$window_size, y=log_df_short_term$fluc, 
                                                initial_confidences=initial_confidences, 
                                                do_plot=doplot, alpha_fac=alpha_fac)
  soft_reg_long_term <- soft_linear_regression(x=log_df_long_term$window_size, y=log_df_long_term$fluc, 
                                               initial_confidences=initial_confidences, 
                                               do_plot=doplot, alpha_fac=alpha_fac)
  
  alpha_evid <- soft_reg$slope_possibility
  alpha_evid_EbetP <- soft_reg$slope_EbetP
  alpha_evid_EbetP_short_term  <- soft_reg_short_term $slope_EbetP
  alpha_evid_EbetP_long_term <- soft_reg_long_term$slope_EbetP
  alpha_max_contour <- soft_reg$slope_max_contour
  alpha_min_content <- soft_reg$slope_min_content
  
  
  if (doplot){
    p <- soft_reg$plot + xlab("log(term)") + ylab("log(fluctuation)")
  } else {
    p <- NULL
  }
  
  return(list(fluctuations=fluctuations, alpha_evid=alpha_evid, 
              alpha_EbetP=alpha_evid_EbetP, alpha_EbetP_short_term=alpha_evid_EbetP_short_term, alpha_EbetP_long_term=alpha_evid_EbetP_long_term, 
              alpha_min_content=alpha_min_content, alpha_max_contour=alpha_max_contour, plot=p))
}
################################
fluctions_cloud <- function(ts, min_win_size=10, max_win_size=length(ts)-2, n_sizes=1000, n_wind=100,
                            transparency=0.2, points_size=0.1, slope=F, loc_slopes=F, log=T, do.plot=F){
  # min_win_size=10; max_win_size=length(ts)-2; n_sizes=1000; n_wind=100; transparency=0.2; points_size=0.1; slope=F; loc_slopes=F; do.plot=F
  fluctuations <- compute_fluctuations(ts, n_sizes=n_sizes, n_wind=n_wind, 
                                       min_win_size=min_win_size, max_win_size=max_win_size)
  if (log){
    log_fluc <- log10(fluctuations)
  } else {
    log_fluc <- fluctuations
  }
  
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
plot_cloud_fluctuations <- function(fluctuations, transparency=0.1, points_size=230*exp(-max(fluctuations$window_size)*0.051), fluc_type="dfa", log=F, title=NULL){
  
  if (log){
    fluctuations <- log10(fluctuations)
    names(fluctuations)[1:2] <- c('log_n', 'log_fluc')
    n <- max(fluctuations$log_n)
    p <- ggplot(fluctuations, aes(x=log_n, y=log_fluc))
  } else{
    n <- max(fluctuations$window_size)
    p <- ggplot(fluctuations, aes(x=window_size, y=fluc))
  }

  p <- p + 
    geom_point(size=points_size, alpha=transparency, color="white") + 
    theme_void() + ggtitle(title) +
    theme(
      plot.background=element_rect(fill="#04304B"),
          #plot.title=element_text(size=40, colour="white", hjust=0.5),
          plot.margin=unit(c(-1, -1, -1, -1), "cm")#, 
          #axis.title=element_text(size=20, colour='white'),
          #axis.text=element_text(colour='white')
      )
  return(p)
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
################################
get_synced_plot <- function(fluc_df, start_idx, stop_idx, ts,
                            size=0.2, alpha=0.2, color="white", bg_color="#04304B", veil_scaling=F, min_fluc=NULL, max_fluc=NULL) {
  
  df_ecg <- data.frame(time = 1:length(ts), val = ts)
  
  p_ecg <- ggplot(df_ecg, aes(x = time, y = val)) +
    geom_line(color = color, linewidth = 0.5) +
    annotate("rect", xmin = start_idx, xmax = stop_idx, 
             ymin = min(ts) - 0.1 * abs(min(ts)), 
             ymax = max(ts) + 0.1 * abs(max(ts)), 
             fill = color, alpha = 0.25) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = bg_color, color = NA),
      panel.background = element_rect(fill = bg_color, color = NA),
      plot.margin      = margin(5, 5, 0, 5)
    )
  
  p_veil <- ggplot(fluc_df, aes(x = window_size, y = fluc)) +
    geom_point(size = size, alpha = alpha, color = color) +
    theme_void() +
    theme(
      plot.background  = element_rect(fill = bg_color, color = NA),
      panel.background = element_rect(fill = bg_color, color = NA),
      plot.margin      = margin(0, 5, 5, 5)
    )
  
  if (veil_scaling){
    p_veil <- p_veil + ylim(c(min_fluc, max_fluc))
  }
  # Combiner avec un fond global sombre sous l'assemblage grid.arrange
  g <- grid.arrange(p_ecg, p_veil, heights = c(1, 2.2))
  
  # Ajout d'un rectangle de fond bleu nuit sous toute la figure
  bg <- rectGrob(gp = gpar(fill = bg_color, col = NA))
  grobTree(bg, g)
}
################################
compute_veils <- function(ts, window_size=300, lag=1, verbose=T, design='sliding'){
  n <- length(ts)
  veil_starts <- seq(1, n - window_size + 1, by = lag)
  
  # Veils computation
  fluctuations_list <- vector("list", length(veil_starts))
  
  veil_start <- 1
  veil_stop  <- veil_start + window_size - 1
  window     <- ts[veil_start:veil_stop]
  
  if (verbose){
    cat("First veil computation...\n")
  }
  fluctuations <- compute_fluctuations(window, windowing="complete", verbose = verbose)
  fluctuations_list[[1]] <- fluctuations
  if (verbose){
    cat("\n")
  }
  
  new_stops <- seq(window_size - lag + 1, window_size)
  
  for (iveil in 2:length(veil_starts)) {
    if (verbose){
      cat('veil', iveil, '-> ')
    }
    
    if (design == 'sliding'){
      v_start <- veil_starts[iveil]
    } else { # 'growing'
      v_start <- 1
    }
    
    v_stop  <- v_start + window_size - 1
    new_ts  <- ts[v_start:v_stop]
    
    fluctuations$start <- fluctuations$start - lag
    fluctuations$stop  <- fluctuations$stop - lag
    fluctuations       <- fluctuations[fluctuations$start >= 1, ]
    
    X <- cumsum(new_ts - mean(new_ts))
    new_rows <- list()
    idx <- 1
    
    cat(length(new_stops), 'new scales to compute: ')
    for (inew_stops in 1:length(new_stops)) {
      if (verbose){
        cat(inew_stops, ' ')
      }
      
      fluc_stop <- new_stops[inew_stops] 
      for (fluc_start in 1:(fluc_stop - 1)) {
        fluc_window <- X[fluc_start:fluc_stop]
        reg_model   <- lm(y ~ x, data.frame(x = 1:length(fluc_window), y = fluc_window))
        fluc        <- sqrt(mean((predict(reg_model) - fluc_window)^2))
        
        new_rows[[idx]] <- data.frame(
          window_size = fluc_stop - fluc_start,
          fluc        = fluc,
          start       = fluc_start,
          stop        = fluc_stop
        )
        idx <- idx + 1
      }
    }
    if (verbose){
      cat("\n")
    }
    
    fluctuations <- rbind(fluctuations, do.call(rbind, new_rows))
    fluctuations_list[[iveil]] <- fluctuations
  }
  return(fluctuations_list)
}
################################
generate_veils_video <- function(fluctuations_list, ts, window_size, lag, save_pdf=F, verbose=T, video_name="fluctuations_veils_video.mp4",
                                 framerate=24, size=0.2, alpha=0.2, color="white", bg_color="#04304B", veil_scaling, design='sliding'){
  n <- length(ts)
  if (design == 'sliding'){
    veil_starts <- seq(1, n - window_size + 1, by = lag)
  } else { # 'growing'
    veil_starts <- rep(1, length(seq(1, n - window_size + 1, by = lag)))
  }
  
  video_output <- paste0("./results/videos/", video_name)
  
  if (veil_scaling){
    # Computation of the most extreme fluctuation values for veils scaling
    min_flucs <- c()
    max_flucs <- c()
    if (veil_scaling){
      for (iveil in 1:length(fluctuations_list)){
        flucs <- fluctuations_list[[iveil]]
        min_1veil <- min(flucs$fluc)
        max_1veil <- max(flucs$fluc)
        min_flucs <- c(min_flucs, min_1veil)
        max_flucs <- c(max_flucs, max_1veil)
      }
      min_fluc <- min(min_flucs)
      max_fluc <- max(max_flucs)
    }
  }
  
  cat("MP4 video generation...\n")
  av_capture_graphics(
    expr = {
      for (iveil in 1:length(veil_starts)) {
        if (verbose){
          cat(iveil, ' ')
        }
        v_start <- veil_starts[iveil]
        v_stop  <- v_start + window_size - 1
        
        # On dessine directement le grobTree complet sans effacement blanc préalable
        grid.draw(get_synced_plot(fluctuations_list[[iveil]], v_start, v_stop, ts, 
                                  size=size, alpha=alpha, color=color, bg_color=bg_color, veil_scaling=veil_scaling, min_fluc=min_fluc, max_fluc=max_fluc))
        if (save_pdf){
          ggsave(file=paste0("./results/videos/veil", iveil, ".png"))
        }
      }
    },
    output    = video_output,
    width     = 900, # 1920
    height    = 500, # 1080
    framerate = framerate,
    res       = 150,
    verbose   = F
  )

  cat("video downloaded as :", video_output, "\n")
  
}
################################
fluc_video <- function(ts, window_size=300, lag=1, save_pdf=F, verbose=T, video_name="fluctuations_veils_video.mp4",
                       framerate=24, size=0.2, alpha=0.2, color="white", bg_color="#04304B", design="sliding", veil_scaling=F){
  
  dir.create("./results/videos", recursive = TRUE, showWarnings = FALSE)
  n <- length(ts)
  veil_starts <- seq(1, n - window_size + 1, by = lag)
  
  if (verbose){
    cat("signal length:", n, "points\n")
    cat("veils dimensions:", window_size, "\n")
    cat("considered lag:", lag, "\n")
    cat("number of veils to compute:", length(veil_starts), "\n")
  }
  
  # Veils computation
  fluctuations_list <- compute_veils(ts, window_size, lag, verbose, design)
  
  # Video rendering
  generate_veils_video(fluctuations_list, ts, window_size, lag, save_pdf=save_pdf, verbose=verbose, 
                       video_name, framerate=framerate, size=size, alpha=alpha, color=color, bg_color=bg_colo, veil_scaling=veil_scaling)
}
################################

