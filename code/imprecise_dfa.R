
# Variables and screen cleaning
graphics.off(); cat("\014"); rm(list=ls()) ;options(warn=-1)

library(ggplot2)

# Pseudo-manual DFA
setwd("/Users/nsc/Google Drive/Recherche/Projets/Escarres")
df <- read.csv("./data/sample_heartbeat.csv", skip=1)
names(df) <- c('time', 'mV') # columns renaming
df$time <- as.numeric(substr(df$time, 4, nchar(df$time)-1)) # data conversion (string to numerics)
freq <- 256
peaks <- rsleep::detect_rpeaks(df$mV, freq)
peaks <- diff(peaks)
dfa_result <- nonlinearTseries::dfa(rsleep::detect_rpeaks(df$mV, freq), do.plot=F)
dfa_estimate <- nonlinearTseries::estimate(dfa_result)
x=attr(dfa_estimate, "fitted")$x
y=attr(dfa_estimate, "fitted")$y
X=log(x)
Y=log(y)
reg_model <- lm(data.frame(y=Y, x=X[1:length(Y)]))
reg_model$coefficients["x"]
dfa_estimate
summary(reg_model)

# Real manual DFA
x <- rsleep::detect_rpeaks(df$mV, freq)#x=df$mV
N <- length(x)
X <- x - mean(x)
times <- 4
while (tail(times, 1) < N){
  times <- c(times, 2 * tail(times, 1))
}
times <- times[-length(times)]
F_n <- c()
n=times[1]
for (n in times){
  segts_mat <- matrix(X, ncol=n, byrow=T)
  rmses <- c()
  for (i in 1:nrow(segts_mat)){
    seg <- segts_mat[i, ]
    reg_model <- lm(data.frame(y=seg, x=1:n))
    predicted_seg <- predict(reg_model, data.frame(x=1:n))
    rmses <- c(rmses, sqrt(mean((seg-predicted_seg)^2)))
  }
  F_n <- c(F_n, sqrt(mean(rmses^2)))
}
reg_model <- lm(data.frame(y=log(times), x=log(F_n)))
alpha <- reg_model$coefficients[['x']]

alpha
as.numeric(nonlinearTseries::estimate(nonlinearTseries::dfa(x, do.plot=F)))

df <- read.csv('/Users/nsc/Downloads/archive-4/ptbdb_normal.csv', header=F)
