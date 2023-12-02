# Functions library for Detrended Distributional Analysis (DDA)
library(nonlinearTseries); library(ggplot2)

url <- "https://raw.githubusercontent.com/sutton-charani/uncertain_DFA/main/code/www/tabeau.jpg"
download.file(url, destfile = "tableau.jpg")
tableau <- readJPEG("tableau.jpg")
###############
tab_plot <- function(x, y, title){
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




