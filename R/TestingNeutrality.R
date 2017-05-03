package_dir <- "/Users/marcwilliams/Google Drive/early_events/testingneutrality/"


load(paste(package_dir, "data/testdata.RData", sep = ""))

library(pracma);
library(dplyr);
library(ggplot2)
library(scales)

source(paste(package_dir, "R/plots.R", sep = ""))
source(paste(package_dir, "R/multiplot.R", sep = ""))
source(paste(package_dir, "R/metrics.R", sep = ""))
source(paste(package_dir, "R/methods.R", sep = ""))
source(paste(package_dir, "R/processing.R", sep = ""))
