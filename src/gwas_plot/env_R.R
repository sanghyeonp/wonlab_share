list.of.packages <- c("data.table", "dplyr", "ggplot2", "ggrepel")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
})

