list.of.packages <- c("readr", "dplyr", "tidyr", "rstudioapi")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tidyr)
    library(rstudioapi)
})
