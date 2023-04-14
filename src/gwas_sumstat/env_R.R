list.of.packages <- c("argparse", "readr", "dplyr", "tidyr", "rstudioapi")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

suppressPackageStartupMessages({
    require(argparse)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(rstudioapi)
})
