list.of.packages <- c("argparse", "fmsb", "dplyr", "data.table", "hash")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

suppressPackageStartupMessages({
    require(argparse)
    library(fmsb)
    library(dplyr)
    library(data.table)
    library(hash)
    library(here)
})