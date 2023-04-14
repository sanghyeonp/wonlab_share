### Load packages from env_R.R
list.of.packages <- c("tidyverse", "rstudioapi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(tidyverse)
getCurrentFileLocation <-  function()
{
    this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
    if (length(this_file)==0){
        this_file <- rstudioapi::getSourceEditorContext()$path
    }
    return(dirname(this_file))
}

env_file <- paste(getCurrentFileLocation(), "env_R.R", sep="/")
source(env_file)


### Define parser arguments
parser <- argparse::ArgumentParser(description=":: Reformat GWAS summary statistics ::", 
                                formatter_class="argparse.ArgumentDefaultsHelpFormatter")

parser$add_argument("--gwas", required=TRUE,
                    help="Path to the GWAS summary statistics.")
parser$add_argument("--col-map", dest="col_map", nargs = "*", default = "NA=NA",
                    help = "Specify the dictionaries of common column names to the columns in GWAS summary statistics. Options for common column name: [CHR, POS, ALT, REF, EAF, MAF, BETA, SE, PVAL, N]")
parser$add_argument("--keep-unspecified-col", dest="keep_unspecified_col", action="store_true", 
                    help="Specify to retain unspecified columns in --col_map. Default=FALSE.")

### Get parser arguments
args <- parser$parse_args()
gwas <- args$gwas
col_map <- args$col_map
keep_unspecified_col <- args$keep_unspecified_col

################################################

### Read GWAS summary statistics
df <- readr::read_table(gwas, col_names = TRUE)

### Parse column mapping from `--col-map` parser argument
old_col_names1 <- NULL
new_col_names1 <- NULL
for (d in col_map) {
    dict <- strsplit(d, "=")
    old_col_names1 <- c(old_col_names1, dict[[1]][1])
    new_col_names1 <- c(new_col_names1, dict[[1]][2])
}

### Make another column name vector -> used for re-naming the input GWAS summary statistics
old_col_names2 <- old_col_names1
new_col_names2 <- new_col_names1
for (col in colnames(df)) {
    # Check if elem is not in vec1
    if (!(col %in% old_col_names2)) {
        # Add elem to vec1
        old_col_names2 <- c(old_col_names2, col)
        new_col_names2 <- c(new_col_names2, col)
    }
}

### Make column name renaming vector
col_map_dict <- setNames(new_col_names2, old_col_names2)

### Rename the columns
names(df) <- col_map_dict

### Subset only the specified columns in `--col-map` if `--keep-unspecified-col` is not specified.
### If `--keep-unspecified-col` is specified, do not subset.
if (keep_unspecified_col == FALSE){
    df <- subset(df, select=new_col_names1)
}

head(df)