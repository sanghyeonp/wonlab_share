### Load packages from env_R.R
list.of.packages <- c("argparse", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

require(argparse)
library(data.table)

### Common files
table_file <- "/data/UKbiobank/ukb_uni20220812.tab"

################################################
### Define parser arguments
parser <- argparse::ArgumentParser(description=":: Extract phenotypes in UK Biobank using field ID ::", 
                                formatter_class="argparse.ArgumentDefaultsHelpFormatter")

parser$add_argument("--field", nargs = "*", default="NA",
                    help="Specify the UK Biobank fields to extract. If you want to extract multipel fields, specify as follows: --field 129 130")
parser$add_argument("--retain-all-instances", dest="retain_all_instances", action="store_true", 
                    help="Specify to retain all the instances of specified fields. Default=only the first instance is retained.")
parser$add_argument("--drop-na", dest="drop_na", action="store_true", 
                    help="Specify to drop the rows if any field has NA. Default=FALSE.")
parser$add_argument("--outf", required=FALSE, default="NA",
                    help="Specify the name of the output file. Default = UKB_phenotype_extract.tsv")
parser$add_argument("--delim-out", dest="delim_out", required=FALSE, default="tab",
                    help="Specify the delimiter for the output file. Options = [tab, comma, whitespace]. Default = tab.")
parser$add_argument("--outd", required=FALSE, default="NA",
                    help="Specify the output directory. Default = current working directory.")
parser$add_argument("--n-cores", dest="n_cores", required=FALSE, default=1,
                    help="Specify the number of cores to use. Default = 1.")
################################################
### Get parser arguments
args <- parser$parse_args()
field_list <- args$field

retain_all_instances <- args$retain_all_instances
drop_na <- args$drop_na

outf <- args$outf
outd <- args$outd
if (outf == "NA"){
    outf <- "UKB_phenotype_extract.tsv"
}
if (outd == "NA"){
    outd <- getwd()
}
out_path <- paste(outd, outf, sep="/")

delim_out <- args$delim_out
delim_map <- list(tab = "\t", comma = ",", whitespace = " ")
################################################

search_datafield <- function(data_fields, foi, retain_all_instances=FALSE){
    # data_fields: data fields in the tab file
    # foi: field of interest
    foi <- paste0("f.", foi)
    data_fields_new <- NULL
    for (field in data_fields){
        data_fields_new <- c(data_fields_new, paste(strsplit(field, "[.]")[[1]][1:2], collapse="."))
    }
    
    retained_fields <- data_fields[grep(paste0("^", foi, "$"), data_fields_new)]
    
    if (!retain_all_instances){
        foi <- paste0(foi, ".0.0")
        retained_fields <- retained_fields[grep(paste0("^", foi, "$"), retained_fields)]
    }
    return(retained_fields)
}


subset_ukb_tab <- function(df_tab, foi_list, retain_all_instances=FALSE){
    data_fields <- NULL
    for (foi in foi_list){
        data_fields <- c(data_fields, search_datafield(colnames(df_tab), foi, retain_all_instances = retain_all_instances))
    }
    
    df_tab_sub <- subset(df_tab, select=c("f.eid", data_fields))
    return(df_tab_sub)
    }

################################################
### Read table file
cat("\n::Run:: Reading UKB phenotype table")
cat(paste0("\n\tCurrent phenotype table being used: ", table_file))
tab_data <- fread(table_file, nThread = args$n_cores)

### Convert fields specified as character
field_list <- as.character(field_list)

### Extraction
cat("\n::Run:: Extract specified fields from UKB phenotype table")
cat(paste0("\n\tSpecified fields: ", paste(field_list, collpase=", ")))
tab_data_sub <- subset_ukb_tab(df_tab=tab_data, 
                                foi_list=field_list, 
                                retain_all_instances=retain_all_instances
                                )

if (!retain_all_instances){
    cat("\n\tOnly the first instance of each field is retained.")
} else {
    cat("\n\tAll instances of corresponding field are retained.")
}


### Remove NA
cat("\n::Run:: Handle rows with NA")
tab_data_sub2 <- tab_data_sub[complete.cases(tab_data_sub), ]
cat(paste0("\n\tNumber of rows before: ", nrow(tab_data_sub)))
cat(paste0("\n\tNumber of rows with NA: ", nrow(tab_data_sub) - nrow(tab_data_sub2)))

if (drop_na){
    cat("\n\tDropping rows with NA...")
    tab_data_sub <- tab_data_sub2
    rm(tab_data_sub2)
    cat(paste0("\n\tNumber of rows after: ", nrow(tab_data_sub)))
} else{
    cat("\n\tKeeping all the rows...")
}

### Save the 
cat("\n::Run:: Save output")
cat(paste0("\n\tSaved at: ", out_path))
write.table(tab_data_sub, 
            out_path, 
            sep = delim_map[[delim_out]],
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
