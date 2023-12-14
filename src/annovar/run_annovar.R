
run_annovar <- function(
    annov_in,
    genome_build, dbsnp_build,
    nthread
){
    dir_annovar <- "/data1/sanghyeon/wonlab_contribute/combined/software/annovar/annovar"
    dir_humandb <- paste0(dir_annovar, "/humandb/")

    cmd_annovar <- paste0(
        dir_annovar, "/annotate_variation.pl",
        " -filter",
        " -build hg", genome_build,
        " -dbtype ", "avsnp", dbsnp_build, " ", annov_in, " ", dir_humandb,
        " -out ", annov_in,
        " -thread ", nthread 
    )

    system(cmd_annovar, wait = TRUE)

    return(paste0(annov_in, ".hg", genome_build, "_avsnp", dbsnp_build, "_dropped"))
}