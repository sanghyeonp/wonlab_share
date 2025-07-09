if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, httr)

###################
# Query setup
###################
### Known drugs 수를 파악하는 query.
query.counts <<- "
query known_drug_count($ensemblId: String!){
    target(ensemblId: $ensemblId){
        knownDrugs{
            count
        }
    }
}
"

### EnsemblID에 대한 drug 정보를 뽑아오는 query. 
query.drugs <<- "
query targetAnnotation($ensemblId: String!, $size: Int!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    symbolSynonyms{
			label
    }
    knownDrugs(size: $size){
        count
    	rows {
            drugId
            ctIds
            approvedName
            drug{
              name
              isApproved
              approvedIndications
            }
            drugType
            disease{
              id
              name
            }
            status
            phase
            mechanismOfAction
            targetClass
      }
    }
	}
}
"

###################
# API setup
###################
# Set base URL of GraphQL API endpoint
base_url <<- "https://api.platform.opentargets.org/api/v4/graphql"

###################
# Functions
###################
parse_drug_count <- function(gene_id){
    variables <- list("ensemblId" = gene_id)
    post_body <- list(query=query.counts, variables=variables)
    r <- POST(url=base_url, body=post_body, encode='json')
    content_list <- content(r, as = "parsed", simplifyVector = T)
    n_size <- content_list$data$target$knownDrugs$count
    return(n_size)
}

parse_known_drugs <- function(gene_id, n_count){
    if (n_count == 0){return(data.frame())}
    
    variables <- list("ensemblId" = gene_id, "size" = n_count)
    post_body <- list(query = query.drugs, variables = variables)
    r <- POST(url=base_url, body=post_body, encode='json')
    
    content_list <- content(r, as = "parsed", simplifyVector = T)
    drug_rows <- content_list$data$target$knownDrugs$rows
    
    # Unnest.
    drug_rows$ctIds <- sapply(drug_rows$ctIds, function(x){paste(x, collapse=";")})
    drug_rows$targetClass <- sapply(drug_rows$targetClass, function(x){paste(x, collapse=";")})
    drug_rows <- drug_rows %>%
        tidyr::unnest_wider(col = drug, names_sep = ".") %>%
        tidyr::unnest_wider(col = disease, names_sep = ".")
    drug_rows$drug.approvedIndications <- sapply(drug_rows$drug.approvedIndications, function(x){paste(x[[1]], collapse=";")})
    
    # Add gene info.
    drug_rows$ensemblID <- content_list$data$target$id
    drug_rows$gene_symbol <- content_list$data$target$approvedSymbol
    drug_rows$gene_symbol_synonym <- paste(content_list$data$targe$symbolSynonyms$label, collapse=";")
    
    return(as.data.frame(drug_rows))
}


OpenTargetsPlatform_KnownDrugs <- function(gene_id){
    n_count <- parse_drug_count(gene_id)
    df.parsed <- parse_known_drugs(gene_id, n_count)
    return(df.parsed)
}

###################
# Parse (Example)
###################
# gene_id <- "ENSG00000110245" # APOC3
# df <- OpenTargetsPlatform_KnownDrugs(gene_id)
# 
# gene_id <- "ENSG00000169174" # PCSK9
# df <- OpenTargetsPlatform_KnownDrugs(gene_id)
# 
# gene_id <- "ENSG00000158186" # MRAS
# df <- OpenTargetsPlatform_KnownDrugs(gene_id)
# 
# gene_id <- "ENSG00000001626" # CFTR
# df <- OpenTargetsPlatform_KnownDrugs(gene_id)
