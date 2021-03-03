fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}

## Functions to map orthologs between human and mouse.

# Function to obtain mapping of orthologs between two gene lists using biomaRt.
# Returns orthologs table with column names HGNC.symbol, MGI.symbol (or vice versa)
#
# @param genes_domain: vector of gene names in the "from" species 
# @param genes_target: vector of gene names in the "to" species
# @param from either 'mouse' or 'human'
# @param to either 'mouse' or 'human'
get_ortho <- function(genes_domain, genes_target, from, to) {
    tryCatch({
        library(biomaRt)
    }, error = function(e) {
        stop('Must either provide orthologs_table or install biomaRt')
    })
    attr_domain <- ifelse(from == 'human', 'hgnc_symbol', 'mgi_symbol')
    attr_target <- ifelse(to == 'human', 'hgnc_symbol', 'mgi_symbol')
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    mart_domain <- list(human, mouse)[[as.integer(factor(from, c('human', 'mouse')))]]
    mart_target <- list(human, mouse)[[as.integer(factor(to, c('human', 'mouse')))]]
    orthologs <- getLDS(
        attributesL = attr_domain, 
        filtersL = attr_domain, 
        valuesL = genes_domain, 
        martL = mart_domain, 
        attributes = attr_target, 
        filters = attr_target,
        values = genes_target,
        mart = mart_target, 
        uniqueRows = FALSE
    )
    return(orthologs)
}

# Function to turn a ortholog table into a mapping matrix.
# Convert ortholog pairs into sparse linear map 
#
# @param orthologs_dict a dictionary of ortholog pairs (as output by get_ortho())
dict_to_map <- function(orthologs_dict) {
    colnames(orthologs_dict) <- c('to', 'from') 
    X <- data.table(orthologs_dict)[
        , w := 1 / .N, by = from
    ][]
    i <- factor(X$to)
    j <- factor(X$from)
    map_mat <- sparseMatrix(
        i = as.integer(i), 
        j = as.integer(j), 
        x = X$w,
        dimnames = list(to = levels(i), from = levels(j))
    )  
    return(map_mat)
}

# Function that converts a gene expression matrix from one species (human or mouse) to its orthologous
# gene expression in the other species
# 
# @param exprs_matrix input expression matrix
# @genes_target gene names for target species
# @orthologs_table if NULL, pulls from biomaRt; else, orthologs_table should have columns MGI.symbol and HGNC.symbol
map_species <- function(exprs_matrix, genes_target, from='mouse', to='human', orthologs_table=NULL, round_fxn=identity) {
    genes_domain <- rownames(exprs_matrix)
    if (is.null(orthologs_table)) {
        message('No orthologues DF provided, pulling data from Biomart')
        orthologs_table <- get_ortho(genes_domain, genes_target, from, to)
    } else {
        message('CAUTION: orthologs_table should have columns MGI.symbol and HGNC.symbol')
    }
    ## convert ortholog pairs into sparse linear map 
    linear_map <- dict_to_map(orthologs_table)
    ## apply linear map 
    genes_use <- intersect(genes_domain, colnames(linear_map))
    exprs_mapped <- linear_map[, genes_use] %*% exprs_matrix[genes_use, ]
    exprs_mapped <- round_fxn(exprs_mapped)
    return(exprs_mapped)
}