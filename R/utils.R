fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}

# Colors for PBMCs
pbmc_colors = c("B" = "#66C2A5", 
              "DC" = "#FC8D62",
              "HSC" = "#8DA0CB",
              "MK" = "#E78AC3", 
              "Mono_CD14" = "#A6D854",
              "Mono_CD16" = "#f2ec72",
              "NK" = "#62AAEA", 
              "T_CD4" = "#D1C656",
              "T_CD8" = "#968763")

# Colors for pancreas
celltype.colors = c('alpha'="#ed2bb1",
                    'beta'="#239eb3",
                    'gamma'="#d1bfec",
                    'delta'= "#FF6347",
                    'stellate'="#11e38c",
                    'immune'="#812050",
                    'ductal'="#b2d27a",
                    'endothelial'="#4e2da6",
                    'acinar'="#f6bb86",
                    'schwann'="#115d52",
                    'epsilon'="#a1def0",
                    'mast'="#8fec2f")

# Colors for pancreas query donors (Baron et al., 2016)
querydonor.colors = c('human1' = '#b9dbf0',
                      'human2' = '#77a1ba',
                      'human3' = '#6c7ca8',
                      'human4' = '#364261',
                      'mouse1' = '#e68c8c',
                      'mouse2' = '#b35757')

# Colors for fetal liver hematopoeisis
group.colors = c(   'B cell'='#f2bd80',
                    'DC precursor'='#1d6d1f',
                    'DC1'='#8c3ba0',
                    'DC2'='#6533ed',
                    'Early Erythroid'='#83e3f0',
                    'Early lymphoid/T'='#fd5917',
                    'Endothelial cell'='#4f8c9d',
                    'Fibroblast'='#eb1fcb',
                    'Hepatocyte'='#f5cdaf',
                    'HSC_MPP'='#9698dc',
                    'ILC precursor'='#20f53d',
                    'Kupffer Cell'='#f283e3',
                    'Late Erythroid'='#ffb2be',
                    'Mast cell'='#f3d426',
                    'Megakaryocyte'='#5ebf72',
                    'MEMP'='#a67649',
                    'Mid Erythroid'='#2f5bb1',
                    'Mono-Mac'='#90a479',
                    'Monocyte'='#f6932e',
                    'Monocyte precursor'='#d59e9a',
                    'Neut-myeloid prog.'='#caf243',
                    'NK'='#38b5fc',
                    'pDC precursor'='#c82565',
                    'Pre pro B cell'='#d6061a',
                    'pre-B cell'='#e36f6f',
                    'pro-B cell'='#1dfee1',
                    'VCAM1+ EI macro.'='#506356',
                    'centroid' ='black')

# Custom ordering to match original author publication ordering of states
group.ordering = c("HSC_MPP", "Pre pro B cell", 'pro-B cell', 'pre-B cell', 'B cell',
            'ILC precursor', 'Early lymphoid/T', 'NK', 'Neut-myeloid prog.',
            'pDC precursor','DC precursor', 'DC1', 'DC2', 'Monocyte precursor', 'Monocyte', 
            'Mono-Mac', 'Kupffer Cell', 'VCAM1+ EI macro.', 'MEMP', 'Mast cell',
            'Megakaryocyte', 'Early Erythroid', 'Mid Erythroid', 'Late Erythroid',
            'Endothelial cell', 'Fibroblast', 'Hepatocyte') 


#' Basic function to plot cells, colored and faceted by metadata variables
#' 
#' @param metadata metadata, with UMAP labels in UMAP1 and UMAP2 slots
#' @param title Plot title
#' @param color.by metadata column name for phenotype labels
#' @param facet.by metadata column name for faceting
#' @param color.mapping custom color mapping
#' @param show.legend Show cell type legend

plotBasic = function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                        title = 'Query',         # Plot title
                        color.by = 'cell_type',  # metadata column name for coloring
                        facet.by = NULL,         # (optional) metadata column name for faceting
                        color.mapping = NULL,    # custom color mapping
                        legend.position = 'right') {  # Show cell type legend
    
    p = umap_labels %>%
            dplyr::sample_frac(1L) %>% # permute rows randomly
            ggplot(aes(x = UMAP1, y = UMAP2)) + 
            geom_point_rast(aes(col = get(color.by)), size = 0.3, stroke = 0.2, shape = 16)
        if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
    
    # Default formatting
    p = p + theme_bw() +
            labs(title = title, color = color.by) + 
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.position=legend.position) +
            theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
            guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = FALSE)

    if(!is.null(facet.by)) {
        p = p + facet_wrap(~get(facet.by)) +
                theme(strip.text.x = element_text(size = 12)) }    
    return(p)
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