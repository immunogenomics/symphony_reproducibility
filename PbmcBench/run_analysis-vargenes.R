# R script that runs benchmarking analysis for cell type identification
# Calculates F1 scores by cell type

suppressPackageStartupMessages({
    library("optparse")
    library(singlecellmethods)
    library(harmony)
    library(symphony)
    library(tidyverse)
    library(data.table)
    library(presto)
    library(pheatmap)
    library(Matrix)
    library(irlba)
    library(matrixStats)
    library(singlecellmethods)
    library(patchwork)
    library(ggthemes)
    library(class) # for knn1
    
    # for cell type F1
    library(caret)
    library(e1071)
})

# Options for bash
option_list = list(
  make_option(c("--analysis"), type = "character", default = NULL, 
              help = "name of analysis", metavar = "character"),
  make_option(c("--rexp"), type = "character", default = NULL, 
              help = "path to reference expression matrix", metavar = "character"),
  make_option(c("--rlab"), type = "character", default = NULL, 
              help = "path to reference cell type labels", metavar = "character"),
  make_option(c("--qexp"), type = "character", default = NULL, 
              help = "path to query expression matrix", metavar = "character"),
  make_option(c("--qlab"), type = "character", default = NULL, 
              help = "path to query cell type labels", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = getwd(), 
              help = "output directory [default = %default]", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
# retrieve individual arguments using opt$argument

# Read in data and label files
ref_exp = t(readRDS(opt$rexp))
ref_labels = readRDS(opt$rlab)
q_exp = t(readRDS(opt$qexp))
q_labels = readRDS(opt$qlab)

colnames(q_labels) = c('cell_type_orig')

# Build reference
reference = buildReference(
    ref_exp, 
    ref_labels,
    vars = NULL, # no Harmony
    K = 100,
    verbose = TRUE,
    do_umap = FALSE,
    do_normalize = TRUE,
    vargenes_method = 'vst',
    topn = 2000
)

# Map query
query = mapQuery(q_exp, q_labels, reference, do_normalize = TRUE, do_umap = FALSE)

# Predict query cell types using knn
set.seed(0)
query = knnPredict(query, reference, reference$meta_data$cell_type, k = 5, confidence = TRUE)

query$meta_data$cell_type_orig = as.character(query$meta_data$cell_type_orig)
res_knn = evaluate(query$meta_data$cell_type_orig, query$meta_data$cell_type_pred_knn)

saveRDS(res_knn, paste(opt$out, "/", opt$analysis, '_res_knn.rds', sep = ''))
saveRDS(reference, paste(opt$out, "/", opt$analysis, '_reference.rds', sep = ''))
saveRDS(query, paste(opt$out, "/", opt$analysis, '_query.rds', sep = ''))