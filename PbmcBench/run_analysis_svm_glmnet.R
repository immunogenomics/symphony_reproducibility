# R script that runs benchmarking analysis for cell type identification
# Calculates F1 scores by cell type

suppressPackageStartupMessages({
    library("optparse")
    library(symphony)
    library(singlecellmethods)
    library(harmony)
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
  make_option(c("-o", "--out"), type = "character", default = getwd(), 
              help = "output directory [default = %default]", metavar = "character")
);

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
# retrieve individual arguments using opt$argument

a = opt$analysis

# VARGENES------------------------------------------------------------------------
reference = readRDS(paste('results/results_vargenes/', a, '_reference.rds', sep = ''))
query = readRDS(paste('results/results_vargenes/', a, '_query.rds', sep = ''))

# Fitting via glmnet
suppressWarnings(mod.glmnet <- glmnet::cv.glmnet(
    x = t(reference$Z_corr), y = as.factor(reference$meta_data$cell_type),
    family = "multinomial",
    alpha = 0,
    type.multinomial = "grouped"
))

predictions = mod.glmnet %>% predict(t(query$Z), type = 'class', s = 'lambda.min')
res = evaluate(predictions, query$meta_data$cell_type_orig)
saveRDS(res, paste(opt$out, "/", opt$analysis, '_res_vargenes_glmnet.rds', sep = ''))

# Fitting via SVM
model = svm(t(reference$Z_corr), as.factor(reference$meta_data$cell_type)) # radial kernel
predictions = predict(model, t(query$Z))
res = evaluate(predictions, query$meta_data$cell_type_orig)
saveRDS(res, paste(opt$out, "/", opt$analysis, '_res_vargenes_svm.rds', sep = ''))

# 20 DEGs------------------------------------------------------------------------
reference = readRDS(paste('results/results_marker-genes-20/', a, '_reference.rds', sep = ''))
query = readRDS(paste('results/results_marker-genes-20/', a, '_query.rds', sep = ''))

# fitting via glmnet
suppressWarnings(mod.glmnet <- glmnet::cv.glmnet(
    x = t(reference$Z_corr), y = as.factor(reference$meta_data$cell_type),
    family = "multinomial",
    alpha = 0,
    type.multinomial = "grouped"
))

predictions = mod.glmnet %>% predict(t(query$Z), type = 'class', s = 'lambda.min')
res = evaluate(predictions, query$meta_data$cell_type_orig)
saveRDS(res, paste(opt$out, "/", opt$analysis, '_res_degs_glmnet.rds', sep = ''))

# fitting via SVM
model = svm(t(reference$Z_corr), as.factor(reference$meta_data$cell_type))
predictions = predict(model, t(query$Z))
res = evaluate(predictions, query$meta_data$cell_type_orig)
saveRDS(res, paste(opt$out, "/", opt$analysis, '_res_degs_svm.rds', sep = ''))

print('all done!')