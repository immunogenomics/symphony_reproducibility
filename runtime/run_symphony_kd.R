.libPaths(rev(.libPaths()))

suppressPackageStartupMessages({
    library(harmony)
    library(symphony)
    library(data.table)
    library(Matrix)
    library(singlecellmethods)
})

# arguments
args = commandArgs(trailingOnly=TRUE)

k = as.numeric(args[1]) # num centroids
d = as.numeric(args[2]) # num dimensions

# name the output
prefix = '/data/srlab2/jkang/symphony_reproducibility/results_kd/'
filename = paste('symphony_', 'k=', k, '_d=', d, sep = '')
print(filename)

# read in data
ref_exp = readRDS(paste0('benchmark_datasets/', 'ref_exp_50000.rds'))
ref_metadata = readRDS(paste0('benchmark_datasets/', 'ref_metadata_50000.rds'))

query_exp = readRDS(paste0('benchmark_datasets/', 'query_exp_10000.rds'))
query_metadata = readRDS(paste0('benchmark_datasets/', 'query_metadata_10000.rds'))

#------------------------------------------------------------

print('Timing reference building...')

system.time({
reference = buildReference(
    ref_exp,
    ref_metadata,
    vars = c('donor'),  # variables to integrate over
    K = k,              # number of Harmony clusters
    verbose = TRUE,
    do_umap = FALSE,       # can set to FALSE if want to run umap separately later
    do_normalize = FALSE, # set to TRUE if input counts are not normalized yet.
    vargenes_method = 'vst', # use vst
    topn = 2000,             # number of variable genes to use
    theta = c(2),
    d = d                  # number of PCs
)})

print('Timing query mapping...')

system.time({
query = mapQuery(query_exp, data.frame(query_metadata), reference, 
                 vars = c("donor"), do_normalize = FALSE, do_umap = FALSE)
})