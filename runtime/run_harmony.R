.libPaths(rev(.libPaths()))

suppressPackageStartupMessages({
    library(harmony)
    library(symphony)
    library(data.table)
    library(Matrix)
    library(singlecellmethods)
})

set.seed(0)

# arguments
args = commandArgs(trailingOnly=TRUE)

m = as.numeric(args[1]) # num query cells
n = as.numeric(args[2]) # num reference cells

# name the output
prefix = '/data/srlab2/jkang/referencemapping/tbru/runtime4_figure/'
filename = paste('harmony_', 'q=', m, '_r=', n, sep = '')
print(filename)

# read in data
ref_exp = readRDS(paste0('benchmark_datasets/', 'ref_exp_', n, '.rds'))
ref_metadata = readRDS(paste0('benchmark_datasets/', 'ref_metadata_', n, '.rds'))

query_exp = readRDS(paste0('benchmark_datasets/', 'query_exp_', m, '.rds'))
query_metadata = readRDS(paste0('benchmark_datasets/', 'query_metadata_', m, '.rds'))

#------------------------------------------------------------
total_exp = cbind(ref_exp, query_exp)
total_metadata = rbind(ref_metadata, query_metadata)

print('Timing de novo Harmony...')

system.time({
    
total = symphony::buildReference(
    total_exp,
    total_metadata,
    vars = c('donor', 'batch'),  # variables to integrate over
    K = 100,                     # number of Harmony clusters
    verbose = TRUE,              # metadata column specifying groups for vargene selection
    do_umap = FALSE,             # can set to FALSE if want to run umap separately later
    do_normalize = FALSE,        # set to TRUE if input counts are not normalized yet.
    vargenes_method = 'vst',     # use vst
    topn = 2000,                 # number of variable genes to use per donor
    theta = c(2,2),
    d = 20                       # number of PCs
)

print('Finished.')
})