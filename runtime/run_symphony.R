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
filename = paste('symphony_', 'q=', m, '_r=', n, sep = '')
print(filename)

# read in data
ref_exp = readRDS(paste0('benchmark_datasets/', 'ref_exp_', n, '.rds'))
ref_metadata = readRDS(paste0('benchmark_datasets/', 'ref_metadata_', n, '.rds'))

query_exp = readRDS(paste0('benchmark_datasets/', 'query_exp_', m, '.rds'))
query_metadata = readRDS(paste0('benchmark_datasets/', 'query_metadata_', m, '.rds'))

#------------------------------------------------------------

print('Timing reference building...')
system.time({
reference = buildReference(
    ref_exp,
    ref_metadata,
    vars = c('donor', 'batch'),  # variables to integrate over
    K = 100,                     # number of Harmony clusters
    verbose = TRUE,
    do_umap = FALSE,         # can set to FALSE if want to run umap separately later
    do_normalize = FALSE,    # set to TRUE if input counts are not normalized yet.
    vargenes_method = 'vst', # use vst
    topn = 2000,             # number of variable genes to use
    theta = c(2,2),
    d = 20                  # number of PCs
)
})

print('Timing query mapping...')
system.time({
    query = mapQuery(query_exp, data.frame(query_metadata), reference, 
                 vars = c("donor", "batch"), do_normalize = FALSE, do_umap = FALSE)
})