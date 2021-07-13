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

qexpfile = args[1]
qmetafile = args[2]

# read in data
ref_exp = readRDS(paste0('benchmark_datasets/', 'ref_exp_50000.rds'))
ref_metadata = readRDS(paste0('benchmark_datasets/', 'ref_metadata_50000.rds'))

query_exp = readRDS(qexpfile)
query_metadata = readRDS(qmetafile)

#------------------------------------------------------------

print('Timing reference building...')
system.time({
reference = buildReference(
    ref_exp,
    ref_metadata,
    vars = c('donor'),  # variables to integrate over
    K = 100,              # number of Harmony clusters
    verbose = TRUE,
    do_umap = FALSE,       # can set to FALSE if want to run umap separately later
    do_normalize = FALSE, # set to TRUE if input counts are not normalized yet.
    vargenes_method = 'vst', # use vst
    topn = 2000,             # number of variable genes to use
    theta = c(2),
    d = 20                  # number of PCs
)})

print('Timing query mapping...')
system.time({
query = mapQuery(query_exp, data.frame(query_metadata), reference, 
                 vars = c("donor"), do_normalize = FALSE, do_umap = FALSE)
})