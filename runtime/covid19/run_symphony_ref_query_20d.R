.libPaths(rev(.libPaths()))

suppressPackageStartupMessages({
    library(harmony)
    library(symphony)
    library(data.table)
    library(Matrix)
    library(singlecellmethods)
})

set.seed(0)

# read in data
ref_exp = readRDS('ref_query_datasets/ref_exp_280samples.rds')
ref_metadata = readRDS('ref_query_datasets/ref_meta_280samples.rds')

query_exp = readRDS('ref_query_datasets/query_exp_14samples.rds')
query_metadata = readRDS('ref_query_datasets/query_meta_14samples.rds')

#------------------------------------------------------------

print('Timing reference building...')
system.time({
reference = buildReference(
    ref_exp,
    ref_metadata,
    vars = c('Sample.name', 'dataset'),  # variables to integrate over
    K = 100,                     # number of Harmony clusters
    verbose = TRUE,
    do_umap = FALSE,         # can set to FALSE if want to run umap separately later
    do_normalize = FALSE,    # set to TRUE if input counts are not normalized yet.
    vargenes_method = 'vst', # use vst
    topn = 1301,             # all var genes original authors used
    theta = c(2.5, 1.5),
    d = 20                  # number of PCs
)
})

print('Running UMAP')
system.time({

# Run UMAP using Seurat default parameters
set.seed(123)
my_umap = uwot::umap(
            t(reference$Z_corr), n_neighbors = 30, learning_rate = 0.5, init = "laplacian", 
            metric = 'cosine', fast_sgd = FALSE, n_sgd_threads = 1, # for reproducibility
            min_dist = 0.2, n_threads = 4, ret_model = TRUE)

# Save UMAP model
save_uwot_path = '/data/srlab2/jkang/symphony_reproducibility/runtime/covid19/results/covid19_reference_280samples_20d.umap'
model = uwot::save_uwot(my_umap, file = save_uwot_path, unload = FALSE, verbose = FALSE)
reference$save_uwot_path = save_uwot_path

# Save UMAP coordinates
colnames(my_umap$embedding) = c('UMAP1', 'UMAP2')
reference$umap = my_umap

})

print('Saving reference')
saveRDS(reference, '/data/srlab2/jkang/symphony_reproducibility/runtime/covid19/results/covid19_reference_280samples_20d.rds')


print('Timing query mapping without UMAP...')
system.time({
    query = mapQuery(query_exp, data.frame(query_metadata), reference, 
                 vars = c('Sample.name'), do_normalize = FALSE, do_umap = FALSE)
})


print('Timing query mapping with UMAP...')
system.time({
    query = mapQuery(query_exp, data.frame(query_metadata), reference, 
                 vars = c('Sample.name'), do_normalize = FALSE, do_umap = TRUE)
})

print('Saving query')
saveRDS(query, '/data/srlab2/jkang/symphony_reproducibility/runtime/covid19/results/covid19_query_14samples_20d.rds')