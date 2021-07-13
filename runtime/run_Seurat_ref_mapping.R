.libPaths(rev(.libPaths()))

## Script runs Seurat reference mapping pipeline
## May 25, 2021
suppressPackageStartupMessages({
    library(Seurat)
})

set.seed(0)

# Arguments
args = commandArgs(trailingOnly=TRUE)

m = as.numeric(args[1]) # num query cells
n = as.numeric(args[2]) # num reference cells

# Name the output
prefix = '/data/srlab2/jkang/referencemapping/tbru/runtime4_figure/'
filename = paste('seurat_ref_mapping_', 'q=', m, '_r=', n, sep = '')
print(filename)

# Read in pre-made datasets
ref_exp = readRDS(paste0('benchmark_datasets/', 'ref_exp_', n, '.rds'))
ref_metadata = readRDS(paste0('benchmark_datasets/', 'ref_metadata_', n, '.rds'))

query_exp = readRDS(paste0('benchmark_datasets/', 'query_exp_', m, '.rds'))
query_metadata = readRDS(paste0('benchmark_datasets/', 'query_metadata_', m, '.rds'))

#------------------------------------------------------------

rownames(ref_metadata) = ref_metadata$cell_id # prevents Seurat "Some cells in meta.data not present in provided counts matrix" error
rownames(query_metadata) = query_metadata$cell_id

print('Timing reference building Seurat...')
system.time({

reference <- CreateSeuratObject(counts = ref_exp, project = "tbru_ref", meta.data = ref_metadata)
    
print('Splitting reference object...')
ref.list <- SplitObject(reference, split.by = "donor")

print('Finding reference variable features...')
for (i in 1:length(ref.list)) {
    # Data already normalized
    ref.list[[i]] = FindVariableFeatures(ref.list[[i]], selection.method = "vst", nfeatures = 2000)
}

print('Finding reference integration anchors...')
# Identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.
set.seed(0)
ref.anchors <- FindIntegrationAnchors(object.list = ref.list, dims = 1:20)

print('Integrating reference data...')
# We then pass these anchors to the IntegrateData function, which returns a Seurat object.
# The returned object will contain a new Assay, which holds an integrated (or 'batch-corrected') expression matrix for all cells, enabling them to be jointly analyzed
ref.integrated <- IntegrateData(anchorset = ref.anchors, dims = 1:20) # default is 30, but use 20 to compare with Symphony
DefaultAssay(ref.integrated) <- "integrated"

print('PCA...')
## PCA to get integrated embedding
ref.integrated <- ScaleData(ref.integrated, verbose = FALSE)
ref.integrated <- RunPCA(ref.integrated, npcs = 20, verbose = FALSE)

print('Finished reference building')
})

#------------------------------------------------------------

print('Timing query mapping Seurat...')
system.time({
    
query <- CreateSeuratObject(counts = query_exp, project = "tbru_query", meta.data = query_metadata)

print('Splitting query object...')
query.list <- SplitObject(query, split.by = "donor")

print('Finding transfer anchors...')
anchors <- list()
for (i in 1:length(query.list)) {
    anchors[[i]] <- FindTransferAnchors(
        reference = ref.integrated,
        query = query.list[[i]],
        k.filter = NA,
        reduction = "pcaproject",       # either pcaproject or cca: only two options
        reference.reduction = "pca",    # spca is for ADT and mRNA
        #     reference.neighbors = "spca.annoy.neighbors", 
        dims = 1:20                     # npc is smaller than defaut 50, so set it to 20
    )
}

print('Mapping query...')
# Individually map each of the datasets.
for (i in 1:length(query.list)) {
    query.list[[i]] <- MapQuery(
        anchorset = anchors[[i]], 
        query = query.list[[i]],
        reference = ref.integrated, 
        refdata = list(
            cell_type = "dummy" # TODO: replace with meta variable from Aparna
            ),
        #   reduction = "umap" # 
        reference.reduction = "pca" # Tried to set this to "pca" or "umap", but just give me pca embeddings
        #   reduction.model = "umap"  # setting "umap" caused this step run forever ...
      )
}
print('Finished query mapping.')

})