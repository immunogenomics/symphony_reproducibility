.libPaths(rev(.libPaths()))

## Script runs Seurat de novo integration pipeline
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
filename = paste('seurat_de_novo_', 'q=', m, '_r=', n, sep = '')
print(filename)

# Read in pre-made datasets
ref_exp = readRDS(paste0('benchmark_datasets/', 'ref_exp_', n, '.rds'))
ref_metadata = readRDS(paste0('benchmark_datasets/', 'ref_metadata_', n, '.rds'))

query_exp = readRDS(paste0('benchmark_datasets/', 'query_exp_', m, '.rds'))
query_metadata = readRDS(paste0('benchmark_datasets/', 'query_metadata_', m, '.rds'))

#------------------------------------------------------------
total_exp = cbind(ref_exp, query_exp)
total_metadata = rbind(ref_metadata, query_metadata)
rownames(total_metadata) = total_metadata$cell_id # prevents Seurat "Some cells in meta.data not present in provided counts matrix" error
print(paste('Read in ', (nrow(total_metadata)), ' metadata rows'))
print(paste('Read in ', (ncol(total_exp)), ' exp cols'))

print('Timing de novo Seurat...')
system.time({

denovo <- CreateSeuratObject(counts = total_exp, project = "tbru", meta.data = total_metadata)
    
print('Splitting object...')
tbru.list <- SplitObject(denovo, split.by = "donor")

print('Finding variable features...')
for (i in 1:length(tbru.list)) {
    # Data already normalized
    tbru.list[[i]] = FindVariableFeatures(tbru.list[[i]], selection.method = "vst", nfeatures = 2000)
}

print('Finding integration anchors...')
# Identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.
set.seed(0)
tbru.anchors <- FindIntegrationAnchors(object.list = tbru.list, dims = 1:20)

print('Integrating data...')
# We then pass these anchors to the IntegrateData function, which returns a Seurat object.
# The returned object will contain a new Assay, which holds an integrated (or 'batch-corrected') expression matrix for all cells, enabling them to be jointly analyzed
tbru.integrated <- IntegrateData(anchorset = tbru.anchors, dims = 1:20) # default is 30, but we use 20 to match Symphony analysis
DefaultAssay(tbru.integrated) <- "integrated"

print('PCA...')
## PCA to get integrated embedding
tbru.integrated <- ScaleData(tbru.integrated, verbose = FALSE)
tbru.integrated <- RunPCA(tbru.integrated, npcs = 20, verbose = FALSE)

print('Finished.')
})