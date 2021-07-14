library(harmony)
library(data.table)
library(matrixStats)
library(patchwork)
library(Matrix)
library(uwot)
library(plyr)
library(dplyr)
library(singlecellmethods)
library(Rcpp)
library(irlba)
library(class) # for knn
library(presto)
library(Seurat)
library(RANN) # for nn
library(reshape2)
library(purrr)
library(symphony)

# Fig. 7 and Supp Fig. 15

k <- 50

exprs_norm = readRDS("/path/exprs_norm.rds")
meta_data = readRDS("/path/meta_data.rds")

meta_data <- meta_data[meta_data$memTgate,]

idx_query = sample(names(table(meta_data$id)), 54)
idx_query = which(meta_data$id %in% idx_query)
ref_exp = exprs_norm[, -idx_query]
ref_metadata = meta_data[-idx_query, ]
query_exp = exprs_norm[, idx_query]
query_metadata = meta_data[idx_query, ]

rm(exprs_norm,meta_data);gc()

reference = buildReference(
    ref_exp,
    ref_metadata,
    vars = c('donor','batch'),  # variables to integrate over
    K = 100,            # number of Harmony clusters
    verbose = TRUE,
    do_umap = TRUE,     # can set to FALSE if want to run umap separately later
    do_normalize = FALSE, # set to TRUE if input counts are not normalized yet.
    weightedPCA = FALSE,
    pca_function = 'irlba', # use irlba
    vargenes_method = 'vst', # use vst
    topn = 2000,             # number of variable genes to use
    theta = c(2,2),
    d = 20,                  # number of PCs
    save_uwot_path = paste0('/path/symphony_umap_54_cca_k', k, '.rds')
)

query = mapQuery(query_exp, data.frame(query_metadata), reference, vars = c("donor", "batch"), do_normalize = FALSE)

exprs_norm = readRDS("/path/adt_exprs_norm.rds")
harmony_res <- readRDS("/path/harmony_res_donor_batch_memT.rds")

library(RANN)
all_nn <- nn2(harmony_res, k = k, eps = 0)
all_pred <- apply(all_nn$nn.idx, 1, function(x) {rowMeans(exprs_norm[,x])})

ref_exp = exprs_norm[, -idx_query]
query_exp = exprs_norm[, idx_query]
rm(exprs_norm);gc()

test_nn <- nn2(t(reference$Z_corr), t(query$Z), k = k, eps = 0)
test_pred <- apply(test_nn$nn.idx, 1, function(x) {rowMeans(ref_exp[,x])})

rm(ref_exp, query_exp);gc()

all_pred <- all_pred[row.names(all_pred) != "MouseIgG",]
test_pred <- test_pred[row.names(test_pred) != "MouseIgG",]

test_cor <- sapply(1:30, function(x){cor(all_pred[x,idx_query], test_pred[x,])})
names(test_cor) <- row.names(all_pred)

all_pred <- Matrix(all_pred, sparse = T)
test_pred <- Matrix(test_pred, sparse = T)

# Fig 7b

centroids = reference$Z_corr %*% t(reference$R)
ref_umap_model = uwot::load_uwot(reference$save_uwot_path, verbose = FALSE)
umap_centroids = uwot::umap_transform(t(centroids), ref_umap_model)
umap_centroids= umap_centroids %>% as.data.frame()
colnames(umap_centroids) = c('UMAP1', 'UMAP2')

colnames(reference$umap$embedding) = c('UMAP1', 'UMAP2')
umap_labels = cbind(reference$meta_data, reference$umap$embedding)

p = umap_labels %>%
    sample_frac(1L) %>% # permute rows randomly
    ggplot(aes(x = UMAP1, y = UMAP2)) +
    stat_density_2d(geom='polygon', aes(alpha = ..level..), 
                    contour_var = "ndensity", bins = 15, fill = "grey40") + lims(x = c(-6,9),y = c(-7.5,5)) +
    geom_point(data = umap_centroids, aes(x = UMAP1, y = UMAP2), size = 0.8, shape = 19) +
    theme_bw() +
    labs(title = 'PCA reference (217 samples, 411,004 cells)', color = '', fill = '') + 
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme(legend.position="none") +
    theme(legend.text = element_text(size=11), axis.text=element_text(size=14))
p

# Fig 7c

plot_shuffled_features <- function(ab, umap, exprs, pct = 0.95) {
    max.cutoff = quantile(exprs[ab,], pct)
    min.cutoff = quantile(exprs[ab,], 1-pct)
    tmp <- sapply(X = exprs[ab,], FUN = function(x) {
        return(ifelse(test = x > max.cutoff, yes = max.cutoff,
            no = x))
    })
    tmp <- sapply(X = tmp, FUN = function(x) {
        return(ifelse(test = x < min.cutoff, yes = min.cutoff,
            no = x))
    })
    umap_res_plot <- cbind(umap, tmp)
    return(ggplot(data = as.data.frame(umap_res_plot)[sample(nrow(umap_res_plot)),] , aes(x = V1, y = V2)) +
      geom_point(mapping = aes(color = tmp), shape = ".") +
      scale_color_viridis(option = "plasma", end = .9) +
      theme_classic() +
      theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank()) +
      labs(title = ab))
}

# Supp Fig 15

mrna_cor_bydonor_agg <- c()
for(x in c(5,10,50)) {
	# Run code above
	# ...
	# ...

	mrna_cor_bydonor <- sapply(1:30, function(x){
		sapply(unique(meta_data$donor[idx_query]), function(y,x){cor(all_pred[x,idx_query][meta_data$donor[idx_query] == y], test_pred[x,][meta_data$donor[idx_query] == y])}, x)})

	mrna_cor_bydonor <- data.frame(mrna_cor_bydonor)
	colnames(mrna_cor_bydonor) <- gsub("\\.1", "", names(test_cor))

	mrna_cor_bydonor_agg <- rbind(mrna_cor_bydonor_agg, mrna_cor_bydonor %>% data.frame(check.names = F) %>% tibble::rownames_to_column("id") %>% 
    gather(marker, correlation, -id) %>% group_by(marker) %>% 
    summarise(mean = mean(correlation), stdev = sd(correlation)) %>% mutate(k = k))
}
options(repr.plot.height = 4,repr.plot.width = 15)
ggplot(mrna_cor_bydonor_agg, aes(x = marker, y = mean, fill = factor(k))) + geom_bar(stat = "identity", width = .7, position=position_dodge(.7)) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
        geom_errorbar(aes(ymin=mean-stdev, ymax=mean+stdev), width=.2, position=position_dodge(.7)) + 
        ylab("Mean Pearson correlation") + xlab("Surface protein") + scale_fill_manual(values = c("red3", "darkgoldenrod1","dodgerblue3"))

