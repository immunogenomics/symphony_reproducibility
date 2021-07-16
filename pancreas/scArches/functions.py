#!/PHShome/ik936/anaconda3/envs/gpu/bin/python

def do_scanvi(source_adata, do_save, suffix, outdir):
    tic = time.perf_counter()
    sca.dataset.setup_anndata(source_adata, batch_key=condition_key, labels_key=cell_type_key)
    vae = sca.models.SCANVI(
        source_adata,
        "Unknown",
        n_layers=2,
        encode_covariates=True,
        deeply_inject_covariates=False,
        use_layer_norm="both",
        use_batch_norm="none",
    )
    vae.train(
        n_epochs_unsupervised=vae_epochs,
        n_epochs_semisupervised=scanvi_epochs,
        n_epochs_kl_warmup=warmup_epochs,
        unsupervised_trainer_kwargs=dict(early_stopping_kwargs=early_stopping_kwargs),
        semisupervised_trainer_kwargs=dict(metrics_to_monitor=["elbo", "accuracy"],
                                           early_stopping_kwargs=early_stopping_kwargs_scanvi),
        frequency=1
    )
    toc = time.perf_counter()
    if do_save:
        reference_latent = sc.AnnData(vae.get_latent_representation())
        reference_latent.obs = source_adata.obs
        sc.pp.neighbors(reference_latent, n_neighbors=8)
        sc.tl.leiden(reference_latent)
        sc.tl.umap(reference_latent)
        np.savetxt('{}/scanvi_{}_umap.csv'.format(outdir, suffix), reference_latent.obsm['X_umap'], delimiter=',')
        np.savetxt('{}/scanvi_{}_X.csv'.format(outdir, suffix), reference_latent.X, delimiter=',')
        reference_latent.obs.to_csv('{}/scanvi_{}_meta.csv'.format(outdir, suffix))
        return(toc - tic)
    else: 
        return([vae, toc - tic])

    

def do_scanvi_scarches(source_adata, target_adata, suffix, outdir):
    adata_full = source_adata.concatenate(target_adata)
    [vae, time_vae] = do_scanvi(source_adata, False, suffix, outdir)
    ## scArches
    tic = time.perf_counter()
    target_adata.obs['orig_cell_types'] = target_adata.obs[cell_type_key].copy()
    target_adata.obs[cell_type_key] = vae.unlabeled_category_
    model = sca.models.SCANVI.load_query_data(target_adata, vae, freeze_dropout = True)
    model._unlabeled_indices = np.arange(target_adata.n_obs)
    model._labeled_indices = []
    model.train(
        n_epochs_unsupervised=surgery_epochs,
        n_epochs_semisupervised=surgery_epochs,
        n_epochs_kl_warmup=warmup_epochs,
        train_base_model=False,
        semisupervised_trainer_kwargs=dict(metrics_to_monitor=["accuracy", "elbo"], 
                                           weight_decay=0,
                                           early_stopping_kwargs=early_stopping_kwargs_surgery
                                          ),
        frequency=1
    )
    toc = time.perf_counter()
    
    ## Latent data 
    full_latent = sc.AnnData(model.get_latent_representation(adata=adata_full))
    full_latent.obs = adata_full.obs
    sc.pp.neighbors(full_latent)
    sc.tl.leiden(full_latent)
    sc.tl.umap(full_latent)
    np.savetxt('{}/scanvi_{}_umap.csv'.format(outdir, suffix), full_latent.obsm['X_umap'], delimiter=',')
    np.savetxt('{}/scanvi_{}_X.csv'.format(outdir, suffix), full_latent.X, delimiter=',')
    full_latent.obs.to_csv('{}/scanvi_{}_meta.csv'.format(outdir, suffix))
    return([time_vae, toc - tic])
    
def do_trvae(source_adata, do_save, suffix, outdir, do_norm=False, find_vargenes=False):
    tic = time.perf_counter()
    print('Ncells: {}', source_adata.shape[0])
    if do_norm:
        print('Start normalize')
        sc.pp.normalize_total(source_adata, target_sum=1e4)
        sc.pp.log1p(source_adata)
        
    if find_vargenes:
        print('Start vargenes')
        sc.pp.highly_variable_genes(source_adata, n_top_genes=2000)
        vargenes=np.where(source_adata.var.highly_variable)[0]
        source_adata = source_adata[:, vargenes].copy()
    
    print('Initialize model')
    vae = sca.models.TRVAE(
        adata=source_adata,
        condition_key=condition_key,
        recon_loss='mse', ## NECESSARY FOR NON-COUNT DATA 
        conditions=source_adata.obs[condition_key].unique().tolist(),
        hidden_layer_sizes=[128, 128],
    )
    print('Train model for {} epochs'.format(trvae_epochs))
    vae.train(
        n_epochs=trvae_epochs,
        alpha_epoch_annea=nepochs_alpha,
        early_stopping_kwargs=early_stopping_kwargs
    )
    toc = time.perf_counter()
    if do_save:
        reference_latent = sc.AnnData(vae.get_latent())
        reference_latent.obs = source_adata.obs
#         sc.pp.neighbors(reference_latent, n_neighbors=8)
#         sc.tl.leiden(reference_latent)
#         sc.tl.umap(reference_latent)
#         np.savetxt('{}/trvae_{}_umap.csv'.format(outdir, suffix), reference_latent.obsm['X_umap'], delimiter=',')
        np.savetxt('{}/trvae_{}_X.csv'.format(outdir, suffix), reference_latent.X, delimiter=',')
        reference_latent.obs.to_csv('{}/trvae_{}_meta.csv'.format(outdir, suffix))
        return(toc - tic)
    else: 
        return([vae, toc - tic])
    
def do_trvae_scarches(source_adata, target_adata, suffix, outdir, do_norm=False, find_vargenes=False):
    if do_norm:
        print('Start normalize')
        sc.pp.normalize_total(source_adata, target_sum=1e4)
        sc.pp.log1p(source_adata)
        sc.pp.normalize_total(target_adata, target_sum=1e4)
        sc.pp.log1p(target_adata)
        
    if find_vargenes:
        print('Start vargenes')
        sc.pp.highly_variable_genes(source_adata, n_top_genes=2000)
        vargenes=np.where(source_adata.var.highly_variable)[0]
        source_adata = source_adata[:, vargenes].copy()
        target_adata = target_adata[:, vargenes].copy()
    
    [vae, time_vae] = do_trvae(source_adata, False, suffix, outdir, False, False)
    tic = time.perf_counter()
    print('Initialize query model')
    model = sca.models.TRVAE.load_query_data(adata=target_adata, reference_model=vae)
    print('Initialize query model mapping')
    model.train(
        n_epochs=surgery_epochs,
        alpha_epoch_anneal=nepochs_alpha,
        early_stopping_kwargs=early_stopping_kwargs,
        weight_decay=0
    )
    toc = time.perf_counter()
    ## Latent data 
    adata_full = source_adata.concatenate(target_adata)
    full_latent = sc.AnnData(model.get_latent(adata_full.X, adata_full.obs[condition_key]))
    full_latent.obs = adata_full.obs
#     sc.pp.neighbors(full_latent)
#     sc.tl.leiden(full_latent)
#     sc.tl.umap(full_latent)
#     np.savetxt('{}/trvae_{}_umap.csv'.format(outdir, suffix), full_latent.obsm['X_umap'], delimiter=',')
    np.savetxt('{}/trvae_{}_X.csv'.format(outdir, suffix), full_latent.X, delimiter=',')
    full_latent.obs.to_csv('{}/trvae_{}_meta.csv'.format(outdir, suffix))    
    return([time_vae, toc - tic])
