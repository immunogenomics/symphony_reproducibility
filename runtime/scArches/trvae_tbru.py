import scarches as sca
exec(open('/data/srlab/ik936/symphony/scArches/python/libs.py').read())
exec(open('/data/srlab/ik936/symphony/scArches/python/trvae_options.py').read())
exec(open('/data/srlab/ik936/symphony/scArches/python/functions.py').read())
outdir = '/data/srlab/ik936/symphony/scArches/results'

condition_key = 'donor'
datadir = '/data/srlab/ik936/symphony/data/TBRU'

fd = open('/data/srlab/ik936/symphony/data/TBRU/experiment_order.txt', 'r')
experiments = [i.strip() for i in fd.readlines()]

for exp_name in experiments:
    adata_fname = '{}/adata_{}.h5ad'.format(datadir, exp_name)
    print(adata_fname)
    adata = sc.read_h5ad(adata_fname)
    adata = remove_sparsity(adata)
    time_mapping = do_trvae_scarches(
        adata[adata.obs.TYPE=='Ref', ].copy(),
        adata[adata.obs.TYPE=='Query', ].copy(),
        'TBRU_mapping_{}'.format(exp_name), outdir,
        do_norm=True, find_vargenes=True 
    )
    
    
    time_denovo = do_trvae(
        adata, True, 
        'TBRU_denovo_{}'.format(exp_name), outdir,
        do_norm=True, find_vargenes=True
    )

    ## Cache those times 
    times = {}
    times['TBRU_mapping_{}'.format(exp_name)] = time_mapping
    times['TBRU_denovo_{}'.format(exp_name)] = time_denovo

    with open('{}/TBRU_{}_time.json'.format(outdir, exp_name), 'w') as fp:
        json.dump(times, fp)


        
        
        
        
        
