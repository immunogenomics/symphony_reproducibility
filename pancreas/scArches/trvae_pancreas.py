import scarches as sca
exec(open('/data/srlab/ik936/symphony/scArches/python/libs.py').read())
exec(open('/data/srlab/ik936/symphony/scArches/python/trvae_options.py').read())
exec(open('/data/srlab/ik936/symphony/scArches/python/functions.py').read())
outdir = '/data/srlab/ik936/symphony/scArches/results'

outdir = '/data/srlab/ik936/symphony/scArches/results'
condition_key = 'batch_nn'
cell_type_key = 'cell_type'

adata = sc.read_h5ad('/data/srlab/ik936/symphony/scArches/data/pancreas_mapping_cache.h5ad')
adata = remove_sparsity(adata)
time_mapping = do_trvae_scarches(
    adata[adata.obs.batch_nn.isin(['celseq', 'celseq2', 'c1', 'smartseq']), :].copy(),
    adata[adata.obs.batch_nn.isin(['mouse1', 'mouse2', 'human1', 'human2', 'human3', 'human4']), :].copy(),
    'pancreasmapindrop', outdir
)

adata = sc.read_h5ad('/data/srlab/ik936/symphony/scArches/data/pancreas_denovo_cache.h5ad')
adata = remove_sparsity(adata)
time_denovo = do_trvae(adata, True, 'pancreasdenovo', outdir)

## Cache those times 
times = {
    'pancreas_trvae_denovo':time_denovo, 
    'pancreas_trvae_mapping':time_mapping
}
with open('{}/pancreas_time.json'.format(outdir), 'w') as fp:
    json.dump(times, fp)

