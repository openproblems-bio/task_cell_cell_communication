import anndata as ad
import liana as li

## VIASH START
par = {
  "input": "resources/datasets/SCP2167/raw_dataset.h5ad",
  "output": "resources/datasets/SCP2167/dataset_with_inferred_truth.h5ad"
}
meta = {
  "resources_dir": "src/dataset_loaders/infer_truth/"
}
## VIASH END

import sys
sys.path.append(meta["resources_dir"])
from funs import non_mirrored_product, onehot_groupby, format_truth

# read the dataset
adata = ad.read_h5ad(par["input"])
adata.uns['dataset_organism'] = 'homo_sapiens' # NOTE:: Delete
# Get needed params
groupby = 'cell_type' 
organism = adata.uns['dataset_organism']

# one hot encode cell types
li.ut.spatial_neighbors(adata, bandwidth=1000, max_neighbours=10)
ctdata = onehot_groupby(adata, groupby=groupby)

organism = adata.uns['dataset_organism']
resource_name = 'mouseconsesnsus' if organism == 'mouse' else 'consensus'

lr = li.mt.bivariate(adata,
                     global_name='morans',
                     local_name=None,
                     use_raw=False,
                     resource_name=resource_name,
                     verbose=True,
                     n_perms=1000)

# Infer Co-localized Cell types
interactions = non_mirrored_product(ctdata.var.index, ctdata.var.index)
ct = li.mt.bivariate(ctdata,
                     global_name='morans',
                     local_name=None,
                     use_raw=False,
                     interactions=interactions,
                     verbose=True,
                     n_perms=1000,
                     x_name='source',
                     y_name='target')


# Format predictions as x, y, boolean
lr = format_truth(lr, x_name='ligand', y_name='receptor')
ct = format_truth(ct, x_name='source', y_name='target')

# cross join
assumed_truth = lr.assign(key=1).merge(ct.assign(key=1), on='key')
# simplify the dataframe
assumed_truth['colocalized'] = assumed_truth['truth_y'] * assumed_truth['truth_x']
assumed_truth = assumed_truth[assumed_truth['colocalized']]
assumed_truth = assumed_truth[['source', 'target', 'ligand', 'receptor', 'colocalized']]

adata.uns["assumed_truth"] = assumed_truth

# save the new dataset
adata.write_h5ad(par["output"], compression="gzip")