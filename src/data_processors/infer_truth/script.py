import anndata as ad
import liana as li

## VIASH START
par = {
  "input": "resources_test/common/singlecell_broadinstitute_scp2167_human_brain/dataset.h5ad",
  "output": "output.h5ad"
}
meta = {
  "resources_dir": "src/data_processors/infer_truth"
}
## VIASH END

import sys
sys.path.append(meta["resources_dir"])
from funs import non_mirrored_product, onehot_groupby, format_truth

# read the dataset
adata = ad.read_h5ad(par["input"])
adata.X = adata.layers["normalized"]
adata.var.index = adata.var["feature_name"]

# one hot encode cell types
li.ut.spatial_neighbors(adata, bandwidth=1000, max_neighbours=10)

# get organism
organism = adata.uns['dataset_organism']
resource_name_map = {
  "homo_sapiens": "consensus",
  "mus_musculus": "mouseconsensus"
}
resource_name = resource_name_map[organism]

# run LR
lr = li.mt.bivariate(adata,
                     global_name='morans',
                     local_name=None,
                     use_raw=False,
                     resource_name=resource_name,
                     verbose=True,
                     n_perms=1000)

# run CP
ctdata = onehot_groupby(adata, groupby='cell_type')
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
