import sys
import anndata as ad
import openproblems as op

## VIASH START
par = {
    'input': 'resources_test/common/cxg_mouse_pancreas_atlas/dataset.h5ad',
    'output_dataset': 'train.h5ad',
    'output_solution': 'test.h5ad'
}
meta = {
    'resources_dir': 'target/executable/data_processors/process_dataset',
    'config': 'target/executable/data_processors/process_dataset/.config.vsh.yaml'
}
## VIASH END

# import helper functions
sys.path.append(meta['resources_dir'])
from subset_h5ad_by_format import subset_h5ad_by_format

# read viash config
config = op.project.read_viash_config(meta["config"])

# read the dataset
input = ad.read_h5ad(par['input'])

print(">> Creating test data", flush=True)
output_dataset = subset_h5ad_by_format(
    input,
    config,
    "output_test"
)

print(">> Creating solution data", flush=True)
output_solution = subset_h5ad_by_format(
    input,
    config,
    "output_solution"
)

print(">> Writing data", flush=True)
output_dataset.write_h5ad(par["output_dataset"], compression="gzip")
output_solution.write_h5ad(par["output_solution"], compression="gzip")
