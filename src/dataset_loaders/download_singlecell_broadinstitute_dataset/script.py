import scanpy as sc
import pandas as pd
import re
import os
import urllib.request

## VIASH START
par = {
  "dataset_curl_config": "resources/datasets/raw/singlecell_broadinstitute_configs/SCP2169.txt",
  "dataset_id": "singlecell_broadinstitute_scp2169",
  "dataset_name": "Human tonsil Slide-tags snRNA-seq",
  "dataset_url": "https://singlecell.broadinstitute.org/single_cell/study/SCP2169/slide-tags-snrna-seq-on-human-tonsil#study-summary",
  "dataset_reference": "doi: 10.1038/s41586-023-06837-4",
  "dataset_summary": "Slide-tags snRNA-seq data on the human tonsil.",
  "dataset_description": "Recent technological innovations have enabled the high-throughput quantification of gene expression and epigenetic regulation within individual cells, transforming our understanding of how complex tissues are constructed. Missing from these measurements, however, is the ability to routinely and easily spatially localise these profiled cells. We developed a strategy, Slide-tags, in which single nuclei within an intact tissue section are ‘tagged’ with spatial barcode oligonucleotides derived from DNA-barcoded beads with known positions. These tagged nuclei can then be used as input into a wide variety of single-nucleus profiling assays. We used Slide-tags to profile two different stages of development in the mouse brain.\n\nOverall design 	Slide-tags was used to spatially barcode nuclei from 20-micron thick fresh frozen tissue sections. These spatially barcoded nuclei were then used as input for 10x Genomics Chromium v3.1 snRNA-seq.",
  "dataset_organism": "homo_sapiens",
  "output": "resources/datasets/SCP2169/raw_dataset.h5ad"
}
meta = {
  "temp_dir": "/tmp"
}
## VIASH END

temp_dir = f'{meta["temp_dir"]}/downloader_singlecell_broadinstitute_dataset'

if not os.path.exists(temp_dir):
  os.makedirs(temp_dir)

def read_typed_csv(path: str) -> pd.DataFrame:
  col_names = pd.read_csv(path, nrows=1).columns.tolist()
  col_types = pd.read_csv(path, skiprows=0, nrows=1).iloc[0].tolist()

  # Convert type names to appropriate Pandas dtypes
  dtype_mapping = {
      "group": "category",
      "category": "category",
      "TYPE": "string",
      "numeric": "float64",
  }
  col_types = [dtype_mapping.get(t, t) for t in col_types]  

  # Read the rest of the CSV, applying the column names and types
  df = pd.read_csv(
    path,
    skiprows=2,
    names=col_names,
    dtype=dict(zip(col_names, col_types)),
    index_col="NAME"
  )

  return df

with open(par["dataset_curl_config"], "r") as file:
  curl_config = file.readlines()

# split config by \n\n
curl_config_split = "".join(curl_config).split("\n\n")

# turn into a list of dicts with 'key="value"\nkey="value"'
curl_config_dict_entry = '([^=]*)="(.*)"'
curl_config_dicts = [
  dict([
    re.match(curl_config_dict_entry, x).groups()
    for x in config.split("\n")
    if re.match(curl_config_dict_entry, x) is not None
  ])
  for config in curl_config_split
]
# remove empty entries
curl_config_dicts = [x for x in curl_config_dicts if x]

def get_download_info(filename_query):
  return next(
    (
      {
        "url": x["url"],
        "output": x["output"],
        "dest": f"{temp_dir}/{x['output']}"
      }
      for x in curl_config_dicts
      if re.match(filename_query, x["output"])
    ),
    None
  )

# fetch 10x_mtx data. get url for row with "matrix.mtx.gz" in output
matrix_mtx = get_download_info(".*matrix.mtx.gz$")
barcodes_tsv = get_download_info(".*barcodes.tsv.gz$")
features_tsv = get_download_info(".*features.tsv.gz$")

def download_file(info):
  os.makedirs(os.path.dirname(info["dest"]), exist_ok=True)
  urllib.request.urlretrieve(info["url"], info["dest"])

# download mtx
download_file(matrix_mtx)
download_file(barcodes_tsv)
download_file(features_tsv)

adata = sc.read_10x_mtx(os.path.dirname(matrix_mtx["dest"]))
adata
# AnnData object with n_obs × n_vars = 14165 × 36601
#     var: 'gene_ids', 'feature_types'

adata.obs_names = adata.obs_names.str.replace("-1", "")

cluster_info = get_download_info(".*_cluster.csv$")
download_file(cluster_info)
cluster = read_typed_csv(cluster_info["dest"])
cluster.index = cluster.index.str.replace("-1", "")
cluster
#                             X          Y    cell_type
# NAME                                                 
# CTACATTCAGCTTTGA-1 -16.746614  -1.999025   Excitatory
# AACCTTTCACTGGATT-1 -19.254899 -17.433621   Excitatory
# CCATCACGTTAGTCGT-1 -16.703370  -1.978855   Excitatory
# CTCCTTTCAGACCATT-1 -16.329322  -1.227991   Excitatory
# GCAACCGCACCAAATC-1 -20.231160  -2.528652   Excitatory
# ...                       ...        ...          ...
# AGTTCGACAGACGCTC-1   4.417137 -22.180945  Endothelial
# CCTCATGTCCTTATCA-1   4.948904 -23.954339  Endothelial
# TACGCTCTCCATCTCG-1   4.965150 -24.211097  Endothelial
# ATTCCTAGTGGAAATT-1  16.381543   9.752836        Oligo
# GCATCGGAGCACCGTC-1  12.886445  11.456965        Oligo

# [4067 rows x 3 columns]

spatial_info = get_download_info(".*_spatial.csv$")
download_file(spatial_info)
spatial = read_typed_csv(spatial_info["dest"])
spatial.index = spatial.index.str.replace("-1", "")
spatial
#                               X            Y    cell_type
# NAME                                                     
# CTACATTCAGCTTTGA-1  4982.182919  2715.304054   Excitatory
# AACCTTTCACTGGATT-1  3877.271881  1197.231379   Excitatory
# CCATCACGTTAGTCGT-1  2977.881263  4087.932046   Excitatory
# CTCCTTTCAGACCATT-1  4738.346031  2468.137565   Excitatory
# GCAACCGCACCAAATC-1  3010.716289  2978.643333   Excitatory
# ...                         ...          ...          ...
# AGTTCGACAGACGCTC-1  1274.691325  1424.097687  Endothelial
# CCTCATGTCCTTATCA-1  3541.376000  3306.541636  Endothelial
# TACGCTCTCCATCTCG-1  3164.842000  5547.069250  Endothelial
# ATTCCTAGTGGAAATT-1  2257.890000  1391.706675        Oligo
# GCATCGGAGCACCGTC-1  2657.236500  2220.441000        Oligo

# [4067 rows x 3 columns]

metadata_info = get_download_info(".*_metadata.csv$")
download_file(metadata_info)
metadata = read_typed_csv(metadata_info["dest"])
metadata.index = metadata.index.str.replace("-1", "")
metadata
#                    biosample_id       donor_id         species  ... library_preparation_protocol__ontology_label     sex      cluster
# NAME                                                            ...                                                                  
# CTACATTCAGCTTTGA-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female   Excitatory
# AACCTTTCACTGGATT-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female   Excitatory
# CCATCACGTTAGTCGT-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female   Excitatory
# CTCCTTTCAGACCATT-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female   Excitatory
# GCAACCGCACCAAATC-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female   Excitatory
# ...                         ...            ...             ...  ...                                          ...     ...          ...
# AGTTCGACAGACGCTC-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female  Endothelial
# CCTCATGTCCTTATCA-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female  Endothelial
# TACGCTCTCCATCTCG-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female  Endothelial
# ATTCCTAGTGGAAATT-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female        Oligo
# GCATCGGAGCACCGTC-1   cortex_rna  human_cortex1  NCBITaxon_9606  ...                                    10x 3' v3  female        Oligo

# [4067 rows x 12 columns]

# >>> metadata.columns
# Index(['biosample_id', 'donor_id', 'species', 'species__ontology_label',
#        'disease', 'disease__ontology_label', 'organ', 'organ__ontology_label',
#        'library_preparation_protocol',
#        'library_preparation_protocol__ontology_label', 'sex', 'cluster'],
#       dtype='object')

# find intersect of indices
index_intersect = adata.obs.index\
  .intersection(metadata.index)\
  .intersection(cluster.index)\
  .intersection(spatial.index)\
  .astype(str)

# copy data to 
output = adata[index_intersect].copy()
for metadata_col in metadata.columns:
  output.obs[metadata_col] = metadata.loc[index_intersect, metadata_col]

# output.obs = metadata.loc[index_intersect, metadata.columns].copy()
output.obs["cell_type"] = list(cluster.loc[index_intersect, "cell_type"].values)
output.obsm["X_umap"] = cluster.loc[index_intersect, ["X", "Y"]].values
output.obsm["spatial"] = spatial.loc[index_intersect, ["X", "Y"]].values

# add uns
cols = ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]
for col in cols:
  output.uns[col] = par[col]

# AnnData object with n_obs × n_vars = 4065 × 36601
#     obs: 'biosample_id', 'donor_id', 'species', 'species__ontology_label', 'disease', 'disease__ontology_label', 'organ', 'organ__ontology_label', 'library_preparation_protocol', 'library_preparation_protocol__ontology_label', 'sex', 'cluster', 'cell_type'
#     var: 'gene_ids', 'feature_types'
#     obsm: 'X_umap', 'spatial'

output.write_h5ad(par["output"], compression="gzip")
