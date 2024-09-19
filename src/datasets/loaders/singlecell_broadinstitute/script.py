import scanpy as sc
import pandas as pd
from pathlib import Path
import numpy as np

## VIASH START
par = {
  "raw_data_dir": "resources/raw_data/singlecell_broadinstitute/SCP2167",
  "dataset_id": "singlecell_broadinstitute_scp2167",
  "dataset_name": "--",
  "dataset_url": "--",
  "dataset_reference": "--",
  "dataset_summary": "--",
  "dataset_description": "--",
  "dataset_organism": "homo_sapiens",
  "output": "resources/datasets/SCP2169/raw_dataset.h5ad"
}
meta = {
  "temp_dir": "/tmp"
}
## VIASH END

# temp_dir = f'{meta["temp_dir"]}/downloader_singlecell_broadinstitute_dataset'

# if not os.path.exists(temp_dir):
#   os.makedirs(temp_dir)

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

def get_counts(raw_data_dir):
  matrix_file = next(raw_data_dir.glob("**/matrix.mtx.gz"))

  # read counts
  adata = sc.read_10x_mtx(matrix_file.parent, gex_only=False)
  # AnnData object with n_obs × n_vars = 14165 × 36601
  #     var: 'gene_ids', 'feature_types'

  # remove "-1" from the index
  adata.obs_names = adata.obs_names.str.replace("-1", "")

  return adata

def get_cluster_info(raw_data_dir):
  cluster_info = next(raw_data_dir.glob("**/*_cluster.csv"))
  cluster = read_typed_csv(cluster_info)
  cluster.index = cluster.index.str.replace("-1", "")
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

  return cluster

def get_spatial_info(raw_data_dir):
  spatial_info = next(raw_data_dir.glob("**/*_spatial.csv"))
  spatial = read_typed_csv(spatial_info)
  spatial.index = spatial.index.str.replace("-1", "")
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

  return spatial

def get_metadata_info(raw_data_dir):
  metadata_info = next(raw_data_dir.glob("**/*_metadata.csv"))
  metadata = read_typed_csv(metadata_info)
  metadata.index = metadata.index.str.replace("-1", "")
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
  return metadata

# search for 'matrix.mtx.gz' in the par['raw_data_dir']
raw_data_dir = Path(par["raw_data_dir"])

# read data
adata = get_counts(raw_data_dir)
cluster = get_cluster_info(raw_data_dir)
spatial = get_spatial_info(raw_data_dir)
metadata = get_metadata_info(raw_data_dir)

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
# output.obs["cell_type"] = list(cluster.loc[index_intersect, "cell_type"].values)
# output.obsm["X_umap"] = cluster.loc[index_intersect, ["X", "Y"]].values
output.obsm["spatial"] = spatial.loc[index_intersect, ["X", "Y"]].values

# add uns
cols = ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]
for col in cols:
  output.uns[col] = par[col]


# Fix formatting of the output

# current format:
# AnnData object with n_obs × n_vars = 4065 × 36601
#     obs: 'biosample_id', 'donor_id', 'species', 'species__ontology_label', 'disease', 'disease__ontology_label', 'organ', 'organ__ontology_label', 'library_preparation_protocol', 'library_preparation_protocol__ontology_label', 'sex', 'cluster', 'cell_type'
#     var: 'gene_ids', 'feature_types'
#     obsm: 'X_umap', 'spatial'

# desired format:
# AnnData object with n_obs × n_vars = 584944 × 28024 with slots:
#     obs: soma_joinid, dataset_id, assay, assay_ontology_term_id, cell_type, cell_type_ontology_term_id, development_stage, development_stage_ontology_term_id, disease, disease_ontology_term_id, donor_id, is_primary_data, self_reported_ethnicity, self_reported_ethnicity_ontology_term_id, sex, sex_ontology_term_id, suspension_type, tissue, tissue_ontology_term_id, tissue_general, tissue_general_ontology_term_id, batch, size_factors
#     var: soma_joinid, feature_id, feature_name, hvg, hvg_score
#     obsp: knn_connectivities, knn_distances
#     obsm: X_pca
#     varm: pca_loadings
#     layers: counts, normalized
#     uns: dataset_description, dataset_id, dataset_name, dataset_organism, dataset_reference, dataset_summary, dataset_url, knn, normalization_id, pca_variance

# fix obs
output.obs.rename(
  columns={
    "donor_id": "batch",
    "species": "species_ontology_term_id",
    "species__ontology_label": "species",
    "disease": "disease_ontology_term_id",
    "disease__ontology_label": "disease",
    "organ": "tissue_ontology_term_id",
    "organ__ontology_label": "tissue",
    "library_preparation_protocol": "assay",
    "library_preparation_protocol__ontology_label": "assay_ontology_term_id",
    "cluster": "cell_type",
  },
  inplace=True
)
output.obs["sex_ontology_term_id"] = np.where(output.obs["sex"] == "female", "PATO:0000383", "PATO:0000384")
cell_type_ontology_mapping = {
  "Astrocyte": "CL:0000127",
  "Endothelial": "CL:0000115",
  "Excitatory": "CL:0008030",
  "Inhibitory": "CL:0008029",
  "Microglia": "CL:0000129",
  "Oligo": "CL:0000128",
  "OPC": "CL:0002453"
}
output.obs["cell_type_ontology_term_id"] = output.obs["cell_type"].map(cell_type_ontology_mapping)

# fix var
output.var.index.name = "feature_name"
output.var.reset_index("feature_name", inplace=True)
output.var.rename(columns={"gene_ids": "feature_id"}, inplace=True)
output.var.set_index("feature_id", drop=False, inplace=True)
output.var.drop("feature_types", axis=1, inplace=True)

# move X to counts layer
output.layers = {
  "counts": output.X
}
del output.X

# result
# AnnData object with n_obs × n_vars = 4065 × 36601
#     obs: 'biosample_id', 'batch', 'species_ontology_term_id', 'species', 'disease_ontology_term_id', 'disease', 'tissue_ontology_term_id', 'tissue', 'assay', 'assay_ontology_term_id', 'sex', 'cell_type', 'sex_ontology_term_id', 'cell_type_ontology_term_id'
#     var: 'feature_name', 'feature_id'
#     uns: 'dataset_id', 'dataset_name', 'dataset_url', 'dataset_reference', 'dataset_summary', 'dataset_description', 'dataset_organism'
#     obsm: 'spatial'
#     layers: 'counts'

# write to file
output.write_h5ad(par["output"], compression="gzip")
