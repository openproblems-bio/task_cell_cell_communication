#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

cat > /tmp/params.yaml << HERE
param_list:
  - id: singlecell_broadinstitute_scp2167_human_brain
    raw_data_dir: resources/raw_data/singlecell_broadinstitute/SCP2167
    dataset_name: "Human brain Slide-tags snRNA-seq"
    dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP2167"
    dataset_reference: "10.1038/s41586-023-06837-4"
    dataset_summary: "Slide-tags snRNA-seq data on the human brain."
    dataset_description: "..."
    dataset_organism: "homo_sapiens"

do_subsample: true
n_obs: 800
n_vars: 1000

output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
output_state: "\$id/state.yaml"
publish_dir: resources_test/common
HERE

nextflow run . \
  -main-script target/nextflow/datasets/workflows/process_singlecell_broadinstitute/main.nf \
  -profile docker \
  -resume \
  -params-file /tmp/params.yaml

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  "resources_test/common/singlecell_broadinstitute_scp2167_human_brain" \
  s3://openproblems-data/resources_test/common/singlecell_broadinstitute_scp2167_human_brain \
  --delete --dryrun
