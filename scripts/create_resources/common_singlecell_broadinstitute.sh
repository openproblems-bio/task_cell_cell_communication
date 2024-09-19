#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# # remove this when you have implemented the script
# echo "TODO: replace the commands in this script with the sequence of components that you need to run to generate test_resources."
# echo "  Inside this script, you will need to place commands to generate example files for each of the 'src/api/file_*.yaml' files."
# exit 1

set -e

RAW_DATA=resources/raw_data/singlecell_broadinstitute
DATASET_ID=singlecell_broadinstitute_scp2167_human_brain
DATASET_DIR=resources_test/common

mkdir -p $DATASET_DIR

# turn raw data into spatial anndata file
viash run src/datasets/loaders/singlecell_broadinstitute/config.vsh.yaml -- \
  --raw_data_dir "$RAW_DATA/SCP2167" \
  --dataset_id $DATASET_ID \
  --dataset_name "Human brain Slide-tags snRNA-seq" \
  --dataset_url "https://singlecell.broadinstitute.org/single_cell/study/SCP2167" \
  --dataset_reference "10.1038/s41586-023-06837-4" \
  --dataset_summary "Slide-tags snRNA-seq data on the human brain." \
  --dataset_description "..." \
  --dataset_organism homo_sapiens \
  --output "$DATASET_DIR/$DATASET_ID/dataset.h5ad"

# write manual state.yaml. this is not actually necessary but you never know it might be useful
cat > $DATASET_DIR/$DATASET_ID/state.yaml << HERE
id: $DATASET_ID
dataset: !file dataset.h5ad
HERE

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  "resources_test/common/singlecell_broadinstitute_scp2167_human_brain" \
  s3://openproblems-data/resources_test/common/singlecell_broadinstitute_scp2167_human_brain \
  --delete --dryrun
