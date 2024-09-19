#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

IN_DIR=resources_test/common/singlecell_broadinstitute_scp2167_human_brain
OUT_DIR=resources_test/task_cell_cell_communication/singlecell_broadinstitute_scp2167_human_brain

mkdir -p $OUT_DIR

# infer truth
viash run src/data_processors/infer_truth/config.vsh.yaml -- \
  --input $IN_DIR/dataset.h5ad \
  --output $OUT_DIR/dataset_with_inferred_truth.h5ad

# process dataset
viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --input $OUT_DIR/dataset_with_inferred_truth.h5ad \
  --output_dataset $OUT_DIR/dataset.h5ad \
  --output_solution $OUT_DIR/solution.h5ad

# run one method
viash run src/methods/logistic_regression/config.vsh.yaml -- \
    --input $OUT_DIR/dataset.h5ad \
    --output $OUT_DIR/prediction.h5ad

# run one metric
viash run src/metrics/accuracy/config.vsh.yaml -- \
    --input_prediction $OUT_DIR/prediction.h5ad \
    --input_solution $OUT_DIR/solution.h5ad \
    --output $OUT_DIR/score.h5ad

# write manual state.yaml. this is not actually necessary but you never know it might be useful
cat > $OUT_DIR/state.yaml << HERE
id: singlecell_broadinstitute_scp2167_human_brain
dataset: !file dataset.h5ad
solution: !file solution.h5ad
prediction: !file prediction.h5ad
score: !file score.h5ad
HERE

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  resources_test/task_cell_cell_communication/singlecell_broadinstitute_scp2167_human_brain \
  s3://openproblems-data/resources_test/task_cell_cell_communication/singlecell_broadinstitute_scp2167_human_brain \
  --delete --dryrun
