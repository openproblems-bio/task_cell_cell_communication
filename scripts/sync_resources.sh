#!/bin/bash

# download data
common/sync_resources/sync_resources \
  --input s3://openproblems-data/resources_test/cell_cell_communication \
  --output resources

aws s3 sync resources \
  s3://openproblems-data/resources/cell_cell_communication \
  --delete \
  --dryrun