#!/bin/bash

# Go to https://www.nature.com/articles/s41586-023-06837-4#data-availability and for each of the SCP21xx datasets, open up the link.
# Click on the 'Bulk download' button and copy the URL.

echo "Run these commands from the terminal line-by-line to download the raw data." && exit 1

CONFIG_DIR=resources/raw_data/singlecell_broadinstitute
RAW_DIR="s3://openproblems-data/resources/raw_data/singlecell_broadinstitute"

curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2162&auth_code=cC5TdSLA&directory=all&context=study" -o $CONFIG_DIR/SCP2162.txt
pushd $CONFIG_DIR; curl -k -K SCP2162.txt; rm SCP2162.txt; popd

curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2170&auth_code=YJ0YWhmE&directory=all&context=study" -o $CONFIG_DIR/SCP2170.txt
pushd $CONFIG_DIR; curl -k -K SCP2170.txt; rm SCP2170.txt; popd

curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2167&auth_code=UMp0hscZ&directory=all&context=study" -o $CONFIG_DIR/SCP2167.txt
pushd $CONFIG_DIR; curl -k -K SCP2167.txt; rm SCP2167.txt; popd

curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2169&auth_code=SHKALG5D&directory=all&context=study" -o $CONFIG_DIR/SCP2169.txt
pushd $CONFIG_DIR; curl -k -K SCP2169.txt; rm SCP2169.txt; popd

curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2171&auth_code=Zupo1qmK&directory=all&context=study" -o $CONFIG_DIR/SCP2171.txt
pushd $CONFIG_DIR; curl -k -K SCP2171.txt; rm SCP2171.txt; popd

curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2176&auth_code=IEarrguR&directory=all&context=study" -o $CONFIG_DIR/SCP2176.txt
pushd $CONFIG_DIR; curl -k -K SCP2176.txt; rm SCP2176.txt; popd

aws s3 sync --profile op \
  resources/raw_data/singlecell_broadinstitute \
  s3://openproblems-data/resources/raw_data/singlecell_broadinstitute \
  --delete --dryrun
