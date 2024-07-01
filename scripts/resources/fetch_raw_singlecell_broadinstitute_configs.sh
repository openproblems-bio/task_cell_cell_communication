#!/bin/bash

# Go to https://www.nature.com/articles/s41586-023-06837-4#data-availability and for each of the SCP21xx datasets, open up the link.
# Click on the 'Bulk download' button and copy the URL.

CONFIG_DIR=resources/datasets_raw/singlecell_broadinstitute_configs

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2162&auth_code=TsNVrxUf&directory=all&context=study" -o $CONFIG_DIR/SCP2162.txt

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2170&auth_code=kEQJBwiY&directory=all&context=study" -o $CONFIG_DIR/SCP2170.txt

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2167&auth_code=Fy69UsMs&directory=all&context=study" -o $CONFIG_DIR/SCP2167.txt

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2169&auth_code=LXe9HEwf&directory=all&context=study" -o $CONFIG_DIR/SCP2169.txt

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2171&auth_code=hUnl6u1f&directory=all&context=study" -o $CONFIG_DIR/SCP2171.txt

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2176&auth_code=RkOJsvTu&directory=all&context=study" -o $CONFIG_DIR/SCP2176.txt
