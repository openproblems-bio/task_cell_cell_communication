#!/bin/bash

# Go to https://www.nature.com/articles/s41586-023-06837-4#data-availability and for each of the SCP21xx datasets, open up the link.
# Click on the 'Bulk download' button and copy the URL.

CONFIG_DIR=resources/datasets/raw/singlecell_broadinstitute_configs

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2162&auth_code=xxxxxxx&directory=all&context=study" -o $CONFIG_DIRSCP2162.txt

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2170&auth_code=xxxxxxx&directory=all&context=study" -o $CONFIG_DIRSCP2170.txt

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2167&auth_code=xxxxxxx&directory=all&context=study" -o $CONFIG_DIRSCP2167.txt

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2169&auth_code=xxxxxxx&directory=all&context=study" -o $CONFIG_DIRSCP2169.txt

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2171&auth_code=xxxxxxx&directory=all&context=study" -o $CONFIG_DIRSCP2171.txt

# curl -k "https://singlecell.broadinstitute.org/single_cell/api/v1/bulk_download/generate_curl_config?accessions=SCP2176&auth_code=xxxxxxx&directory=all&context=study" -o $CONFIG_DIRSCP2176.txt
