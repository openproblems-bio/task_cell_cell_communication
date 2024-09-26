#!/bin/bash

RAW_DIR="s3://openproblems-data/resources/raw_data/singlecell_broadinstitute"
OUT_DIR="s3://openproblems-data/resources/datasets"

cat > /tmp/params.yaml << HERE
param_list:
  # - id: singlecell_broadinstitute_scp2162
  #   dataset_curl_config: resources/raw_data/singlecell_broadinstitute_configs/SCP2162.txt
  #   dataset_id: singlecell_broadinstitute_scp2162
  #   dataset_name: "Mouse brain snRNA-seq"
  #   dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP2162"
  #   dataset_reference: "doi: 10.1038/s41586-023-06837-4"
  #   dataset_summary: "snRNA-seq data on the mouse brain."
  #   dataset_description: "Recent technological innovations have enabled the high-throughput quantification of gene expression and epigenetic regulation within individual cells, transforming our understanding of how complex tissues are constructed. Missing from these measurements, however, is the ability to routinely and easily spatially localise these profiled cells. We developed a strategy, Slide-tags, in which single nuclei within an intact tissue section are tagged with spatial barcode oligonucleotides derived from DNA-barcoded beads with known positions. These tagged nuclei can then be used as input into a wide variety of single-nucleus profiling assays. We used Slide-tags to profile two different stages of development in the mouse brain. Overall design 	Slide-tags was used to spatially barcode nuclei from 20-micron thick fresh frozen tissue sections. These spatially barcoded nuclei were then used as input for 10x Genomics Chromium v3.1 snRNA-seq."
  #   dataset_organism: "mus_musculus"
  - id: singlecell_broadinstitute/scp2167_human_brain
    dataset_id: singlecell_broadinstitute_scp2167
    raw_data_dir: $RAW_DIR/SCP2167
    dataset_name: "Human brain snRNA-seq"
    dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP2167"
    dataset_reference: "doi: 10.1038/s41586-023-06837-4"
    dataset_summary: "snRNA-seq data on the human brain."
    dataset_description: "Recent technological innovations have enabled the high-throughput quantification of gene expression and epigenetic regulation within individual cells, transforming our understanding of how complex tissues are constructed. Missing from these measurements, however, is the ability to routinely and easily spatially localise these profiled cells. We developed a strategy, Slide-tags, in which single nuclei within an intact tissue section are tagged with spatial barcode oligonucleotides derived from DNA-barcoded beads with known positions. These tagged nuclei can then be used as input into a wide variety of single-nucleus profiling assays. We used Slide-tags to profile two different stages of development in the mouse brain. Overall design 	Slide-tags was used to spatially barcode nuclei from 20-micron thick fresh frozen tissue sections. These spatially barcoded nuclei were then used as input for 10x Genomics Chromium v3.1 snRNA-seq."
    dataset_organism: "homo_sapiens"
  - id: singlecell_broadinstitute/scp2169_human_tonsil
    dataset_id: singlecell_broadinstitute/scp2169_human_tonsil
    raw_data_dir: $RAW_DIR/SCP2169
    dataset_name: "Human tonsil Slide-tags snRNA-seq"
    dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP2169"
    dataset_reference: "doi: 10.1038/s41586-023-06837-4"
    dataset_summary: "Slide-tags snRNA-seq data on the human tonsil."
    dataset_description: "Recent technological innovations have enabled the high-throughput quantification of gene expression and epigenetic regulation within individual cells, transforming our understanding of how complex tissues are constructed. Missing from these measurements, however, is the ability to routinely and easily spatially localise these profiled cells. We developed a strategy, Slide-tags, in which single nuclei within an intact tissue section are tagged with spatial barcode oligonucleotides derived from DNA-barcoded beads with known positions. These tagged nuclei can then be used as input into a wide variety of single-nucleus profiling assays. We used Slide-tags to profile two different stages of development in the mouse brain. Overall design 	Slide-tags was used to spatially barcode nuclei from 20-micron thick fresh frozen tissue sections. These spatially barcoded nuclei were then used as input for 10x Genomics Chromium v3.1 snRNA-seq."
    dataset_organism: "homo_sapiens"
  - id: singlecell_broadinstitute/scp2170_mouse_embryonic_brain
    dataset_id: singlecell_broadinstitute/scp2170_mouse_embryonic_brain
    raw_data_dir: $RAW_DIR/SCP2170
    dataset_name: "Mouse embryonic brain Slide-tags snRNA-seq"
    dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP2170"
    dataset_reference: "doi: 10.1038/s41586-023-06837-4"
    dataset_summary: "Slide-tags snRNA-seq data on the mouse embryonic brain."
    dataset_description: "Recent technological innovations have enabled the high-throughput quantification of gene expression and epigenetic regulation within individual cells, transforming our understanding of how complex tissues are constructed. Missing from these measurements, however, is the ability to routinely and easily spatially localise these profiled cells. We developed a strategy, Slide-tags, in which single nuclei within an intact tissue section are tagged with spatial barcode oligonucleotides derived from DNA-barcoded beads with known positions. These tagged nuclei can then be used as input into a wide variety of single-nucleus profiling assays. We used Slide-tags to profile two different stages of development in the mouse brain. Overall design 	Slide-tags was used to spatially barcode nuclei from 20-micron thick fresh frozen tissue sections. These spatially barcoded nuclei were then used as input for 10x Genomics Chromium v3.1 snRNA-seq."
    dataset_organism: "mus_musculus"
  - id: singlecell_broadinstitute/scp2171_human_melanoma
    dataset_id: singlecell_broadinstitute/scp2171_human_melanoma
    raw_data_dir: $RAW_DIR/SCP2171
    dataset_name: "Human melanoma Slide-tags snRNA-seq"
    dataset_url: "https://singlecell.broadinstitute.org/single_cell/study/SCP2171"
    dataset_reference: "doi: 10.1038/s41586-023-06837-4"
    dataset_summary: "Slide-tags snRNA-seq data on the human melanoma."
    dataset_description: "Recent technological innovations have enabled the high-throughput quantification of gene expression and epigenetic regulation within individual cells, transforming our understanding of how complex tissues are constructed. Missing from these measurements, however, is the ability to routinely and easily spatially localise these profiled cells. We developed a strategy, Slide-tags, in which single nuclei within an intact tissue section are tagged with spatial barcode oligonucleotides derived from DNA-barcoded beads with known positions. These tagged nuclei can then be used as input into a wide variety of single-nucleus profiling assays. We used Slide-tags to profile two different stages of development in the mouse brain. Overall design 	Slide-tags was used to spatially barcode nuclei from 20-micron thick fresh frozen tissue sections. These spatially barcoded nuclei were then used as input for 10x Genomics Chromium v3.1 snRNA-seq."
    dataset_organism: "homo_sapiens"

output: "\$id/dataset.h5ad"
output_state: "\$id/state.yaml"
# for when this becomes a proper wf:
output_dataset: "\$id/dataset.h5ad"
output_meta: "\$id/dataset_meta.yaml"
publish_dir: $OUT_DIR
HERE

nextflow run . \
  -main-script target/nextflow/dataset_loaders/workflow/main.nf \
  -latest \
  -resume \
  -profile docker \
  -params-file /tmp/params.yaml

tw launch https://github.com/openproblems-bio/task_cell_cell_communication.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/datasets/loaders/singlecell_broadinstitute/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --entry-name auto \
  --config common/nextflow_helpers/labels_tw.config \
  --labels task_template,full