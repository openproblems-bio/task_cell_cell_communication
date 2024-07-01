workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

      | download_singlecell_broadinstitute_dataset.run(
        fromState: [
          "dataset_curl_config",
          "dataset_id",
          "dataset_name",
          "dataset_url",
          "dataset_reference",
          "dataset_summary",
          "dataset_description",
          "dataset_organism"
        ],
        toState: [
          "raw_dataset": "output"
        ]
      )

      | infer_truth.run(
        fromState: [
          "input": "raw_dataset"
        ],
        toState: [
          "output"
        ]
      )

  emit:
    output_ch
}