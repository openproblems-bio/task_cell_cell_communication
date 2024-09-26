workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch

      | infer_truth.run(
        fromState: [
          "input"
        ],
        toState: [
          "dataset_with_groundtruth"
        ]
      )

      | process_datasets.run(
        fromState: [
          "input": "dataset_with_groundtruth"
        ],
        toState: [
          "output_dataset",
          "output_solution"
        ]
      )

      | setState(["output_dataset", "output_solution"])

  emit:
    output_ch
}