//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process bed_to_span {

  tag "$sampleID"
  publishDir "$params.out_dir/spans", mode: "copy"

  input:
  tuple val(sampleID), file(in_bed)

  output:
  tuple val(sampleID),
    file("${sampleID}_spans.txt"),
    emit: spans

  script:
  """
  python $workflow.projectDir/bin/python/bed_to_span.py \
  -b ${in_bed} \
  -o ${sampleID}_spans.txt
  """
}
