//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process bed_extract_representatives {

  tag "$sampleID"
  publishDir "$params.out_dir/bed_representatives", mode: "copy"

  input:
  tuple val(sampleID), file(in_bed), file(in_spans)

  output:
  tuple val(sampleID),
    file("${sampleID}_representatives.bed"),
    emit: representatives_bed

  script:
  """
  python $workflow.projectDir/bin/python/bed_extract_representatives.py \
  -i ${in_bed} \
  -o ${sampleID}_representatives.bed \
  -s ${in_spans}
  """
}
