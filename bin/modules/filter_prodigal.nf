//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.filter_prodigal_allowable_strands = "+"

//============================================================================//
// Define process
//============================================================================//
process filter_prodigal {

  tag "$sampleID"
  publishDir "$params.out_dir/filter_prodigal", mode: "copy"

  input:
  tuple val(sampleID), file(in_prodigal), file(in_bed)

  output:
  tuple val(sampleID),
    file("${sampleID}_prodigal_filtered.txt.gz"),
    emit: filtered_prodigal

  script:
  """
  python $workflow.projectDir/bin/python/filter_prodigal.py \
  -p ${in_prodigal} \
  -b ${in_bed} \
  -o ${sampleID}_prodigal_filtered.txt.gz \
  -s ${params.filter_prodigal_allowable_strands}
  """
}
