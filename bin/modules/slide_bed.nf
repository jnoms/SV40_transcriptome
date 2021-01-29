//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.genome_length = 5243
params.slide_bed_only_keep_wraparound = "no"

//============================================================================//
// Define process
//============================================================================//
process slide_bed {

  tag "$sampleID"
  publishDir "$params.out_dir/slide_bed", mode: "copy"

  // default memory
  memory "10G"

  input:
  tuple val(sampleID), file(in_bed)

  output:
  tuple val(sampleID),
    file("${sampleID}_slid.bed"),
    emit: slid_bed


  script:
  """
  python $workflow.projectDir/bin/python/bed_slide_wraparound_reads.py \
  -i ${in_bed} \
  -o ${sampleID}_slid.bed \
  -g ${params.genome_length} \
  -w ${params.slide_bed_only_keep_wraparound}
  """
}
