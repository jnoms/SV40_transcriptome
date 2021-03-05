//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process reverse_complement_read1 {

  tag "$sampleID"
  publishDir "$params.out_dir/reverse_complement_read1", mode: "copy"

  input:
  tuple val(sampleID), file(read1), file(read2)

  output:
  tuple val(sampleID),
    file("${sampleID}_r1_rev_comp.fastq.gz"),
    file("${read2}"),
    emit: reads

  script:
  """
  python $workflow.projectDir/bin/python/fastq_reverse_complement.py \
  -i $read1 \
  -o ${sampleID}_r1_rev_comp.fastq.gz
  """
}
