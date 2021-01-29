//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.only_primary_alignments = "yes"
params.STAR_reference_index = "$baseDir/resources/ref/SV40/STAR_SV40_100_DOUBLED"

//============================================================================//
// Define process
//============================================================================//
process STAR {

  tag "$sampleID"
  publishDir "$params.out_dir/STAR", mode: "copy"

  // default memory
  memory "10G"

  input:
  tuple val(sampleID), file(in_fastq)

  output:
  tuple val(sampleID),
    file("${sampleID}.bam"),
    emit: bam

  tuple val(sampleID),
    file("${sampleID}.bed"),
    emit: bed

  script:
  """
  $workflow.projectDir/bin/bash/star.sh \
  -i ${in_fastq} \
  -r ${params.STAR_reference_index} \
  -o ${sampleID}.bam \
  -b ${sampleID}.bed \
  -f ${sampleID}_DERVIED.fasta \
  -p ${params.only_primary_alignments} \
  -t ${task.cpus}
  """
}
