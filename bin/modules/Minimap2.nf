//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.reference_genome_fasta = "$baseDir/resources/ref/SV40/SV40_100_DOUBLED.fasta"
params.only_primary_alignments = "yes"

//============================================================================//
// Define process
//============================================================================//
process Minimap2 {

  tag "$sampleID"
  publishDir "$params.out_dir/Minimap2", mode: "copy"

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

  tuple val(sampleID),
    file("${sampleID}_DERVIED.fasta"),
    emit: fasta

  script:
  """
  $workflow.projectDir/bin/bash/minimap2.sh \
  -i ${in_fastq} \
  -r ${params.reference_genome_fasta} \
  -o ${sampleID}.bam \
  -b ${sampleID}.bed \
  -f ${sampleID}_DERVIED.fasta \
  -P ${params.only_primary_alignments} \
  -t ${task.cpus}
  """
}
