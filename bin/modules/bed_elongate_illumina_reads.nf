//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.reference_length = 10486
params.reference_genome_fasta = "$baseDir/resources/ref/SV40/SV40_100_DOUBLED.fasta"
params.illumina_bed_elongate_only_keep_spliced_reads = "yes"

//============================================================================//
// Define process
//============================================================================//
process bed_elongate_illumina_reads {

  tag "$sampleID"
  publishDir "$params.out_dir/elongated", mode: "copy"

  input:
  tuple val(sampleID), file(in_bed)

  output:
  tuple val(sampleID),
    file("${sampleID}_elongated.bed"),
    emit: elongated_bed
  tuple val(sampleID),
    file("${sampleID}_elongated.fasta"),
    emit: elongated_fasta

  script:
  """
  # Determine spliced_reads setting
  if [[ $params.illumina_bed_elongate_only_keep_spliced_reads == "yes" ]] ; then
    SPLICED_READS_SETTING="-s"
  else
    SPLICED_READS_SETTING=""
  fi

  python $workflow.projectDir/bin/python/bed_elongate_illumina_reads.py \
  \$SPLICED_READS_SETTING \
  -i ${in_bed} \
  -o ${sampleID}_elongated.bed \
  -r ${params.reference_length}

  # Use bedtools to generate elongated fasta
  bedtools getfasta \
  -split \
  -name \
  -s \
  -fi ${params.reference_genome_fasta} \
  -bed ${sampleID}_elongated.bed > ${sampleID}_elongated.fasta
  """
}
