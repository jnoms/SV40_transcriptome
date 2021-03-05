//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process label_illumina_reads {

  tag "$sampleID"
  publishDir "$params.out_dir/processed_reads", mode: "copy"

  input:
  tuple val(sampleID), file(read1), file(read2)

  output:
  tuple val(sampleID),
    file("${sampleID}.fq.gz"),
    emit: processed_reads

  script:
  """
  python $workflow.projectDir/bin/python/fastq_label_read_number.py \
  -i $read1 \
  -o ${sampleID}_r1_labeled.fq.gz \
  -n 1

  python $workflow.projectDir/bin/python/fastq_label_read_number.py \
  -i $read2 \
  -o ${sampleID}_r2_labeled.fq.gz \
  -n 2

  # Concatenate them - this is the output!
  cat ${sampleID}_r1_labeled.fq.gz ${sampleID}_r2_labeled.fq.gz  >
    ${sampleID}.fq.gz
  """
}
