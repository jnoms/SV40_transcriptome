//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.genome_length = 5243

//============================================================================//
// Define process
//============================================================================//
process bam_coverage {

  tag "$sampleID"
  publishDir "$params.out_dir/bam_coverage", mode: "copy"

  input:
  tuple val(sampleID), file(in_bam)

  output:
  tuple val(sampleID),
    file("${sampleID}_slid.bam"),
    file("${sampleID}_slid.bam.bai"),
    emit: slid_bam

  tuple val(sampleID),
    file("${sampleID}_slid_f.cov"),
    file("${sampleID}_slid_r.cov"),
    emit: coverage


  script:
  """

  # Slide the bam
  python $workflow.projectDir/bin/python/bam_slide_wraparound_reads.py \
  -i ${in_bam} \
  -o ${sampleID}_slid.bam \
  -g ${params.genome_length}

  # Split to forward and reverse
  samtools view -bh -F 16 -@ ${task.cpus} \
  -o ${sampleID}_slid_f.bam_tmp \
  -U ${sampleID}_slid_r.bam_tmp

  # Sort and index
  samtools sort -@ ${task.cpus} \
    ${sampleID}_slid_r.bam_tmp > ${sampleID}_slid_r.bam && \
    rm ${sampleID}_slid_r.bam_tmp
  samtools sort -@ ${task.cpus} \
    ${sampleID}_slid_f.bam_tmp > ${sampleID}_slid_f.bam && \
    rm ${sampleID}_slid_f.bam_tmp

  samtools index ${sampleID}_slid_r.bam
  samtools index ${sampleID}_slid_r.bam

  # Calculate depth
  samtools depth -aa -d0 ${sampleID}_slid_r.bam > ${sampleID}_slid_r.cov && \
    rm ${sampleID}_slid_r.bam
  samtools depth -aa -d0 ${sampleID}_slid_f.bam > ${sampleID}_slid_f.cov && \
    rm ${sampleID}_slid_f.bam

  """
}
