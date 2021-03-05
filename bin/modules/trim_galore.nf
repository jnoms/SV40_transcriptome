//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.trim_galore_run_fastqc = "T"
params.trim_galore_minlength = 40

//============================================================================//
// Define process
//============================================================================//
process trim_galore {

  tag "$sampleID"
  publishDir "$params.out_dir/trim_galore", mode: "copy"

  input:
  tuple val(sampleID), file(read1), file(read2)

  output:
  tuple val(sampleID),
    file("trimmed/${sampleID}_r1.fastq.gz"),
    file("trimmed/${sampleID}_r2.fastq.gz"),
    emit: reads

  tuple val(sampleID),
    file("trimmed/reports/*"),
    emit: reports

  script:
  """
  $workflow.projectDir/bin/bash/trim_galore.sh \
  -1 $read1 \
  -2 $read2 \
  -3 trimmed/${sampleID}_r1.fastq.gz \
  -4 trimmed/${sampleID}_r2.fastq.gz \
  -j ${task.cpus} \
  -f ${params.trim_galore_run_fastqc} \
  -l ${params.trim_galore_minlength}
  """
}
