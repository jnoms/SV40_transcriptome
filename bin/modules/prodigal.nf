//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process prodigal {

  tag "$sampleID"
  publishDir "$params.out_dir/prodigal", mode: "copy"

  // default memory
  memory "10G"

  input:
  tuple val(sampleID), file(in_fasta)

  output:
  tuple val(sampleID),
    file("${sampleID}_prodigal.txt.gz"),
    emit: prodigal_out


  script:
  """
  $workflow.projectDir/bin/bash/prodigal.sh \
  -i ${in_fasta} \
  -o ${sampleID}_prodigal.txt.gz \
  -c "TRUE"
  """
}
