//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"

//============================================================================//
// Define process
//============================================================================//
process prodigal_to_orfs_direct {
  tag "$sampleID"
  publishDir "$params.out_dir/prodigal_to_orfs_direct", mode: "copy"

  input:
  tuple val(sampleID), file(in_fasta), file(prodigal_file)

  output:
  tuple val(sampleID),
    file("${sampleID}_pr.fasta"),
    emit: pr_orfs

  tuple val(sampleID),
    file("${sampleID}_nt.fasta"),
    emit: nt_orfs

  script:
  """
  python $workflow.projectDir/bin/python/prodigal_to_orfs_direct.py \
  -f ${in_fasta} \
  -p ${prodigal_file} \
  -N ${sampleID}_nt.fasta \
  -P ${sampleID}_pr.fasta
  """
}
