//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.diamond_database = "$baseDir/resources/ref/SV40/SV40_canonical_pr.dmnd"
params.diamond_outfmt = "6 qseqid qlen sseqid evalue bitscore pident length qstart qend sstart send"
params.diamond_temp_dir = "temp"
params.diamond_evalue = "10"
params.diamond_max_alignments = "0" // 0 - all alignments
params.diamond_include_unaligned = "TRUE" // options TRUE or FALSE

//============================================================================//
// Define process
//============================================================================//
process diamond {
  tag "$sampleID"
  publishDir "$params.out_dir/diamond", mode: "copy"

  // default memory
  memory "10G"

  input:
  tuple val(sampleID), file(pr_fasta)

  output:
  tuple val(sampleID),
    file("${sampleID}_diamond.out"),
    emit: diamond_out

  script:
  """
  $workflow.projectDir/bin/bash/diamond.sh \
  -d ${params.diamond_database} \
  -q ${pr_fasta} \
  -o ${sampleID}_diamond.out \
  -k ${params.diamond_max_alignments} \
  -m ${task.memory.toGiga()} \
  -t ${params.diamond_temp_dir} \
  -e ${params.diamond_evalue} \
  -f "${params.diamond_outfmt}" \
  -s ${sampleID} \
  -u "${params.diamond_include_unaligned}"
  """
}
