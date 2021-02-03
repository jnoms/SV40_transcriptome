//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.diamond_outfmt = "6 qseqid sseqid evalue bitscore pident length qlen slen qstart qend sstart send"

//============================================================================//
// Define process
//============================================================================//
process characterize_ORFs {

  tag "$sampleID"
  publishDir "$params.out_dir/characterize_ORFs", mode: "copy"

  input:
  tuple val(sampleID),
    file(in_diamond),
    file(in_bed),
    file(in_spans)

  output:
  tuple val(sampleID),
    file("${sampleID}_ORFs.txt"),
    emit: ORF_report

  script:
  """
  python $workflow.projectDir/bin/python/characterize_ORFs.py \
  -d ${in_diamond} \
  -b ${in_bed} \
  -s ${in_spans} \
  -o ${sampleID}_ORFs.txt \
  -D "${params.diamond_outfmt}"
  """
}
