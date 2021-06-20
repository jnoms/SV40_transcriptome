//============================================================================//
// Default params
//============================================================================//
params.out_dir = "output"
params.bed_to_span_junc_reads_only = "no"
params.bed_to_span_illumina_spans = ""
params.n_genomes = 1
params.bed_to_span_n_illumina_junctions = 5
params.genome_length = 0
params.bed_to_span_keep_potential_duplicates = "no"

//============================================================================//
// Define process
//============================================================================//
process bed_to_span {

  tag "$sampleID"
  publishDir "$params.out_dir/spans", mode: "copy"

  input:
  tuple val(sampleID), file(in_bed)

  output:
  tuple val(sampleID),
    file("${sampleID}_spans.txt"),
    emit: spans

  tuple val(sampleID),
    file("${sampleID}_UNSUPPORTED_spans.txt"),
    emit: UNSUPPORTED_spans

  script:
  """
  python $workflow.projectDir/bin/python/bed_to_span.py \
  -b ${in_bed} \
  -o ${sampleID}_spans.txt \
  -j ${params.bed_to_span_junc_reads_only} \
  -i "${params.bed_to_span_illumina_spans}" \
  -g ${params.n_genomes} \
  -u ${params.bed_to_span_n_illumina_junctions} \
  -l ${params.genome_length} \
  -x ${sampleID}_UNSUPPORTED_spans.txt \
  -d ${params.bed_to_span_keep_potential_duplicates}

  # Make an empty output file if not specified
  # Ideally, would make this an optional output...
  if [[ ! -f ${sampleID}_UNSUPPORTED_spans.txt ]] ; then
    touch ${sampleID}_UNSUPPORTED_spans.txt
  fi
  """
}
