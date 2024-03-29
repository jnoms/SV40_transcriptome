#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//============================================================================//
// Set up modules
//============================================================================//
include { Minimap2 } from './bin/modules/Minimap2'
include { STAR } from './bin/modules/STAR'
include { slide_bed } from './bin/modules/slide_bed'
include { bed_elongate_illumina_reads } from './bin/modules/bed_elongate_illumina_reads'
include { prodigal } from './bin/modules/prodigal'
include { filter_prodigal } from './bin/modules/filter_prodigal'
include { prodigal_to_orfs_direct } from './bin/modules/prodigal_to_orfs_direct'
include { diamond } from './bin/modules/diamond'
include { bed_extract_representatives } from './bin/modules/bed_extract_representatives'
include { characterize_ORFs } from './bin/modules/characterize_ORFs'
include { bam_coverage } from './bin/modules/bam_coverage'

// Modules with different param settings
include { bed_to_span as bed_to_span_ILLUMINA } from './bin/modules/bed_to_span' \
  addParams(bed_to_span_junc_reads_only: params.bed_to_span_junc_reads_only_ILLUMINA)
include { bed_to_span as bed_to_span_NANOPORE } from './bin/modules/bed_to_span' \
  addParams(bed_to_span_junc_reads_only: params.bed_to_span_junc_reads_only_NANOPORE)


// Note: Add addParams(param_name: 'desired value') to end of include line to
// specify a param value for that module. OTHERWISE, each param is taken
// from the value in this script automatically (or presumably the config).
// E.g. by default it does similar behavior as params(params)
// EXAMPLE: include { slide_bed as slide_bed_illumina } from './bin/modules/slide_bed' \
//  addParams(slide_bed_only_keep_wraparound: params.ILLUMINA_slide_bed_only_keep_wraparound)

// NOTE - use .view{"$it"} to inspect a chanel.. this will print out the tuples.


//============================================================================//
// Defining functions
//============================================================================//
def sampleID_set_from_infile(input) {
  // The purpose of this function is to take an input list of file paths and
  // to return a channel of sets of structure basename:file_path.
  sample_set = []
  for (String item : file(input)) {
    file = file(item)
    name = file.baseName

    // This will handle .fastq.gz and .fq.gz... baseName above only removes
    // the last suffix.
    if (name.endsWith(".fastq") | name.endsWith(".fq") )
      name = name.take(name.lastIndexOf('.'))

    sample_set.add([name, file])
  }
  ch = Channel.from(tuple(sample_set))
  return ch
}

//============================================================================//
// Define workflows
//============================================================================//
workflow illumina {

  take: input_ch
  main:

  // Map using minimap or star
  if (params.aligner == "Minimap2")
    aligned = Minimap2(input_ch)

  else if (params.aligner == "STAR")
    aligned = STAR(input_ch)

  // Slide the bed
  slide_bed(aligned.bed)

  // Generate spans
  bed_to_span_ILLUMINA(slide_bed.out.slid_bed)

  // Calculate coverage
  bam_coverage(aligned.bam)

  // ------------------------------------------------------------ //
  // ORF ANALYSIS
  // ------------------------------------------------------------ //

  // Extract representative transcripts from the slid bed
  slide_bed.out.slid_bed
    .join(bed_to_span_ILLUMINA.out.spans) |\
    bed_extract_representatives

  // Elongate bed
  bed_elongate_illumina_reads(bed_extract_representatives.out.representatives_bed)

  // Predict ORFs with Prodigal
  prodigal(bed_elongate_illumina_reads.out.elongated_fasta)

  // Filter prodigal for only spliced ORFs - need to merge with bed file first.
  prodigal.out.prodigal_out
    .join(bed_elongate_illumina_reads.out.elongated_bed) |\
    filter_prodigal

  // Extract ORFs
  bed_elongate_illumina_reads.out.elongated_fasta
    .join(filter_prodigal.out.filtered_prodigal) |\
    prodigal_to_orfs_direct

  // Align with diamond
  prodigal_to_orfs_direct.out.pr_orfs |\
    diamond

  // Generate ORF report
  diamond.out.diamond_out
    .join(bed_elongate_illumina_reads.out.elongated_bed)
    .join(bed_to_span_ILLUMINA.out.spans) |\
    characterize_ORFs
}

workflow nanopore {

  take: input_ch
  main:

  // Map using minimap
  if (params.aligner == "Minimap2")
    aligned = Minimap2(input_ch)
  else if (params.aligner == "STAR")
    error "For Nanopore reads, MUST use Minimap2 aligner."

  // Slide the bed
  slide_bed(aligned.bed)

  // Generate spans
  bed_to_span_NANOPORE(slide_bed.out.slid_bed)

  // Calculate coverage
  bam_coverage(aligned.bam)

  // ------------------------------------------------------------ //
  // ORF ANALYSIS
  // ------------------------------------------------------------ //

  // Predict ORFs with prodigal
  prodigal(aligned.fasta)

  // Extract ORFs
  aligned.fasta
    .join(prodigal.out.prodigal_out) |\
    prodigal_to_orfs_direct

  // Align with diamond
  prodigal_to_orfs_direct.out.pr_orfs |\
    diamond

  // Generate ORF report
  diamond.out.diamond_out
    .join(slide_bed.out.slid_bed)
    .join(bed_to_span_NANOPORE.out.spans) |\
    characterize_ORFs
}



//============================================================================//
// Validate inputs
//============================================================================//
if( (params.in_fastq_type != "Illumina") && (params.in_fastq_type != "Nanopore")) {
  error "params.datatype must be set to 'Illumina' or 'Nanopore'."
}

if( (params.aligner != "Minimap2") && (params.aligner != "STAR")) {
  error "params.aligner must be set to 'Minimap2' or 'STAR'."
}

if( (params.only_primary_alignments != "yes") && (params.only_primary_alignments != "no")) {
  error "params.aligner must be set to 'yes' or 'no'."
}

if( (params.slide_bed_only_keep_wraparound != "yes") && (params.slide_bed_only_keep_wraparound != "no")) {
  error "params.slide_bed_only_keep_wraparound must be set to 'yes' or 'no'."
}

if( (params.illumina_bed_elongate_only_keep_spliced_reads != "yes") && (params.illumina_bed_elongate_only_keep_spliced_reads != "no")) {
  error "params.bed_elongate_only_keep_spliced_reads must be set to 'yes' or 'no'."
}

if( (params.bed_to_span_junc_reads_only_ILLUMINA != "yes") && (params.bed_to_span_junc_reads_only_ILLUMINA != "no")) {
  error "params.bed_to_span_junc_reads_only_ILLUMINA must be set to 'yes' or 'no'."
}

if( (params.bed_to_span_junc_reads_only_NANOPORE != "yes") && (params.bed_to_span_junc_reads_only_NANOPORE != "no")) {
  error "params.bed_to_span_junc_reads_only_NANOPORE must be set to 'yes' or 'no'."
}

//============================================================================//
// Define main workflow
//============================================================================//
workflow {

  main:
    input_ch = sampleID_set_from_infile(params.in_fastq)

    if ( params.in_fastq_type == "Illumina" )
      illumina(input_ch)

    else if ( params.in_fastq_type == "Nanopore" )
      nanopore(input_ch)
}
