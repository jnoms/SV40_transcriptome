#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//============================================================================//
// Set up modules
//============================================================================//
include { Minimap2 } from './bin/modules/Minimap2'
include { STAR } from './bin/modules/STAR'
include { slide_bed } from './bin/modules/slide_bed'
include { bed_to_span } from './bin/modules/bed_to_span'
include { bed_elongate_illumina_reads } from './bin/modules/bed_elongate_illumina_reads'
include { prodigal } from './bin/modules/prodigal'

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
  bed_to_span(slide_bed.out.slid_bed)

  // ------------------------------------------------------------ //
  // ORF ANALYSIS
  // ------------------------------------------------------------ //

  // Elongate bed
  bed_elongate_illumina_reads(slide_bed.out.slid_bed)

  // Predict ORFs with Prodigal
  prodigal(bed_elongate_illumina_reads.out.elongated_fasta)

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
  bed_to_span(slide_bed.out.slid_bed)

  // ------------------------------------------------------------ //
  // ORF ANALYSIS
  // ------------------------------------------------------------ //

  // Predict ORFs with prodigal
  prodigal(aligned.fasta)

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
