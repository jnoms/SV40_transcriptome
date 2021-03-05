#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//============================================================================//
// Set up modules
//============================================================================//
include { trim_galore } from './bin/modules/trim_galore'
include { STAR } from './bin/modules/STAR'


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
// Define main workflow
//============================================================================//
workflow {

  main:

  // Parse csv input
  Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sample_name, file(row.read1), file(row.read2)) }
    .set { input_ch }

  // Trim with trim_galore
  trim_galore(input_ch)


}
