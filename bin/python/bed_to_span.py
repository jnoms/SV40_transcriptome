#!/usr/bin/env python3

import argparse
import pathlib
import os
from collections import Counter
import sys

def get_args():
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to generate a tidy 'span' file from a
            bed file. This span file will detail, for each transcript, the start
            and stop coodinates of each exon/intron (one per read) as well as
            the tx_start and tx_end (tx == transcript). Furthermore, this script
            will cluster reads by their junctions into tx_classess, and will
            report the counts of reads in each tx_class.

            NOTE - generally you should disable supplementary read alignments
            for the alignment step used to generate the BED file. Take notice of
            the paired_end_input input parameter, as this will influence how
            multiple identical read names are handled.
            """)

    # Required arguments
    parser.add_argument(
        '-b',
        '--bed_path',
        type=str,
        required=True,
        help='''
        Path to the BED file. This file should be generated in the following
        way:
        - Map RNAseq reads (illumina or nanopore dRNAseq) against the viral
        transcriptome using minimap2.
            - Generally good idea to not report secondary alignments, else
              they will end up in the bed file.
        - Convert the bam file to BED using bedtools.
        '''
    )
    parser.add_argument(
        '-o',
        '--output_spans',
        type=str,
        required=True,
        help='''
        Path to a tidy output file containing information on read spans. If a
        transcript has multiple exons it will occupy multiple lines, one line
        per exon or intron. The file will have the following columns:
        - name: name of the transcript
        - start: start of the exon or intron
        - end: end of the exon or intron
        - strand: strand of the transcript
        - span_type: either "exon" or "intron"
        - tx_start: start of the transcript
        - tx_end: end of the transcript
        - tx_class: interger value of transcript group. This is ascending in
        order of greater number of transcripts in the class
        - tx_class_count: number of transcripts in the transcript class
        '''
    )

    # Optional arguments
    parser.add_argument(
        '-i',
        '--invert_strand',
        type=str,
        required=False,
        default="no",
        help='''
        str, options are yes or no, default is no.

        If specified as yes, will invert the strand of each transcript in
        the input BED file. This is useful because sometimes in strand-specific
        sequencing one of the read pairs is inverted relative to the RNA
        molecule.
        '''
    )
    parser.add_argument(
        '-p',
        '--paired_end_input',
        type=str,
        required=False,
        default="yes",
        help='''
        str, options are yes or no, default is no.

        This influences how multiple lines in a bed file with the same read name
        are handled.

        If this is set to "yes", IT IS ASSUMED THAT THERE ARE NO SUPPLEMENTARY
        OR SECONDARY ALIGNMENTS IN THE BED FILE. With this assumption, if a
        second instance of a read name is encountered it will be assumed
        the bed file is formed from mapping of paired-end reads. Therefore, when
        a read name is encountered for a second time the second read name will
        be converted to read_name_2.

        If this is set to "no", the assumption is that there are NOT paired end
        reads in the input. If this is the case, it's assumed that the second
        instance of the read name is a supplementary alignment. In this case,
        the primary and secondary alignments are both tossed.
        '''
    )
    parser.add_argument(
        '-j',
        '--junc_reads_only',
        type=str,
        required=False,
        default="no",
        help='''
        str, options are yes or no, default is no.

        If this is set to yes, only reads that have junctions will be output
        to the span file. This is helpful to reduce file sizes when only
        concerned about junctions.
        '''
    )

    args = parser.parse_args()

    # Validate args
    if args.invert_strand not in ["yes", "no"]:
        msg = "invert_strand must be yes or no. You input {}.".format(args.invert_strand)
        raise ValueError(msg)
    if args.paired_end_input not in ["yes", "no"]:
        msg = "paired_end_input must be yes or no. You input {}.".format(args.paired_end_input)
        raise ValueError(msg)
    if args.junc_reads_only not in ["yes", "no"]:
        msg = "junc_reads_only must be yes or no. You input {}.".format(args.junc_reads_only)
        raise ValueError(msg)

    return args

#------------------------------------------------------------------------------#
# Constants
#------------------------------------------------------------------------------#
bed_cols = "chrom tx_start tx_end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts"
bed_cols = bed_cols.split(" ")

#------------------------------------------------------------------------------#
# Transcript object and associated functions
#------------------------------------------------------------------------------#
class tx_object:
    """
    Holds information for a transcript.
    """

    def __init__(self, name, tx_start, tx_end, strand, blockSizes, blockStarts):
        self.name = name
        self.tx_start = tx_start
        self.tx_end = tx_end
        self.strand = strand
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts

        # Solve for spans and juncs
        self.spans, self.juncs = find_spans(tx_start, tx_end, blockSizes, blockStarts, self.name)

def find_spans(tx_start, tx_end, blockSizes, blockStarts, tx_name):
    """
    This function parses several BED file fields and produces a spans and juncs output.

    Inputs:
    - tx_start: start of the transcript (e.g. bed file chromStart)
    - tx_end: end of the transcript (e.g. the bed file chromEnd)
    - blockSizes: blockSizes field from the bed file. May end with a free comma.
    - blockStarts: blockStarts field from the bed file. May end with a free comma.
    - tx_name: Name of the transcript. This is just used for the error messages.

    blockSizes and blockStarts should be the same length. These items describe each EXON, e.g.
    each exon starts at position mentioned in index i of blockStarts and has a lenght defined at
    index i of blockSizes.

    Output: tuple (spans, juncs)
    - spans: List of tuples that are alternating between exons and introns. This object shows the coordinates
             for each span. Ex:
             [(56, 194, 'exon'),
             (194, 334, 'intron'),
             (334, 426, 'exon')]

    - juncs: list of tuples that simply has the 5' and 3' position of each junction - e.g. only the coordinates
             from the intron spans. Ex:
             [(194, 334)]
    """

    # Make sure tx_start and tx_end are integers
    tx_start = int(tx_start)
    tx_end = int(tx_end)

    # Parse the blocks into lists
    blockSizes = blockSizes.rstrip(",").split(",")
    blockStarts = blockStarts.rstrip(",").split(",")

    # Sanity check
    if len(blockSizes) != len(blockStarts):
        msg = "Found instance where the number of blockSizes and number of blockStarts is not equal!"
        msg = msg + " Transcript is {}".format(tx_name)
        raise ValueError(msg)

    # Iterate through each block and fill out spans and junctions
    spans = []
    juncs = []
    for block in range(len(blockSizes)):

        # Find their relative start and length
        block_start = int(blockStarts[block])
        block_len = int(blockSizes[block])

        # Adjust start to tx_start
        block_start = block_start + tx_start

        # Find current exon
        current_exon = (block_start, block_len)

        # If this is not the first block, first load the preceeding intron
        if block > 0:
            preceeding_intron_start = spans[-1][1]
            preceeding_intron_end = block_start
            spans.append((preceeding_intron_start, preceeding_intron_end, "intron"))

            # Also load into juncs
            juncs.append((preceeding_intron_start, preceeding_intron_end))

        # Load the current exon
        spans.append((block_start, block_start + block_len, "exon"))


    ### Sanity checks
    # Make sure first and last spans are exons
    if spans[0][2] != "exon":
        msg = "First span is not an exon!"
        msg = msg + " Transcript is {}".format(tx_name)
        raise ValueError(msg)
    if spans[-1][2] != "exon":
        msg = "Last span is not an exon!"
        msg = msg + " Transcript is {}".format(tx_name)
        raise ValueError(msg)

    # Make sure last exon ends at tx_end
    if spans[-1][1] != tx_end:
        msg = "Last exon doesn't end at tx_end!"
        msg = msg + " Transcript is {}".format(tx_name)
        raise ValueError(msg)

    # Return the answer!
    return spans, juncs

#------------------------------------------------------------------------------#
# Functions for this script
#------------------------------------------------------------------------------#
def parse_bed_to_tx_objects(bed_path, invert_strand="no", paired_end_input="yes", junc_reads_only="no"):
    """
    Given a path to a bed file, will parse each line and load them into a transcript
    object (one per line). Thus, returns a dictionary of tx_name:tx_object.
    """

    tx_objects = dict()
    with open (bed_path) as infile:

        for line in infile:
            line = line.rstrip("\n").split("\t")

            # Parse fields
            name = line[bed_cols.index("name")]
            tx_start = line[bed_cols.index("tx_start")]
            tx_end = line[bed_cols.index("tx_end")]
            strand = line[bed_cols.index("strand")]
            blockCount = line[bed_cols.index("blockCount")]
            blockSizes = line[bed_cols.index("blockSizes")]
            blockStarts = line[bed_cols.index("blockStarts")]

            # Skip non-junction reads if specified
            if blockCount < 2 and junc_reads_only == "yes":
                continue

            # Invert strand if specified
            if invert_strand == "yes":
                if strand == "+":
                    strand = "-"
                elif strand == "-":
                    strand = "+"

            # Load tx_object
            tx = tx_object(name, tx_start, tx_end, strand, blockSizes, blockStarts)

            # Handle read names that occur more than once.
            # If not paired end, assume this is a secondary alignment and ditch
            # the reads entirely.
            if paired_end_input == "no" and name in tx_objects:
                print("{} has a supplementarty alginment, so not keeping it..".format(name))
                del tx_objects[name]
                continue

            # Otherwise assume it is a second read of a read pair
            if paired_end_input == "yes" and name in tx_objects:
                name += "_2"
                tx.name = name

            tx_objects[name] = tx

    return tx_objects

def count_and_rank_juncs(tx_objects):
    """
    Takes in a dictionary of tx_objects, and returns a dictionary of structure
    juncs: (rank, count).

    Each key is simply the junction of a tx_object converted to a string.
    """

    tx_objects = tx_objects.copy()

    # Count and rank junction classes
    junc_counter = Counter()
    for tx in tx_objects.values():

        junc_counter[str(tx.juncs) + tx.strand] += 1

    # Convert counter to a dictionary of structure - juncs: (rank, count)

    junc_ranked = {pair[0]: (rank, pair[1])
    for rank, pair in enumerate(junc_counter.most_common())}

    return junc_ranked

def load_rank_and_count(tx_objects, junction_ranks):
    """
    Loads the rank (aka tx_class) and count of each tx_class to
    each transcript object.
    """
    tx_objects = tx_objects.copy()
    for name, tx in tx_objects.items():

        rank, count = junction_ranks[str(tx.juncs) + tx.strand]

        tx.tx_class = rank
        tx.tx_class_count = count

        tx_objects[name] = tx

    return tx_objects

def write_tx_spans(in_tx_dict, outfile_path):

    # Make output directory if needed
    out_dir = os.path.dirname(outfile_path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Initiate with colheaders
    output = "name start end strand span_type tx_start tx_end tx_class tx_class_count\n".replace(" ", "\t")

    with open(outfile_path, "w") as outfile_handle:


        for name, tx in in_tx_dict.items():

            for span in tx.spans:

                start, end, span_type = span

                entry = "{} {} {} {} {} {} {} {} {}\n".format(
                    tx.name, start, end, tx.strand, span_type,
                    tx.tx_start, tx.tx_end, tx.tx_class,
                    tx.tx_class_count
                ).replace(" ", "\t")

                output += entry

        outfile_handle.write(output)


def main():

    # Sort out the arguments
    #--------------------------------------------------------------------#
    args = get_args()
    bed_path = args.bed_path
    output_spans = args.output_spans
    invert_strand = args.invert_strand
    paired_end_input = args.paired_end_input
    junc_reads_only = args.junc_reads_only

    # Main
    #--------------------------------------------------------------------#
    print("{}: Starting script".format(sys.argv[0]))

    # Parse the bed file
    tx_objects = parse_bed_to_tx_objects(bed_path, invert_strand, paired_end_input, junc_reads_only)

    # Count and rank the junctions... make new dict of structure
    # junc: (rank, count)
    junction_ranks = count_and_rank_juncs(tx_objects)

    # Load rank (I'll call it tx_class) and count (tx_class_count) to each tx
    tx_objects = load_rank_and_count(tx_objects, junction_ranks)

    # Write to output file.
    write_tx_spans(tx_objects, output_spans)

    print("{}: Finished".format(sys.argv[0]))

if __name__ == '__main__':
    main()
