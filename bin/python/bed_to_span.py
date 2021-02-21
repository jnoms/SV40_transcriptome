#!/usr/bin/env python3

import argparse
import pathlib
import os
from collections import Counter
import sys
from statistics import median
from copy import deepcopy # lets you copy an object

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
    parser.add_argument(
        '-r',
        '--tx_representative_spans',
        type=str,
        required=False,
        default="",
        help='''
        Path to a spans file containing representatives of each tx_class. For
        each tx_class, the representative contains the defining junctions, but
        the tx_start and tx_end are the MEDIAN for the transcript class.

        Default: ""
        '''
    )
    parser.add_argument(
        '-s',
        '--tx_end_binsize',
        type=int,
        required=False,
        default=200,
        help='''
        The end position of each transcript is factored into placing a tx into
        a tx_class. However, dRNAseq reads may have some messiness to ending
        sites. The genome is split into bins of size tx_end_binsize - a tx
        must end within the same bin to be classified as the same tx_class.

        For Illumina reads where tx_end doesn't matter, or if you simply dont
        want to factor in tx_end, set this number to any large number bigger
        than the reference length.

        Default: 200
        '''
    )

    args = parser.parse_args()

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
def parse_bed_to_tx_objects(bed_path, junc_reads_only="no"):
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
            if int(blockCount) < 2 and junc_reads_only == "yes":
                continue

            # Load tx_object
            tx = tx_object(name, tx_start, tx_end, strand, blockSizes, blockStarts)

            # Handle read names that occur more than once.
            # Assume this is a secondary alignment and ditch
            # the reads entirely.
            if name in tx_objects:
                print("{} has a supplementarty alginment, so not keeping it..".format(name))
                del tx_objects[name]
                continue

            tx_objects[name] = tx

    return tx_objects

def get_tx_end_bin(tx_objects, tx_end_binsize = 100):
    """
    Input:
        - dictionary of structure tx_name:tx_object

    Output:
        - same dictionary, but each tx_object has an added .tx_end_bin attribute.

    The tx_end_bin is determined as the numeric bin that the tx_end falls within. For
    example, if the tx_end_binsize = 100 and the transcript ends at 450, the bin would be 4.
    If the transcript ends at 14, the bin would be 0. Bin is determined as
    int(tx_end/tx_end_binsize) (int will round down to nearest whole number).

    Note that this function modifies the input dictionary in place and doesn't return anything.
    """
    for tx in tx_objects.values():
        tx_end_bin = int(int(tx.tx_end)/tx_end_binsize)
        tx.tx_end_bin = tx_end_bin

def count_and_rank_juncs(tx_objects):
    """
    Takes in a dictionary of tx_objects, and returns a dictionary of structure
    juncs: (rank, count).

    Each key tx.juncs + tx.stand + tx.tx_end_bin, meaning that it takes into account
    both the junctions, the strand, and the end bin.
    """

    # Count and rank junction classes
    junc_counter = Counter()
    unspliced_counter = Counter()
    for tx in tx_objects.values():

        # If unspliced, reserve so it doesn't get ranked
        if tx.juncs == []:
            unspliced_counter["u" + tx.strand] += 1
            continue

        junc_counter[str(tx.juncs) + tx.strand + str(tx.tx_end_bin)] += 1

    # Convert counter to a dictionary of structure - juncs: (rank, count)
    junc_ranked = {pair[0]: (rank + 1, pair[1])
    for rank, pair in enumerate(junc_counter.most_common())}

    # Add the unspliced to the junc_ranked
    for unspliced_type, count in unspliced_counter.items():
        #unspliced_type is u + strand (e.g. u+ or u-)
        junc_ranked[unspliced_type] = (unspliced_type, count)

    return junc_ranked

def load_rank_and_count(tx_objects, junction_ranks):
    """
    Loads the rank (aka tx_class) and count of each tx_class to
    each transcript object.

    Note that this modifies the input dictionary in place and doesn't
    return anything.
    """
    for name, tx in tx_objects.items():

        # If the tx isn't unspliced:
        if tx.juncs != []:
            rank, count = junction_ranks[str(tx.juncs) + tx.strand + str(tx.tx_end_bin)]
            tx.tx_class = rank
            tx.tx_class_count = count
            tx_objects[name] = tx
        # No juncs, so just look up count based on strand
        else:
            rank, count = junction_ranks["u" + tx.strand]
            tx.tx_class = rank
            tx.tx_class_count = count
            tx_objects[name] = tx

def alter_spans_begining_end(spans, new, begining=True):
    """
    Inputs:
        - spans: list of tuples, where each tupe is (start, end, span_type)
        - new: New start of the first span or end of the last span
        - begining: if True, modifies start of first span. If False,
        modifies to the end of the last span.

    Output:
        - modified spans

    Note that this copies the input and outputs a new list
    """
    spans = spans.copy()
    if begining == True:
        start, end, span_type = spans[0]
        spans[0] = (new, end, span_type)
        return spans

    elif begining == False:
        start, end, span_type = spans[-1]
        spans[-1] = (start, new, span_type)
        return spans

def get_tx_class_reps(tx_objects):
    """
    Input:
        - tx_objects: Dictionary of structure tx_name:tx_object

    Output:
        -tx_class_reps: One representative per tx_class. The representative
        will have the MEDIAN tx_start and MEDIAN tx_end of the tx_class, with
        the spans corrected accordingly.
    """

    # Find starts/ends for each tx_class, and save one tx as representative.
    # The representative will be altered to have the median start/end later.
    tx_class_reps = dict()
    tx_starts = dict()
    tx_ends = dict()
    for tx_name, tx in tx_objects.items():

        # Starts
        if not tx.tx_class in tx_starts:
            tx_starts[tx.tx_class] = []
        tx_starts[tx.tx_class].append(int(tx.tx_start))

        # Ends
        if not tx.tx_class in tx_ends:
            tx_ends[tx.tx_class] = []
        tx_ends[tx.tx_class].append(int(tx.tx_end))

        # Representative
        if not tx.tx_class in tx_class_reps:
            tx_class_reps[tx.tx_class] = deepcopy(tx) # avoid altering the original


    # Find median start/end for each tx_class
    tx_start_medians = dict()
    tx_end_medians = dict()
    for tx_class, starts in tx_starts.items():
        med = int(median(starts))
        tx_start_medians[tx_class] = med
    for tx_class, ends in tx_ends.items():
        med = int(median(ends))
        tx_end_medians[tx_class] = med

    # Adjust median for each representative
    output_tx = dict()
    for tx_class in tx_start_medians:
        med_start = tx_start_medians[tx_class]
        med_end = tx_end_medians[tx_class]
        tx_rep = tx_class_reps[tx_class]

        # Correct the start
        tx_rep.tx_start = med_start
        tx_rep.spans = alter_spans_begining_end(tx_rep.spans, med_start, begining=True)

        # Correct the end
        tx_rep.tx_end = med_end
        tx_rep.spans = alter_spans_begining_end(tx_rep.spans, med_end, begining=False)

        # Remove the now-innacurate other fields
        tx_rep.blockSizes = ""
        tx_rep.blockStarts = ""
        output_tx[tx_class] = tx_rep

    return output_tx

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
    junc_reads_only = args.junc_reads_only
    tx_representative_spans = args.tx_representative_spans
    tx_end_binsize = args.tx_end_binsize

    # Main
    #--------------------------------------------------------------------#
    print("{}: Starting script".format(sys.argv[0]))

    # Parse the bed file
    tx_objects = parse_bed_to_tx_objects(bed_path, junc_reads_only)

    # Find tx_end bin for each transcript
    get_tx_end_bin(tx_objects, tx_end_binsize)

    # Count and rank the junctions... make new dict of structure
    # junc: (rank, count)
    junction_ranks = count_and_rank_juncs(tx_objects)

    # Load rank (I'll call it tx_class) and count (tx_class_count) to each tx
    load_rank_and_count(tx_objects, junction_ranks)
    write_tx_spans(tx_objects, output_spans)

    # Find representative transcripts if specified
    if tx_representative_spans != "":
        tx_reps = get_tx_class_reps(tx_objects)
        write_tx_spans(tx_reps, tx_representative_spans)

    print("{}: Finished".format(sys.argv[0]))

if __name__ == '__main__':
    main()
