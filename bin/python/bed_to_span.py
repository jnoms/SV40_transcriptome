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
        '-i',
        '--illumina_spans',
        type=str,
        required=False,
        default="",
        help='''
        Path to the spans file, created by this script, from illumina reads.
        If this is specified, this script will consider junctions in the input BED
        by illumia confidence.
        '''
    )
    parser.add_argument(
        '-g',
        '--n_genomes',
        type=int,
        required=False,
        default=20,
        help='''
        Only required if incorporating illumina junctions support.

        This should specify the number of genomes the input bed file was mapped
        against. Default is set to 20, so as long as the bed is from a reference
        that is 20 genomes or less you are all set.
        '''
    )
    parser.add_argument(
        '-u',
        '--n_junctions',
        type=int,
        required=False,
        default=5,
        help='''
        Only required if incorporating illumina junctions support.
        Default: 5

        This is an interger value detailing the minimum number of
        illumina reads with all junctions within a transcript of
        the input BED to designate the BED as high-confidence.
        '''
    )
    parser.add_argument(
        '-l',
        '--genome_length',
        type=int,
        required=False,
        default=0,
        help='''
        Only required if incorporating illumina junctions support.

        This is the genome length of the viral reference genome. Used to propogate
        illumina junctions forward to n_genomes.

        SV40: 5243
        BK: 5153
        '''
    )
    parser.add_argument(
        '-x',
        '--output_spans_unsupported',
        type=str,
        required=False,
        default="",
        help='''
        Only required if incorporating illumina junctions support.

        Path to the output file containing the spans from the unsupported transcripts.
        These transcripts contain at least one junction with fewer illumina reads than
        n_junctions supporting them.
        '''
    )

    parser.add_argument(
        '-d',
        '--keep_potential_duplicates',
        type=str,
        required=False,
        default="no",
        help='''
        Options: yes, no. Default: no

        If set to no (the default), if the script encounters a read with a read
        name that has been seem before, it considers it to be a non-primary
        alignment (secondary or supplementary) and discards both alignments.

        If set to yes, it does NOT care if it encounters the same read name
        more than once. This is helpful if the input reads were from e.g. paired
        end sequencing and the pairs weren't labeled with e.g. _1 and _2. When
        you use this option, make sure that the alignment did NOT allow non-
        primary alignments else they will be represented more than once in the
        span output.
        '''
    )

    args = parser.parse_args()

    # Validate args
    if args.junc_reads_only not in ["yes", "no"]:
        msg = "junc_reads_only must be yes or no. You input {}.".format(args.junc_reads_only)
        raise ValueError(msg)

    return args

#------------------------------------------------------------------------------#
# Constants
#------------------------------------------------------------------------------#
bed_cols = "chrom tx_start tx_end name score strand thickStart thickEnd itemRgb blockCount blockSizes blockStarts"
bed_cols = bed_cols.split(" ")
span_columns = "name start end strand span_type tx_start tx_end tx_class tx_class_count tx_class_illumina_support".split(" ")


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

    blockSizes and blockStarts should be the same length. These items describe
    each EXON, e.g. each exon starts at position mentioned in index i of
    blockStarts and has a lenght defined at index i of blockSizes.

    Output: tuple (spans, juncs)
    - spans: List of tuples that are alternating between exons and introns. This
      object shows the coordinates
             for each span. Ex:
             [(56, 194, 'exon'),
             (194, 334, 'intron'),
             (334, 426, 'exon')]

    - juncs: list of tuples that simply has the 5' and 3' position of each
      junction - e.g. only the coordinates
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
def parse_bed_to_tx_objects(bed_path, junc_reads_only="no", keep_potential_duplicates="no"):
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
            # Assume this is a secondary alignment and ditch the reads entirely.
            # This can be disabled with keep_potential_duplicates != "no"
            if name in tx_objects and keep_potential_duplicates == "no":
                print("{} has a supplementarty alginment, so not keeping it..".format(name))
                del tx_objects[name]
                continue

            tx_objects[name] = tx

    return tx_objects

def illumina_spans_to_junction_counts(illumina_spans, n_genomes, genome_length, span_columns):
    """
    This function reads in and counts the junctions from illumian spans. It first "slides"
    all junctions/introns so they start in the first genome. It then counts instances of
    each junction. Finally, it copies the same count for each position in all other subsequent
    n_genomes (e.g. start + genome_length*i for i in range(n_genomes) and same for end)

    Input:
        - illumina_spans: path to the illumian spans file
        - n_genomes: Number of genomes to propogate forward
        - genome_length: Length of a single viral reference.
        - span_columns: colnames of the span file.

    Output:
        - propogated_junction_counts: Counter of structure
        (start, end, strand): count

    """

    # If no illumina_spans provided, just return an empty counter
    if illumina_spans == "":
        return Counter()

    # (start, end, strand):count
    junction_counts = Counter()

    with open(illumina_spans) as infile:

        for line in infile:

            if line.startswith("name"): continue

            line = line.split()

            # Keep only introns
            if line[span_columns.index("span_type")] == "exon":
                continue

            # Count start, stop, and strand
            start = int(line[span_columns.index("start")])
            end = int(line[span_columns.index("end")])
            strand = line[span_columns.index("strand")]

            # Slide junctions to consolidate in the first genome - will propogate
            # forward in the next loop
            while start > genome_length:
                start = start - genome_length
                end = end - genome_length

            # Add to count
            junction_counts[(start, end, strand)] += 1

    # Propogate to n genomes
    propogated_junction_counts = Counter()
    for (start, end, strand), count in junction_counts.items():

        for i in range(n_genomes):
            new_start = start + (genome_length * i)
            new_end = end + (genome_length * i)

            propogated_junction_counts[(new_start, new_end, strand)] = count

    return propogated_junction_counts


def split_tx_by_support(tx_objects, illumina_junction_counts, n_junctions):
    """
    This function takes in the tx_objects dictionary and separates transcripts
    that are supported by > n_junctions illumina junctions and those that are
    not. It returns two dictionary containing the reads in question.

    If no illumina spans were provided, illumina_junction_counts will be an
    empty Counter(). All tx_objects will then immedaitely be added to the
    supported_tx_objects output dictionary.

    Input:
        - tx_objects: Dict (tx_name:tx) for each transcript. tx is a transcript
        object with all of the transcript information.
        - illumina_junctions_counts: Counter() of structure (start, end, strand):count
        that details number of supporting illumina reads for each junction.
        - n_junctions: The minimum number of illumina reads to support ALL JUNCTIONS
        in a tx for the tx to be considered supported. AGAIN, if a single junction
        within a multi-junction tx is not-supported, the entire transcript is labled
        not-supported.

    Output:
        - tuple of (supported_tx_objects, unsupported_tx_objects)
        - Each is a dictionary, containing supported vs unsupported transcripts.
        - If the input illumina_junction_counts is an empty counter, all transcripts
        will automatically be in supported_tx_objects

    """

    supported_tx_objects = dict()
    unsupported_tx_objects = dict()

    for tx_name, tx in tx_objects.items():

        # If there are no illumina counts, spand weren't provied. Just label
        # the tx accordingly and put them as supported.
        if illumina_junction_counts == Counter():
            tx.illumina_support = "NA"
            supported_tx_objects[tx_name] = tx
            continue

        # Determine if any junctions don't have > n illumina reads
        illumina_support = "supported"
        for junction in tx.juncs:
            lookup = (junction[0], junction[1], tx.strand)

            illumina_count = illumina_junction_counts[lookup]

            if illumina_count < n_junctions:
                illumina_support = "unsupported"
                break

        # Label tx will it's status
        tx.illumina_support = illumina_support

        # Split to new dictionary accordingly
        if illumina_support == "supported":
            supported_tx_objects[tx_name] = tx
        else:
            unsupported_tx_objects[tx_name] = tx


    print("There are {} supported and {} unsupported transcripts.".format(
        len(supported_tx_objects), len(unsupported_tx_objects))
         )

    return supported_tx_objects, unsupported_tx_objects


def count_and_rank_juncs(tx_objects):
    """
    Takes in a dictionary of tx_objects, and returns a dictionary of structure
    juncs: (rank, count).

    Each key is simply the junction of a tx_object converted to a string.
    """

    tx_objects = tx_objects.copy()

    # Count and rank junction classes
    junc_counter = Counter()
    unspliced_counter = Counter()
    for tx in tx_objects.values():

        # If unspliced, reserve so it doesn't get ranked
        if tx.juncs == []:
            unspliced_counter[str(tx.juncs) + tx.strand] += 1
            continue

        junc_counter[str(tx.juncs) + tx.strand] += 1

    # Convert counter to a dictionary of structure - juncs: (rank, count)
    junc_ranked = {pair[0]: (rank + 1, pair[1])
    for rank, pair in enumerate(junc_counter.most_common())}

    # Add the unspliced to the junc_ranked
    for unspliced, count in unspliced_counter.items():
        junc_ranked[unspliced] = ("u" + unspliced[-1], count)

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

def write_tx_spans(in_tx_dict, outfile_path, span_columns):

    # Make output directory if needed
    out_dir = os.path.dirname(outfile_path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Initiate with colheaders
    output = "\t".join(span_columns) + "\n"

    with open(outfile_path, "w") as outfile_handle:


        for name, tx in in_tx_dict.items():

            for span in tx.spans:

                start, end, span_type = span

                entry = "{} {} {} {} {} {} {} {} {} {}\n".format(
                    tx.name, start, end, tx.strand, span_type,
                    tx.tx_start, tx.tx_end, tx.tx_class,
                    tx.tx_class_count, tx.illumina_support
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
    illumina_spans = args.illumina_spans
    n_genomes = args.n_genomes
    n_junctions = args.n_junctions
    genome_length = args.genome_length
    keep_potential_duplicates = args.keep_potential_duplicates
    output_spans_unsupported = args.output_spans_unsupported

    # Main
    #--------------------------------------------------------------------#
    print("{}: Starting script".format(sys.argv[0]))

    # Parse the bed file
    tx_objects = parse_bed_to_tx_objects(bed_path, junc_reads_only, keep_potential_duplicates)

    # Extract junction counts from the illumian spans file
    illumina_junction_counts = illumina_spans_to_junction_counts(illumina_spans, n_genomes, genome_length, span_columns)

    # Split tx_objects into supported/not-support transcripts. If no illumina_spans were provided,
    # just put all into the "supported" one.
    supported_tx_objects, unsupported_tx_objects = split_tx_by_support(tx_objects, illumina_junction_counts, n_junctions)

    # Count and rank the junctions... make new dict of structure
    # junc: (rank, count)
    supported_junction_ranks = count_and_rank_juncs(supported_tx_objects)
    unsupported_junction_ranks = count_and_rank_juncs(unsupported_tx_objects)

    # Load rank (I'll call it tx_class) and count (tx_class_count) to each tx
    supported_tx_objects = load_rank_and_count(supported_tx_objects, supported_junction_ranks)
    unsupported_tx_objects = load_rank_and_count(unsupported_tx_objects, unsupported_junction_ranks)

    # Write to output file.
    write_tx_spans(supported_tx_objects, output_spans, span_columns)
    if output_spans_unsupported != "":
        write_tx_spans(unsupported_tx_objects, output_spans_unsupported, span_columns)

    print("{}: Finished".format(sys.argv[0]))

if __name__ == '__main__':
    main()
