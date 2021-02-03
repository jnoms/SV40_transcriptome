#!/usr/bin/env python3

import pathlib
import os
import argparse
import gzip
import sys
import numpy as np
import io

from utils.misc_utils import open_file, parse_bed_blocks


def get_args():
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to filter an input prodigal file
            and only keep ORFs that overlap splice points of a transcript.
            The two inputs for this script are the prodigal file and the bed
            file that was initially used to generate the fasta that was
            processed with prodigal. This script infers junction positions
            in the derived fasta simply by identifying the end of each exon
            block in the BED file (using cumulative sums of the blockSizes).

            It is assumed that the bedtools getfasta was used on the bed file
            to genera the fasta which was processed with prodigal.

            By default this script will only allow ORFs on the + strand, but
            you can specify otherwise (although that wouldn't make much sense).
            """)

    # Required arguments
    parser.add_argument(
        '-p',
        '--input_prodigal',
        type=str,
        required=True,
        help='''
        Path to the input prodigal file. Gzipped can be handled, as can non-
        gzipped.
        '''
    )
    parser.add_argument(
        '-b',
        '--input_bed',
        type=str,
        required=True,
        help='''
        Path to the BED file used to generate the fasta which was processed
        with prodigal.
        '''
    )
    parser.add_argument(
        '-o',
        '--output_prodigal',
        type=str,
        required=True,
        help='''
        Path to the output prodigal file. If you specify the output file ending
        in .gz this script will gzip it automatically.
        '''
    )

    # Optional arguments
    parser.add_argument(
        '-s',
        '--allowable_strands',
        type=str,
        required=False,
        default="+",
        help='''
        Space-delimited list detailing the allowable strands. Only valid options
        are + and -. Default is + only.
        '''
    )

    args = parser.parse_args()
    args.allowable_strands = args.allowable_strands.split(" ")

    # Validate args
    for strand in args.allowable_strands:
        if not strand in ["+", "-"]:
            msg = "Possible space-delimited entries for allowable_strands are + or -. "
            msg += "You entered {}".format(strand)
            raise ValueError(msg)

    return args


def main():

    # Sort out the arguments
    #--------------------------------------------------------------------#
    args = get_args()
    input_prodigal = args.input_prodigal
    input_bed = args.input_bed
    output_prodigal = args.output_prodigal
    allowable_strands = args.allowable_strands

    # Main
    #--------------------------------------------------------------------#
    print("{}: Starting script".format(sys.argv[0]))

    # Make output file directory if needed
    out_dir = os.path.dirname(output_prodigal)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Parse bed to find blocks
    splice_sites = parse_bed_blocks(input_bed)

    # Filter and write out prodigal file
    with open_file(input_prodigal) as infile:

        with open_file(output_prodigal, "w") as outfile:

            # Track current transcript_name and splices
            current_transcript_name = ""
            current_splice_sites = ""

            # Iterate over each line
            for line in infile:

                # Skip empty lines
                if line == "\n":
                    continue

                # If it's a header line, find the header and then solve
                # for the current splice sites
                if line.startswith("# Sequence Data:"):

                    # Isolate the header
                    header = line.split("seqhdr=")[1].rstrip('\n').strip('"')

                    # Currently has some ORF information from bedtools... remove
                    header = header.split("::")[0]
                    header = "".join(header)

                    current_header = header
                    current_splice_sites = splice_sites[header]

                    # Write the header to the output and continue
                    outfile.write(line)
                    continue

                # If the line isn't blank and doesnt start with a # or "Beg" it is an ORF, so parse
                # it and determine if it overlaps a splice
                if not line.startswith("#") and line != "\n" and not line.startswith("Beg"):

                    line = line.split("\t")
                    start = int(line[0]) -1 # Adjust to 1-indexing
                    end = int(line[1])
                    strand = line[2]

                    # Determine if the strand is allowed
                    if not strand in allowable_strands:
                        continue

                    # Determine if it contains a splice
                    for splice_site in current_splice_sites:
                        if start <= splice_site <= end:
                            outfile.write("\t".join(line))
                            break

                # Still want to keep the headers and such, so write them out also
                else:
                    outfile.write(line)


    print("{}: Finished".format(sys.argv[0]))


if __name__ == '__main__':
    main()
