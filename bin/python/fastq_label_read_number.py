#!/usr/bin/env python3

import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pathlib
import os
import io
from Bio.Seq import Seq
import argparse
import sys
import random

from utils.misc_utils import open_file

def get_args():
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to add _1 or _2 to each read name in
            an input fastq file. Must specify if the input file is read1 or
            read2.

            NOTE - the lable will be appended to the first string if the read
            title is a space-delimited entry of strings.
            """)

    # Required arguments
    parser.add_argument(
        '-i',
        '--infile',
        type=str,
        required=True,
        help='''
        Path to the input fastq file. If ends in .gz will be treated as a
        gzipped file, otherwise will be treated as plaintext.
        '''
    )
    parser.add_argument(
        '-o',
        '--outfile',
        type=str,
        required=True,
        help='''
        Path to the output fastq file that will contain labed read names.
        If ends in .gz the output will be gzipped. Otherwise, the
        output will be plaintext.
        '''
    )
    parser.add_argument(
        '-n',
        '--read_number',
        type=int,
        required=True,
        help='''
        Will add _n to the end of each read name.
        '''
    )

    parser.add_argument(
        '-r',
        '--add_random',
        type=str,
        required=False,
        default="no",
        help='''
        Default is no. If set to anything but no, adds a random number to the
        end of each read name. It's totally random per read.
        '''
    )

    args = parser.parse_args()

    return args


def main():

    # Sort out the arguments
    #--------------------------------------------------------------------#
    args = get_args()
    infile = args.infile
    outfile = args.outfile
    read_number = args.read_number
    add_random = args.add_random

    # Main
    #--------------------------------------------------------------------#
    print("{}: Starting script".format(sys.argv[0]))

    # Make output file directory if needed
    out_dir = os.path.dirname(outfile)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    with open_file(infile) as infile_handle:

        with open_file(outfile, "w") as outfile_handle:

            for title, seq, qual in FastqGeneralIterator(infile_handle):

                # Label the read name
                title = title.split(" ")

                # either add n or random number
                if add_random == "no":
                    title[0] += "_" + str(read_number)
                else:
                    title[0] += "_" + str(random.randint(1,10000))
                title = " ".join(title)

                # Write it out
                outfile_handle.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))


    print("{}: Finished".format(sys.argv[0]))


if __name__ == '__main__':
    main()
