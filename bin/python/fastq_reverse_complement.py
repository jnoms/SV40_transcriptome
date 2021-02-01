#!/usr/bin/env python3

import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pathlib
import os
import io
from Bio.Seq import Seq
import argparse
import sys

from utils.misc_utils import open_file

def get_args():
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to reverse complement all of the
            sequences present in an input fastq. This file is capable of reading
            and writing standard or gzipped fastqs, and will do so depending
            on if .gz is included at the end of the input or output file strings
            entered to this script.
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
        Path to the output fastq file that will contain reverse-complemented
        sequences. If ends in .gz the output will be gzipped. Otherwise, the
        output will be plaintext.
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

    # Main
    #--------------------------------------------------------------------#
    print("{}: Starting script".format(sys.argv[0]))

    # Make output file directory if needed
    out_dir = os.path.dirname(outfile)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    with open_file(infile) as infile_handle:

        with open_file(outfile, "w") as outfile_handle:

            for title, seq, qual in FastqGeneralIterator(infile_handle):

                # Rev comp the seq
                seq = str(Seq(seq).reverse_complement())

                # Write it out
                outfile_handle.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))


    print("{}: Finished".format(sys.argv[0]))


if __name__ == '__main__':
    main()
