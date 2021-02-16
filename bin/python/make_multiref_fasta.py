#!/usr/bin/env python3

import pathlib
import os
import argparse
import gzip
import sys
from Bio import SeqIO

from utils.misc_utils import open_file

def get_args():
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to concatenate the sequence of an
            input fasta n times (keeping the one input header as the header).
            """)

    # Required arguments
    parser.add_argument(
        '-i',
        '--in_fasta',
        type=str,
        required=True,
        help='''
        Path to the input fasta. There should be only one entry/header.
        '''
    )
    parser.add_argument(
        '-o',
        '--out_fasta',
        type=str,
        required=True,
        help='''
        Path to the output fasta, which will contain the sequence n times.
        '''
    )
    parser.add_argument(
        '-n',
        '--n',
        type=int,
        required=True,
        help='''
        Number of times to repeat the sequence.
        '''
    )

    args = parser.parse_args()
    return args

def read_fasta_to_memory(input_fasta):
    """
    Reads fasta into a memory as a dictionary with header:sequence.
    This function can handle .gzip files, but input_fasta needs to
    end with .gz
    """
    fasta_dict = dict()

    if not input_fasta.endswith(".gz"):
        for seq_record in SeqIO.parse(input_fasta, "fasta"):
            fasta_dict[seq_record.id] = str(seq_record.seq)

    elif input_fasta.endswith(".gz"):
        with gzip.open(input_fasta, "rt") as handle:
            for seq_record in SeqIO.parse(handle, "fasta"):
                fasta_dict[seq_record.id] = str(seq_record.seq)

    return fasta_dict




def main():

    # Sort out the arguments
    #--------------------------------------------------------------------#
    args = get_args()
    in_fasta = args.in_fasta
    out_fasta = args.out_fasta
    n = args.n

    # Main
    #--------------------------------------------------------------------#
    print("{}: Starting script".format(sys.argv[0]))

    fasta = read_fasta_to_memory(in_fasta)

    if len(fasta.keys()) > 1:
        msg = "This script is designed for a single fasta entry as input."
        msg += " Detected multiple headers in the input fasta."
        raise ValueError(msg)

    with open_file(out_fasta, "w") as outfile:
        out = ""
        for header, seq in fasta.items():
            out_seq = seq * n
            out += ">" + header + "\n" + out_seq
        outfile.write(out)


    print("{}: Finished".format(sys.argv[0]))

if __name__ == '__main__':
    main()
