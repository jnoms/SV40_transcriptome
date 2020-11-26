#!/usr/bin/env python3

# CLASS AND FUCTION FOR TRANSCRIPT ANALYSIS
# THIS SCRIPT IS NOT MEANT TO BE RUN
# DIRECTLY.

# REQUIRED PACKAGES FOR THESE UTILS:
import pathlib
import os
import sys

# to flatten complex list/tuple/etc objects
from pandas.core.common import flatten

def import_fasta(path):
    """
    Reads in a a fasta as a string. The fasta should have ONLY ONE
    entry/header.
    """

    header_count = 0
    seq = ""

    with open(path) as infile:
        for line in infile:
            if line.startswith(">"):
                header_count += 1
                continue

            seq += line.rstrip("\n")

    if header_count > 1:
        msg = "import_fasta: The fasta must contain only one entry/header.\n"
        msg = msg + "Detected {} headers.".format(header_count)
        raise ValueError(msg)

    return seq

def write_output(msg, path):

    # Make output directory if needed
    out_dir = os.path.dirname(path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    with open(path, "w") as outfile:
        outfile.write(msg)



if __name__ == '__main__':
    msg = "This script contains functions and classes, and is not meant to be \
    run directly."
    raise ValueError(msg)
