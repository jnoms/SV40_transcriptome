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


def write_output(msg, path, append_if_exists=False, append_if_exists_remove_header=False):
    """
    This function generates the output file, but has a bit extra complexity.

    If append_if_exists is set to True, it will append to the desired output file rather
    than generating a new one/overwriting any existing file.

    if append_if_exists_remove_header is set to True and we are appending,
    this function will remove the first line from msg prior to appending. This should be used
    if there is a redunant header in the msg that is no longer needed consdiering the extant
    file being appended to should have it.
    """

    # Make output directory if needed
    out_dir = os.path.dirname(path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    # If we're not appending no matter what, or the file doesn't exist, proceed accordingly
    if append_if_exists == False or not os.path.isfile(path):
        with open(path, "w") as outfile:
            outfile.write(msg)
            return

    # Otherwise, will be appending
    print("write_output: Because {} already exists, am appending to it.".format(path))

    # Now, check if we need to remove the header
    if append_if_exists_remove_header == True:
        msg = "\n".join(msg.split("\n")[1:])

    with open(path, "a") as outfile:
        outfile.write(msg)

def write_output_simple(msg, path):

    # Make output directory if needed
    out_dir = os.path.dirname(path)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    with open(path, "w") as outfile:
        outfile.write(msg)



if __name__ == '__main__':
    msg = "This script contains functions and classes, and is not meant to be \
    run directly."
    raise ValueError(msg)
