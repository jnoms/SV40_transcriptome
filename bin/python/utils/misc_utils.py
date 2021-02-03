#!/usr/bin/env python3

# CLASS AND FUCTION FOR TRANSCRIPT ANALYSIS
# THIS SCRIPT IS NOT MEANT TO BE RUN
# DIRECTLY.

# REQUIRED PACKAGES FOR THESE UTILS:
import pathlib
import os
import sys
import gzip
import io
import pathlib
import numpy as np

#------------------------------------------------------------------------------#
# Constants
#------------------------------------------------------------------------------#

def open_file(path, read_or_write="r"):
    """
    Wrapper for opening files in read/write, but handles if the path ends in .gz.
    Function also will make output directories if needed.

    - if path ends with .gz, will open with the gzip package. Otherwise, will use
    standard open function.
    - if read_or_write == r, will open in read mode. If it == w, will open in write mode.
      I think if you set it to 'a' it'll append.
    """

    # Validate the input parameter
    if read_or_write not in ["r", "w", "a"]:
        msg = "read_or_write must be 'r' or 'w'. You entered {}".format(read_or_write)
        raise ValueError(msg)

    # Make sure output file dir is made if using function in write or append mode
    if read_or_write != "r":
        out_dir = os.path.dirname(path)
        pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Open with correct package
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.GzipFile(path, read_or_write))
    else:
        return open(path, read_or_write)

def parse_bed_blocks(input_bed):
    """
    Takes in the path to a bed file, and returns a dictionary
    containing information on the end of each exon block
    for each transcript (relative to the transcript). An exon block
    end here is simply the cumulative sum of the blockSizes, with the
    exception of the last cumsum blockSize which is the end of the transcript.

    Output dictionary is of structure
    tx_name:[pos1, pos2, ...] where each pos is the position of a splice.

    Note - bed file is assumed to be 0-indexed, prodigal 1-indexed.
    """

    bed_splice_sites = dict()
    with open_file(input_bed) as infile:

        for line in infile:

            line = line.rstrip("\n").split("\t")

            # Parse bed
            blockSizes = line[10].rstrip(",").split(",")
            blockSizes = [int(n) for n in blockSizes]
            tx_name = line[3]
            n_blocks = int(line[9])

            # If there is only one block there are no junctions, so continue
            if n_blocks == 0:
                continue

            # The end of each block is the cumulative sum.
            # Exception is the last cumsum, which is the end of the transcript
            # rather than a splice site.
            splice_sites = list(np.cumsum(blockSizes))[:-1]

            # Handle repetative entries (which shouldn't occur...)
            if tx_name in bed_splice_sites:
                msg = "Found a duplicate transcript name in the bed. "
                msg += "Make sure you disabled supplemental alignments "
                msg += "and that, if using paired end illumina reads, the "
                msg += "reads are labeled with _1 or _2 prior to mapping. "
                msg += "Problematic transcript is {}".format(tx_name)
                raise ValueError(msg)

            # Write to dictionary
            bed_splice_sites[tx_name] = splice_sites

    return bed_splice_sites

def span_file_to_tx_class_dict(in_spans, span_columns):
    """
    Input:
        - path to span file

    Output:
        - dict of structure tx_name:(tx_class, tx_class_count)
    """

    tx_class_dict = dict()
    with open_file(in_spans) as infile:

        for line in infile:

            # Skip header and newlines
            if line.startswith(span_columns[0]):
                continue
            if line == "\n":
                continue

            # Parse
            line = line.rstrip("\n").split("\t")
            tx_name = line[span_columns.index("name")]
            tx_class = line[span_columns.index("tx_class")]
            tx_class_count = line[span_columns.index("tx_class_count")]

            # Write to output dictionary
            tx_class_dict[tx_name] = (tx_class, tx_class_count)

    return tx_class_dict


if __name__ == '__main__':
    msg = "This script contains functions and classes, and is not meant to be \
    run directly."
    raise ValueError(msg)
