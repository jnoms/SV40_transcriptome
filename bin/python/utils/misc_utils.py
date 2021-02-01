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


if __name__ == '__main__':
    msg = "This script contains functions and classes, and is not meant to be \
    run directly."
    raise ValueError(msg)
