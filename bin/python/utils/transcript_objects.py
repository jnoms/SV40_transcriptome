#!/usr/bin/env python3

# CLASS AND FUCTION FOR TRANSCRIPT ANALYSIS
# THIS SCRIPT IS NOT MEANT TO BE RUN
# DIRECTLY.

# REQUIRED PACKAGES FOR THESE UTILS:
import pysam
import argparse
import pathlib
import os
import sys
import pandas as pd
from io import StringIO
from Bio.Seq import Seq

# to flatten complex list/tuple/etc objects
from pandas.core.common import flatten

#--------------------------------#
# Transcript utils objects
#--------------------------------#
class transcript_object:
    """
    This class is used to contain information about a transcript, and can be
    used to identify various characteristics. Note that start, end, and junction
    positions are all 0-indexed.

    Initializing information:
    - transcript_ID: Name of the transcript/read
    - start: Start position relative to the genome_seq
    - junctions: A list of tuples of structure. Each tuple dictates the
      5' and 3' positions of each junction, e.g. (5' position, 3' position)
    - end: The end of the transcript alignment relative to the genome_seq
    - is_supplementary: Boolean. Dictates if the transcript object is the primary or
      supplementary alignment.
    """

    def __init__(self, transcript_ID, start, junctions, end, strand, genome_seq, is_supplementary):

        # Set initial values
        self.transcript_ID = transcript_ID
        self.start = start
        self.junctions = junctions
        self.end = end
        self.strand = strand
        self.genome_seq = genome_seq
        self.is_supplementary = is_supplementary

    # Methods
    def print_info(self):
        """
        Simply prints out the information in the transcript object.
        """
        msg = """
        transcript_ID: {}
        start: {}
        junctions: {}
        end: {}
        strand: {}
        is_supplementary: {}
        """.format(self.transcript_ID,
                   self.start,
                   self.junctions,
                   self.end,
                   self.strand,
                   self.is_supplementary)
        print(msg)


    def find_seq(self):
        seq = ""
        for pos1, pos2, region_type in self.get_region_list():

            if region_type == "intron": continue

            if pos1 < pos2:
                region_sequence = self.genome_seq[pos1:pos2]
            else:
                region_sequence = self.genome_seq[pos2:pos1]
                region_sequence = str(Seq(region_sequence).reverse_complement())

            seq += region_sequence

        return seq

    def export_fasta_seq(self):
        """
        Returns fasta-formated string containing the corrected sequence.
        """
        return ">{}\n{}\n".format(self.transcript_ID, self.find_seq())

    def get_region_list(self):
        """
        This function processes sthe start, junctions, and end coordinates
        to generate a list of spans of exonic regions and intronic regions.
        As a reminder, coordinates are 0-indexed (similar to standard python
        indexing conventions).

        Output is a list of lists, always begining and ending with an exon.
        For example, given the input:
        - start: 142
        - junctions: [(294, 434), (526, 1462)]
        - end: 2666

        The expected output is:
        [[142, 294, 'exon'],
         [294, 434, 'intron'],
         [434, 526, 'exon'],
         [526, 1462, 'intron'],
         [1462, 2666, 'exon']]
        """

        # Define inputs
        start = self.start
        junctions = self.junctions
        end = self.end
        strand = self.strand

        coords = [start, junctions, end]
        region_list = []

        # Flatten this to a simple list
        coords = list(flatten(coords))

        # Sanity check - this should be even
        if len(coords)%2 != 0:
            msg = "Counting the start, end, and the 5' and 3' site of each junction "
            msg += "The flattened coordinates should be even. This message means it was odd!"
            raise ValueError(msg)

        # Iterate through each pair of numbers. Switch between exons and introns.
        # Zip generates pairwise tuples. Enumerate just gives them an index.
        for (i, (first, second)) in enumerate(zip(coords, coords[1:])):
            # Evens are exons, odds are introns
                if i%2 == 0:
                    region = [first, second, "exon"]
                    region_list.append(region)
                else:
                    region = [first, second, "intron"]
                    region_list.append(region)


        # Sanity checks - make sure that the first and last regions are exons
        if region_list[0][2] != "exon":
            msg = "The first region is not an exon. Something broke."
            raise ValueError(msg)
        if region_list[-1][2] != "exon":
            msg = "The last region is not an exon. Something broke."
            raise ValueError(msg)

        return region_list

    def export_region_report(self):
        """
        This function calculates the regions using the get_region_list()
        of this object. It then produces a tab-delimited output with the
        following fields:

        transcript_ID region_start region_end region_type strand transcript_start transcript_end

        This is a *tidy* output! So, each region will have it's own line.
        """
        regions = self.get_region_list()

        output_report = ""
        for region in regions:
            region_start = region[0]
            region_end = region[1]
            region_type = region[2]

            report_field = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(self.transcript_ID,
                                                         region_start,
                                                         region_end,
                                                         region_type,
                                                         self.strand,
                                                         self.start,
                                                         self.end)
            output_report += report_field

        return output_report


if __name__ == '__main__':
    msg = "This script contains functions and classes, and is not meant to be \
    run directly."
    raise ValueError(msg)
