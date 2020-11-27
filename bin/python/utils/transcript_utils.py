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

from utils.transcript_objects import *
from utils.misc_utils import *

#--------------------------------#
# Transcript utils functions
#--------------------------------#
def convert_cigar(cigar, cigar_key):
    result = []

    for tup in cigar:
        symbol = cigar_key[tup[0]]
        result.append((symbol, tup[1]))

    return result

def get_end_coordinate(cigar, start):

    current_pos = start

    for tup in cigar:

        if tup[0] == "I":
            continue
        elif tup[0] in {"M", "N", "D", "X", "="}:
            current_pos += tup[1]

    return current_pos

def get_junctions(cigar, start, genome_seq, min_intron_length=100):
    """
    This function reads in a cigar sequence and outputs all junctions
    (maked with N (deletion) or N (intron)) that are at least the
    required length. The junction position is based on the start
    position, meaning it is syncronized with the genome.
    The output is a list of tuples, where each tuple is
    (junction_start, junction_end).
    The resultant junctions are 0-indexed. Therefore, if a junction is
    i.e. (11, 30) of a sequence that starts at 1 and ends at 55,
    the resultant sequence ignoring the junction is simply
    seq[1:11] + seq[30:55]
    """


    junctions = []
    pos = start

    for entry in cigar:

        sign, length = entry

        # ignore softclipping and insertions
        if sign in {"S", "I"}:
            continue

        # matches, = (equals), and X (mismatch) advance the position
        elif sign in {"M", "=", "X"}:
            pos += length
            continue

        # if deletion or intron, action depends if it passes
        # the length requirement
        elif sign in {"D", "N"}:

            # If under the length requirement, don't consider it a
            # junction but do advance position
            if length < min_intron_length:
                pos += length
                continue

            # Otherwise, we have a junction
            elif length >= min_intron_length:

                # Add junction to output
                junc = (pos, pos+length) # this is 0-indexed
                junctions.append(junc)

                # still need to advance the position for the next item(s)
                pos += length

    return junctions

def flip_antisense_transcript_positions(transcript_dict):
    """
    Takes ina dictionary of the structure:
    transcript_ID: [transcript_object1, transcript_object....]
    Usually there is only one transcript_object per list. If there are
    multiple, they're from supplemental or secondary alignments.

    Currently, the bam contains the start, stop and 5'/3' junction positions
    of each antisense transcript relative to the sense of the reference.
    This means the antisense transcript reported start is actually
    the end of the transcript and vice versa. The 5' and 3' junction sites
    are likewise switched. This function swaps them to be correct.
    """
    transcript_dict = transcript_dict.copy()

    for transcript_ID, transcript_objects in transcript_dict.items():
        for transcript_object in transcript_objects:

            if transcript_object.strand == "+":
                continue


            new_end = transcript_object.start
            new_start = transcript_object.end

            new_junctions = []
            for junction in transcript_object.junctions:
                new_three = junction[0]
                new_five = junction[1]
                new_junction = (new_five, new_three)
                new_junctions.append(new_junction)

            # Also need to reverse the order of the list for when generating the regions
            new_junctions.reverse()

            transcript_object.start = new_start
            transcript_object.end = new_end
            transcript_object.junctions = new_junctions

    return transcript_dict

def parse_bam(bam_path, cigar_key, genome_seq, min_intron_length=100, invert_strand=False, illumina_read_number=""):
    """
    The purpose of this function is to parse in a bam
    and generate a dictionary oftranscript_objects with
    alignment information.
    read_name = [transcript_object1, transcript_object2]

    Generally, there is only one transcript object. However, if there are
    *supplemental alignments* each discrete supplemental alignment will be
    added as an independent transcript object. These objects can be merged
    later.

    Generally this occurs rarely and is due to wraparound of the transcript
    on the circular genome.
    """

    transcript_dict = dict()

    # parse bam
    for read in pysam.AlignmentFile(bam_path, 'rb'):

        if read.is_unmapped:
            continue

        # convert the cigar tuple list to their single-letter
        # designations.
        cigar = convert_cigar(read.cigar, cigar_key)

        # Get start and end coordinates
        start = int(read.reference_start)
        end = int(get_end_coordinate(cigar, start))

        # Find junctions from cigar and start position
        junctions = get_junctions(cigar, start, min_intron_length)

        # Determine strand
        if read.is_reverse == False:
            strand = "+"
        else:
            strand = "-"

        # Inverse strand if called for
        if invert_strand == True:
            if strand == "+":
                strand = "-"
            elif strand == "-":
                strand = "+"

        # Determine if secondary or primary
        if read.is_supplementary == True or read.is_secondary == True:
            is_supplementary = True
        else:
            is_supplementary = False

        # Name of transcript
        transcript_ID = read.query_name

        # Label with read number if specified
        if illumina_read_number != "":
            transcript_ID += "_{}".format(illumina_read_number)

        # Generate the transcript object
        transcript = transcript_object(transcript_ID,
                                       start,
                                       junctions,
                                       end,
                                       strand,
                                       genome_seq,
                                       is_supplementary)

        # Add information from this read to the output
        if transcript_ID in transcript_dict:
            transcript_dict[transcript_ID].append(transcript)
        else:
            transcript_dict[transcript_ID] = [transcript]


    # Format antisense transcripts - need to swap 5' and 3' junction positions,
    # as well as start and stop.
    transcript_dict = flip_antisense_transcript_positions(transcript_dict)

    return transcript_dict


def limit_homology_ends(transcript_dict, min_end_homology):
    """
    Sometimes, minimap calls a big junction that ends at
    a region of homology that is really small. In these cases,
    I believe it is best to end the alignment instead at the end
    of the previous region of good homology.

    This function eliminates junctions that end within
    min_end_homology length of the end of the alignment.

    E.g. if min_end_homology is set to 50, if the
    last junction ends within 50 bases of the end of the alignment
    the last junction will be deleted, and teh end will be re-set
    to the previous start of the deleted junction.
    """
    transcript_dict = transcript_dict.copy()

    for transcript_ID, transcript_object_list in transcript_dict.items():

        for transcript_object in transcript_object_list:

            # Determine if length between end of final junction is
            # close to the end of the sequence
            junctions = transcript_object.junctions
            end = transcript_object.end

            if junctions == []:
                continue

            final_junction = junctions[-1]
            final_junction_start = final_junction[0]
            final_junction_end = final_junction[1]

            # Only continue for those alignments that have below
            # threshold
            if abs(end - final_junction_end) > min_end_homology:
                continue

            # For these, set the new end to the start of the last junction
            # and delete the last junction
            transcript_object.end = final_junction_start
            transcript_object.junctions = transcript_object.junctions[:-1]

    return transcript_dict







if __name__ == '__main__':
    msg = "This script contains functions and classes, and is not meant to be \
    run directly."
    raise ValueError(msg)
