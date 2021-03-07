#!/usr/bin/env python3

import argparse
import pathlib
import os
import sys
from copy import deepcopy

def get_args():
    parser = argparse.ArgumentParser(description="""
            The purpose of this script is to propogate the junction counts
            generated from STAR mapping of illumina reads against a doubled
            reference genome. This is necessary when correcting dRNAseq reads
            using FLAIR if the dRNAseq reads were mapped against many reference
            genomes. Here, dRNAseq reads that wraparound and map to multiple
            genomes need to also be corrected. With the illumina reads, the
            junctions can support the same location amongst/between all genomes
            of this extended reference.

            This script first consolidates junction counts to the first
            reference genome, althuogh the cross-genome-spanning junctions are
            kept as-is. It then propogates them forward to each similar position
            in the desired n genome lengths.
            """)

    # Required arguments
    parser.add_argument(
        '-i',
        '--in_star_file',
        type=str,
        required=True,
        help='''
        Path to the input star junction file.
        '''
    )
    parser.add_argument(
        '-o',
        '--out_star_file',
        type=str,
        required=True,
        help='''
        Path to the output, propogated star junction file.
        '''
    )
    parser.add_argument(
        '-g',
        '--genome_length',
        type=int,
        required=True,
        help='''
        Length of the viral genome. SV40: 5243 , BKPyV: 5153
        '''
    )
    parser.add_argument(
        '-n',
        '--n',
        type=int,
        required=True,
        help='''
        Number of genomes to propogate to.
        '''
    )

    # Optional arguments
    parser.add_argument(
        '-m',
        '--min_count',
        type=int,
        required=False,
        default = 5,
        help='''
        Minimum number of reads supporting the junction. If a junction is
        supported by fewer reads than this value, it is not output.
        '''
    )

    args = parser.parse_args()
    return args

#------------------------------------------------------------------------------#
# Objects
#------------------------------------------------------------------------------#
class star_junc_object:
    """
    Holds junction information from a star file.
    """

    def __init__(self, chrom, start, end, strand, motif, annotation, uniq_count, multi_count, overhang):

        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.motif = motif
        self.annotation = annotation
        self.uniq_count = uniq_count
        self.multi_count = multi_count
        self.overhang = overhang


    def output_line(self):
        output = "{} {} {} {} {} {} {} {} {}\n".format(self.chrom, self.start, self.end,
                                                       self.strand, self.motif, self.annotation,
                                                       self.uniq_count, self.multi_count, self.overhang)
        output = output.replace(" ", "\t")
        return output

#------------------------------------------------------------------------------#
# Functions for this script
#------------------------------------------------------------------------------#
def load_junction_objects(in_star_file, genome_length):
    """
    Loads each line of the input star junction file into an object. If
    the junction starts not in the first genome it is slid until it starts
    in the first genome. If a junction position, when slid, is already present
    in the dictionary, the uniq_count is added to the exiting entry. Note that
    the multi_count is roughly equal between the two because the same reads
    mapped to both positions (although the counts might be slightly different).

    Input:
        - in_star_file: Path to the input star file.

    Output:
        - Dictionary of structure dict_lookup: object, where dict_lookup is a
        tuple of struicture (start, end, strand)
    """

    star_juncs = dict()

    with open(in_star_file) as infile:
        for line in infile:

            line = line.rstrip("\n").split("\t")

            # Parse the line
            chrom, start, end, strand, motif, annotation, uniq_count, multi_count, overhang = line
            start = int(start)
            end = int(end)
            uniq_count = int(uniq_count)
            multi_count = int(multi_count)

            # Slide junctions
            while start > genome_length:
                start = start - genome_length
                end = end - genome_length

            # Generate object
            obj = star_junc_object(chrom, start, end, strand, motif, annotation, uniq_count, multi_count, overhang)

            dict_lookup = (start, end, strand)

           # Initiate dictionary if needed
            if not dict_lookup in star_juncs:
                star_juncs[dict_lookup] = obj
                continue

            # Otherwise, add uniq_count
            # NOTE - the multi_count in STAR is redundant... e.g. the two offset positions
            # will *both* get the redundant counts. So, the first mutli_count should be
            # similar (but not exactly equal) to the second multi_count of the offset position
            star_juncs[dict_lookup].uniq_count += uniq_count

    return star_juncs

def propogate_star_junctions(star_juncs, n, genome_length):
    """
    Inputs:
        - star_juncs: Dictionary of structure dict_lookup: object, where dict_lookup is a
        tuple of struicture (start, end, strand)
        - n: Interger value detailing the number of genomes to propogate into
        - genome_length: Interger value detailing the genome length.

    Outputs:
        - propogated_star_juncs: Dict of same structure as star_juncs but has propogated
        the junctions across n genomes.
    """

    propogated_star_juncs = dict()
    for dict_lookup, junction_object in star_juncs.items():

        # Propogate - non-genome spanning
        if junction_object.end < genome_length:

            for i in range(n):
                to_add = genome_length * i

                new_junction_object = deepcopy(junction_object)
                new_junction_object.start += to_add
                new_junction_object.end += to_add

                new_dict_lookup = (new_junction_object.start, new_junction_object.end, new_junction_object.strand)
                propogated_star_juncs[new_dict_lookup] = new_junction_object

        # Propogate - genome spanning (these need to go to n-1)
        else:

            for i in range(n-1):
                to_add = genome_length * i

                new_junction_object = deepcopy(junction_object)
                new_junction_object.start += to_add
                new_junction_object.end += to_add

                new_dict_lookup = (new_junction_object.start, new_junction_object.end, new_junction_object.strand)
                propogated_star_juncs[new_dict_lookup] = new_junction_object

    return propogated_star_juncs


def main():

    # Sort out the arguments
    #--------------------------------------------------------------------#
    args = get_args()
    in_star_file = args.in_star_file
    out_star_file = args.out_star_file
    genome_length = args.genome_length
    n = args.n
    min_count = args.min_count

    # Main
    #--------------------------------------------------------------------#
    print("{}: Starting script".format(sys.argv[0]))

    # Load the junctions into a dictionary, and remove redundancy
    star_juncs = load_junction_objects(in_star_file, genome_length)

    # Propogate to n genomes
    star_juncs = propogate_star_junctions(star_juncs, n, genome_length)


    # Write to output

    # Make directory if needed
    out_dir = os.path.dirname(out_star_file)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)

    with open(out_star_file, "w") as outfile:

        for dict_lookup, junction_object in star_juncs.items():

            # Filter junctions with fewer supporting reads than min_count
            if (junction_object.uniq_count + junction_object.multi_count) < min_count:
                continue

            # Write to output
            outfile.write(junction_object.output_line())

    print("{}: Finished".format(sys.argv[0]))

if __name__ == '__main__':
    main()
