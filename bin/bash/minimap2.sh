#!/bin/bash

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#
usage() {
        echo "
        The purpose is script is to map spliced RNAseq reads against a reference
        genome. This script works for Nanopore dRNAseq and should work well
        for Illumina short reads as long as the read lenght is above ~100nt.

        USAGE: $0

        Required:
        -i <INPUT_FASTX>
            Input file containing RNAseq reads. Can be in either fastq or
            fasta format, doesn't matter.
        -r <REF_FASTA>
            Path to the fasta containing the reference genome.
        -o <OUT_BAM>
            Path to the resultant output bam. Will contain only the mapped reads.

        Optional:
        -b <OUT_BED> ['']
            Path to the output BED file. If specified, bedtools bamtobed will be
            used to generate a bed12 file from OUT_BAM. This process will
            generate bed blocks dictated by minimap2-called introns, and will
            effectively ignore sequencing errors.
            (WARN: If ONLY_PRIMARY_READS==no, the bed
            will contain information for supplementary alignments.)

        -f <OUT_FASTA> ['']
            Path to the output fasta file. If specified, bedtools getfasta will
            be used to generate a fasta file using the OUT_BED. This is an
            effective method of pulling out 'corrected' sequencing reads. This
            fasta will be based on the BED, which only takes into account read
            start, end, and minimap2-called introns! So, no sequencing errors.
            The fasta will be of coordinate-derived transcripts.
            (WARN: If ONLY_PRIMARY_READS==no, the bed will
            contain information for supplementary alignments. These will then
            make it into the fasta, meaning some fasta entries will be
            portions of a read from supplementary alignments.)

        - P <ONLY_PRIMARY_READS> [yes]
            Options are yes or no. If set to yes, this script will not output
            secondary or supplementary alignments. If set to no, supplementary
            and secondary alignments will be allowed with the minimap2 default
            of 5 allowable secondary alignments.

        -t <THREADS> [1]
            Number of computing threads available.
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts i:r:o:b:f:P:t: option ; do
        case "${option}"
        in
                i) INPUT_FASTX=${OPTARG};;
                r) REF_FASTA=${OPTARG};;
                o) OUT_BAM=${OPTARG};;
                b) OUT_BED=${OPTARG};;
                f) OUT_FASTA=${OPTARG};;
                P) ONLY_PRIMARY_READS=${OPTARG};;
                t) THREADS=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Set defaults
#------------------------------------------------------------------------------#
OUT_BED=${OUT_BED:-""}
OUT_FASTA=${OUT_FASTA:-""}
ONLY_PRIMARY_READS=${ONLY_PRIMARY_READS:-"yes"}
THREADS=${THREADS:-1}

#------------------------------------------------------------------------------#
# Configure settings for excluding secondary/supplementary reads.
#------------------------------------------------------------------------------#
if [[ $ONLY_PRIMARY_READS == "yes" ]] ; then
  NONPRIMARY_SETTINGS="--secondary no"
elif [[ $ONLY_PRIMARY_READS == "no" ]] ; then
  NONPRIMARY_SETTINGS="--secondary yes"
else
  echo "NONPRIMARY_SETTINGS must be set to 'yes' or 'no'."
  echo "You entered $NONPRIMARY_SETTINGS"
  exit 1
fi

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
echo "
INPUTS -

INPUT_FASTX: $INPUT_FASTX
REF_FASTA: $REF_FASTA
OUT_BAM: $OUT_BAM
OUT_BED: $OUT_BED
OUT_FASTA: $OUT_FASTA
ONLY_PRIMARY_READS: $ONLY_PRIMARY_READS
THREADS: $THREADS
"


echo "Starting script."

# Make output dir if needed
mkdir -p $(dirname $OUT_BAM)

# Run minimap
echo "Running minimap2."
minimap2 \
-ax splice \
-uf \
-k14 \
--sam-hit-only \
-t $THREADS \
$NONPRIMARY_SETTINGS \
$REF_FASTA $INPUT_FASTX | samtools view -bh -@ $THREADS | samtools sort -@ $THREADS > $OUT_BAM
samtools index $OUT_BAM

# Geneate BED if desired.
if [[ $OUT_BED != "" ]] ; then
    mkdir -p $(dirname $OUT_BED)
    echo "Generating BED file."
    bedtools bamtobed \
    -bed12 \
    -i $OUT_BAM > $OUT_BED
fi

# Geneate fasta if desired.
if [[ $OUT_FASTA != "" ]] ; then
  mkdir -p $(dirname $OUT_FASTA)
  bedtools getfasta \
  -split \
  -name \
  -s \
  -fi $REF_FASTA \
  -bed $OUT_BED > $OUT_FASTA
fi

echo "Finished."
