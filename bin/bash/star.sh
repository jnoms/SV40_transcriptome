#!/bin/bash

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#
usage() {
        echo "
        The purpose is to map input short reads to a reference index using
        the STAR aligner. Settings are adapted from Kim et al's settings and
        is desired to limit penalties for non-canonical read junctions.
        USAGE: $0
        Required:
        -i <INPUT_FASTQ>
            Input file containing reads in fastq format. If this ends in .gz,
            will be handled accordingly.
        -r <REF_INDEX>
            Path to the star index directory.
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

        - P <ONLY_PRIMARY_READS> [yes]
            Options are yes or no. If set to yes, this script will not output
            secondary or supplementary alignments. If set to no, supplementary
            and secondary alignments will be allowed with the minimap2 default
            of 5 allowable secondary alignments.

        -t <THREADS> [1]
            Number of computing threads available.

        -j STAR_junction_file [""]
            If the STAR-generated junction file is desired, will extract it
            from the working directory rather than deleting it. Enter the
            desired file name.
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts i:r:o:b:P:t:j: option ; do
        case "${option}"
        in
                i) INPUT_FASTQ=${OPTARG};;
                r) REF_INDEX=${OPTARG};;
                o) OUT_BAM=${OPTARG};;
                b) OUT_BED=${OPTARG};;
                P) ONLY_PRIMARY_READS=${OPTARG};;
                t) THREADS=${OPTARG};;
                j) STAR_junction_file=${OPTARG};;
        esac
done

#------------------------------------------------------------------------------#
# Set defaults
#------------------------------------------------------------------------------#
OUT_BED=${OUT_BED:-""}
ONLY_PRIMARY_READS=${ONLY_PRIMARY_READS:-"yes"}
THREADS=${THREADS:-1}
STAR_junction_file=${STAR_junction_file:-""}

#------------------------------------------------------------------------------#
# Validate inputs
#------------------------------------------------------------------------------#
if [[ $ONLY_PRIMARY_READS == "yes" ]] ; then
  :
elif [[ $ONLY_PRIMARY_READS == "no" ]] ; then
  :
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

INPUT_FASTQ: $INPUT_FASTQ
REF_INDEX: $REF_INDEX
OUT_BAM: $OUT_BAM
OUT_BED: $OUT_BED
ONLY_PRIMARY_READS: $ONLY_PRIMARY_READS
THREADS: $THREADS
"

echo "Starting script."

# Make directories if necessary
mkdir -p $(dirname $OUT_BAM)

SAMPLE=$(basename ${OUT_BAM%.bam})

# Check if ends in .gz... need to unzip if so
# Can't use more typical approach of --readFilesCommand zcat
# because that doesn't handle symlinks generated by nextflow
if [[ $INPUT_FASTQ == *gz ]] ; then
  echo "Detected gzipped input file. Unzipping."
  gunzip -c $INPUT_FASTQ > ${SAMPLE}_unzipped
  INPUT_FASTQ=${SAMPLE}_unzipped
fi


STAR \
--runThreadN $THREADS \
--genomeDir $REF_INDEX \
--readFilesIn $INPUT_FASTQ \
--outFileNamePrefix $(dirname $OUT_BAM)/${SAMPLE}_working/${SAMPLE}_ \
--outSAMtype BAM Unsorted \
$UNZIP_COMMAND \
--outFilterType BySJout

# Extract just the bam
mv $(dirname $OUT_BAM)/${SAMPLE}_working/${SAMPLE}*bam ${OUT_BAM}_unsorted

# If the STAR-generated junction file is desired, also extract it
if [[ $STAR_junction_file != "" ]] ; then
  mkdir -p $(dirname $STAR_junction_file)
  mv $(dirname $OUT_BAM)/${SAMPLE}_working/${SAMPLE}_SJ.out.tab $STAR_junction_file
fi

# If no secondary reads are desired, remove non-primary alignments
if [[ $ONLY_PRIMARY_READS == "yes" ]] ; then
  mv ${OUT_BAM}_unsorted ${OUT_BAM}_unsorted_TMP
  samtools view -bh -F 2304 ${OUT_BAM}_unsorted_TMP > ${OUT_BAM}_unsorted && rm ${OUT_BAM}_unsorted_TMP
elif [[ $ONLY_PRIMARY_READS == "no" ]] ; then
  :
else
  echo "NONPRIMARY_SETTINGS must be set to 'yes' or 'no'."
  echo "You entered $NONPRIMARY_SETTINGS"
  exit 1
fi

# Sort and index
samtools sort ${OUT_BAM}_unsorted > ${OUT_BAM} && rm ${OUT_BAM}_unsorted
samtools index ${OUT_BAM}

# Geneate BED if desired.
if [[ $OUT_BED != "" ]] ; then
    mkdir -p $(dirname $OUT_BED)
    echo "Generating BED file."
    bedtools bamtobed \
    -bed12 \
    -i $OUT_BAM > $OUT_BED
fi

# If made an unzipped file prior to processing, delete it
if [[ -f ${SAMPLE}_unzipped ]] ; then
  rm ${SAMPLE}_unzipped
fi

echo "Finished."
