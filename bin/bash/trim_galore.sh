#!/bin/bash

#SBATCH -t 0-5:00:0
#SBATCH -p short
#SBATCH --mem=10GB
#SBATCH -c 4

#------------------------------------------------------------------------------#
#Processing input and setting usage
#------------------------------------------------------------------------------#
usage() {
        echo "
        The purpose of this script is to run trimming and quality control, using
        trimgalore, on input reads. This script executes trimgalore with
        the following options:
        --stringency 1 (Even a single bp on the 3' end of a read that overlapps
          the adapter sequence will be trimmed)

        Trimming reports and fastqc reports will be in the directory of READ1_out
        in the reports/ subdirectory... e.g. $(dirname $READ1_out)/reports/*{html,zip,txt}

        IMPORTANT NOTE IMPORTANT NOTE READ ME!!!
        - READ1 AND READ2 MUST END IN .fastq.gz or .fq.gz
        - Outputs will be gzipped, so recommend ending e.g. READ1_out with .fastq.gz

        USAGE:
        $0
        -1 <READ1>
            Path to first fastq. Prefix should be .fq.gz.
        -3 <READ1_out>
            Name of the output read1 file.

        #OPTIONAL
        -2 <READ2>
            Path to second fastq, should be paired with READ1.
        -4 <READ2_out>
            Name of the output read2 file.
        -j <CORES> [1]
            Desired number of cores. Don't go past 4.
        -f <FASTQC> [T]
            Set to anything else besides default of "T" and fastqc won't be run.
        -l <MIN_LENGTH> [36]
            Minimum read length to keep
        "
}

#Setting input
while getopts 1:2:3:4:j:f:l: option ; do
        case "${option}"
        in
                1) READ1=${OPTARG};;
                2) READ2=${OPTARG};;
                3) READ1_out=${OPTARG};;
                4) READ2_out=${OPTARG};;
                j) CORES=${OPTARG};;
                f) FASTQC=${OPTARG};;
                l) MIN_LENGTH=${OPTARG};;
        esac
done
#If less than 3 option are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#------------------------------------------------------------------------------#
# Setting defaults
#------------------------------------------------------------------------------#
CORES=${CORES:-1}
FASTQC=${FASTQC:-T}
MIN_LENGTH=${MIN_LENGTH:-36}

if [[ $FASTQC == "T" ]] ; then
  FASTQC="--fastqc"
elif [[ $FASTQC != "T" ]] ; then
  FASTQC=""
fi
#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#

echo "Inputs:
READ1: $READ1
READ1_out: $READ1_out
READ2: $READ2
READ2_out: $READ2_out
CORES: $CORES
FASTQC: $FASTQC
MIN_LENGTH: $MIN_LENGTH
"

# Make output directories if needed
mkdir -p trim_galore_working $(dirname $READ1_out) $(dirname $READ2_out)

trim_galore \
  --output_dir trim_galore_working \
  --length $MIN_LENGTH \
  --stringency 1 \
  --cores $CORES \
  $FASTQC \
  $READ1 $READ2

# Rename R1
R1_PREFIX=trim_galore_working/${READ1%.fastq.gz}
R1_PREFIX=${R1_PREFIX%.fq.gz}
mv ${R1_PREFIX}_trimmed.fq.gz ${READ1_out}

# Rename R2 if there is one
if [[ $READ2 != "" ]] ; then
  R2_PREFIX=trim_galore_working/${READ2%.fastq.gz}
  R2_PREFIX=${R2_PREFIX%.fq.gz}
  mv ${R2_PREFIX}_trimmed.fq.gz ${READ2_out}
fi

# Move reports to the directory that holds R1
mkdir -p $(dirname $READ1_out)/reports
mv trim_galore_working/*{html,zip,txt} $(dirname $READ1_out)/reports && rmdir trim_galore_working
