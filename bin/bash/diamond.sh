#!/bin/bash

#SBATCH -t 0-1:00:0
#SBATCH -p priority
#SBATCH --mem=10GB
#SBATCH -c 2

#------------------------------------------------------------------------------#
# Defining usage and setting input
#------------------------------------------------------------------------------#

usage() {
        echo "
        This purpose of this file is to use DIAMOND to map a file containing
        query sequences against a pre-determined protein database. DIAMOND will
        automatically detect and use # of threads available.

        USAGE: $0

        Required:
        -d <DATABASE>    Path to the diamond database.
        -q <QUERY>       Path to the query file containing nucleotide sequences.
        -o <OUTPUT_PATH> Path to the output file.
        -k <MAX_TARGETS> Maximum number of subjects to report. Setting to 0
                         reports all alignments. Otherwise, 25 is a good default.

        Optional:
        -m <MEMORY>      Memory available in GB. DIAMOND blocksize is set to
                          Default: 10
        -t <TEMP_DIR>    Path to the temporary directory
                          Default: ./diamond_temp
        -e <EVALUE>      Max evalue.
                          Default: 10
        -f <OUT_FORMAT>  Output format.
                          Default: '6 qseqid stitle sseqid evalue bitscore
                                    pident length qlen slen qstart qend sstart
                                    send'
                          - length == alignment length
        -l <LOG_FILE>   Path to the log file
                           Default: No log file.
        -s <SAMPLE_ID>  Name of the sample, used for logging purposes.
                           Default: Basename of query.
        -T <DIAMOND_TYPE> In what mode to run DIAMOND. Options are 'blastx' or
                        'blastp'.
                           Default: 'blastp'
        -u <INCLUDE_UNALIGNED> Options are 'TRUE' or 'FALSE'. If TRUE, will
                               output a mostly-blank line for queries that
                               didnt have an alignment.
                               Default: TRUE
        "
}

#If less than 3 options are input, show usage and exit script.
if [ $# -le 3 ] ; then
        usage
        exit 1
fi

#Setting input
while getopts d:q:o:m:t:e:f:l:s:T:k:u: option ; do
        case "${option}"
        in
                d) DATABASE=${OPTARG};;
                q) QUERY=${OPTARG};;
                o) OUTPUT_PATH=${OPTARG};;
                m) MEMORY=${OPTARG};;
                t) TEMP_DIR=${OPTARG};;
                e) EVALUE=${OPTARG};;
                f) OUT_FORMAT=${OPTARG};;
                l) LOG_FILE=${OPTARG};;
                s) SAMPLE_ID=${OPTARG};;
                T) DIAMOND_TYPE=${OPTARG};;
                k) MAX_TARGETS=${OPTARG};;
                u) INCLUDE_UNALIGNED=${OPTARG};;
        esac
done


#------------------------------------------------------------------------------#
# Setting defaults
#------------------------------------------------------------------------------#
MEMORY=${MEMORY:-10}
TEMP_DIR=${TEMP_DIR:-./diamond_temp}
EVALUE=${EVALUE:-10}
OUT_FORMAT=${OUT_FORMAT:-"6 qseqid sseqid evalue bitscore pident length qlen slen qstart qend sstart send"}
SAMPLE_ID=${SAMPLE_ID:-$(basename $QUERY)}
DIAMOND_TYPE=${DIAMOND_TYPE:-blastp}
INCLUDE_UNALIGNED=${INCLUDE_UNALIGNED:-TRUE}

# Calculate DIAMOND block size
if [[ $MEMORY < 10 ]] ; then
  BLOCK=1
else
  BLOCK=$((MEMORY/10))
fi

# Set DIAMOND type
if [[ $DIAMOND_TYPE != "blastx" ]] && [[ $DIAMOND_TYPE != "blastp" ]] ; then
  echo "DIAMOND type must be 'blastx' or blastp.''"
  echo "You entered $DIAMOND_TYPE"
  exit 1
fi

if [[ $MAX_TARGETS == '' ]] ; then
    echo "MAX_TARGETS, -k, is a required param."
    exit 1
fi

# Optional unaligned setting. Make a switch if necessary.
if [[ $INCLUDE_UNALIGNED == "FALSE" ]] ; then
  UNALIGNED_SWITCH_LINE=""
elif [[ $INCLUDE_UNALIGNED == "TRUE" ]] ; then
  UNALIGNED_SWITCH_LINE="--unal 1"
else
  echo "INCLUDE_UNALIGNED must be set to TRUE or FALSE."
  echo "You entered $INCLUDE_UNALIGNED"
  exit 1
fi

#------------------------------------------------------------------------------#
# Main
#------------------------------------------------------------------------------#
write_log() {
  # The purpose of this function is to append to a log, which is CSV output.
  # Structure is sampleID, file_name/process, message, error(0 or 1). Log file
  # is optional - if not input, will print to screen.

  # Usage:
  # write_log <sampleID> <message> <error - 0 or 1> <log_file>
  sampleID=$1
  FILE_NAME=$(basename $0)
  message=$2
  error=$3
  log_file=$4

  OUTPUT="${sampleID},${FILE_NAME},${message},${error}"
  echo $OUTPUT
  date

  # If no log file, leave
  if [[ $log_file == '' ]] ; then
    return
  fi

  # Write to log
  mkdir -p $(dirname $log_file)
  echo $OUTPUT >> $log_file
}

echo "
INPUTS -

DATABASE: $DATABASE
QUERY: $QUERY
OUTPUT_PATH: $OUTPUT_PATH
MEMORY: $MEMORY
TEMP_DIR: $TEMP_DIR
EVALUE: $EVALUE
OUT_FORMAT: $OUT_FORMAT
LOG_FILE: $LOG_FILE
SAMPLE_ID: $SAMPLE_ID
DIAMOND_TYPE: $DIAMOND_TYPE
MAX_TARGETS: $MAX_TARGETS
INCLUDE_UNALIGNED: $INCLUDE_UNALIGNED
"

#Make sure all necesssary directories have been made
mkdir -p $( dirname $OUTPUT_PATH )
mkdir -p $TEMP_DIR

#Write to log
write_log \
$SAMPLE_ID \
"Starting at $(date)." \
0

#Make sure inputs exist exists
if [ ! -f $DATABASE ] ; then
  write_log \
  $SAMPLE_ID \
  "Cannot find the input database, $DATABASE. Exiting." \
  1 \
  $LOG_FILE

  exit 1
elif [ ! -f $QUERY ] ; then
  write_log \
  $SAMPLE_ID \
  "Cannot find the query file, $QUERY. Exiting." \
  1 \
  $LOG_FILE

  exit 1
fi

#Running DIAMOND
diamond $DIAMOND_TYPE \
-d $DATABASE \
-q $QUERY \
--more-sensitive \
-o $OUTPUT_PATH \
--tmpdir $TEMP_DIR \
--evalue $EVALUE \
--outfmt $OUT_FORMAT \
--index-chunks 1 \
--block-size $BLOCK \
$UNALIGNED_SWITCH_LINE \
--max-target-seqs $MAX_TARGETS \
--max-hsps 10

# Remove temp directory
rm -r $TEMP_DIR

#Check that the output file exists. It no matches it should exist, though it will be empty.
if [ ! -f $OUTPUT_PATH ] ; then
  write_log \
  $SAMPLE_ID \
  "CANNOT find DIAMOND output file. Exiting script." \
  1 \
  $LOG_FILE
  exit 1
fi

# If there is no contents in the OUTPUT_FILE don't quit, but report warning
if [[ ! -s $OUTPUT_PATH ]] ; then
  write_log \
  $SAMPLE_ID \
  "WARN: Output file $OUTPUT_PATH doesn't have any contents." \
  1 \
  $LOG_FILE
fi

write_log \
$SAMPLE_ID \
"Finished at $(date)." \
0
