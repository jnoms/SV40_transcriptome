# SV40_transcriptome

The code in this repo is sufficient to reproduct the analyses present in the manuscript "Long-read sequencing reveals complex patterns of wraparound transcription in polyomaviruses", PLOS Pathogens 2022 - https://doi.org/10.1371/journal.ppat.1010401

## About this repo  
There are several major components to this repo:
  1. A nextflow pipeline to process short-read (Illumina) or long-read (Nanopore and PacBio) RNAseq data generated from polyomavirus-infected cells. 
  2. Stand-alone bash and python scripts for individual steps.
  3. Google colab documents for reproducible downstream analysis and plotting in R - these documents are sufficient to reproduce all main figures in the manuscript.  

In general, I consider the nextflow pipeline to be present for transparency of pipeline steps rather than for reuse. However, the Nextflow pipeline is reproducible and comes with a pre-configured conda environment containing all dependencies. Furthermore, virus-specific config files are present in resources/virus_configs and associated references genomes are located in resources/ref. There are likewise three Nextflow processes - one for running the pipeline on linux using slrum (`-process o2_slurm`), one for running it on linux but not using slurm (`--process o2_local`), and one for running the pipeline locally on mac (`-process mac`). **In general, I recommend using the flexible stand-alone bash and python scripts for your analyses.** Therefore, I won't describe the nextflow pipeline parameters in more depth (although a skilled user should be able to get it working with minimal trouble).

## Description of the main stand-alone modules  
The heart of this repository is a workflow to map short- or long-read RNAseq to a reference genome and produce a tidy output file that defines, for each read, what splice junctions are present and to group reads into transcript classes based on shared splice junctions.  

The basic pipeline goes like this:  
  1. Map the reads using `minimap.sh` or `star.sh`. Typically I recommend using Minimap unless your reads are too short (<100bp), in which case Minimap does not behave well and you should use STAR. These scripts can output a bed file using the `-b` switch is specified. This bed file is the goal, rather than the bam file the script will also produce. **The bed file will list where on the reference genome each read starts and end, as well as list the exons present in the read (which can be inferred from the blockStarts and blockSize columns in the bed file).** Notably, exons/introns are called based on minimap- or STAR-called introns.
  2. Convert the BED file from step 1 into a tidy "span" file using `bed_to_span.py`. **The output file will contain an individual line for each exon or intron present in each read as defined from the bed file, and will indiciate their start and end coordinates.** Thus, a spliced read will be present on multiple lines. This is a really convenient file for using R and ggplot (or python, if that's your cup of tea) for analyzing and plotting the data. If you're using short-RNAseq reads annd you're just looking for splice sites, the start and end position of each read don't really matter, so to vastly reduce the output file size you can use the `--junc_reads_only` flag to only output information on the introns.    

Detailed help information will be produced by the bash scripts when executed with no arguments and by the python scripts when the -h flag is used. 

You'll notice that there are many other python scripts. These are mostly for polyomavirus-analysis-specific tasks (e.g. correcting/sliding BED files when mapped against a reference genome that contains multiple copies of itself). The scripts have a detailed help description.  

## Dependencies  
There are pretty few dependencies - just STAR or Minimap2, bedtools, samtools, and python 3. There are several conda yaml files present in resources/conda which will generate conda environments that definitely work - they come with a lot of additional software besides these core dependencies, so I recommend making your own environment.  
