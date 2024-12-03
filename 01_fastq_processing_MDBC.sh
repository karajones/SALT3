#!/bin/bash
####################################################################################################
### This script will run through the steps needed to clean raw fastq files in preparation for    
### running DADA2 in R. It's also just a good baseline script for cleaning sequencing data       
### before any analysis. This script should work on both Linux and Mac, but has only been
### tested on a Mac.
### Kara Jones
### ksjones@usgs.gov
### October 2023
###
###
### Install the following dependencies before running:
### - FastQC: https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc
### - cutadapt: https://cutadapt.readthedocs.io/en/stable/installation.html
### - fastp: https://github.com/OpenGene/fastp#get-fastp
### - seqkit: https://bioinf.shenwei.me/seqkit/
### - rename: brew install rename
###
### Instructions for running the script:
### 1. Make sure the dependencies are installed (see list and links below)
### 2. Make the indicated changes to variables (e.g., primer sequences)
### 3. Run the script in the same directory that contains the folder of FASTQ files
### 
### whatever_directory <- run the script in this directory
### ├── fastqs
### │   ├── 18S-SA22-036
### │   │   ├── 18S-SA22-036_S169_L001_R1_001.fastq.gz
### │   │   └── 18S-SA22-036_S169_L001_R2_001.fastq.gz
### │   ├── 18S-SA22-037
### │   │   ├── 18S-SA22-037_S170_L001_R1_001.fastq.gz
### │   │   └── 18S-SA22-037_S170_L001_R2_001.fastq.gz
### ... etc ...
###
### Note: Names of the files and directories mostly do not matter, but the script makes two assumptions:
### 1. Fastqs have the extension .fastq.gz
### 2. The extra junk we want to get rid of in the names is separated from the part of the name 
### we want to keep by an underscore (i.e. 18S-SA22-036_S169_L001_R1_001.fastq.gz --> 18S-SA22-036.fq.gz)
###
### OUTPUT:
### - raw: all raw .fastq.gz files (unmodified)
### - deduplicated: deduplicated and error-correct reads (intermediate step)
### - trimmed: reads trimmed of primers (intermediate step)
### - clean: final reads
### - stats: all files with information on read quality
###			- fastp_reports: html reports on intial (raw) read quality for each sample
###			- fastqc_clean: html reports on clean reads for each sample
###			- read_stats: tsv files with read stats for each stage of processing
###
### If there are still quality issues after processing, make changes to the settings and rerun the script.
### However, old files WILL be overwritten by any new data!!
###
####################################################################################################

# Primer sequences:
# CO1 primers: 'F-GGWACWGGWTGAACWGTWTAYCCYCC' and 'R-TAAACTTCAGGGTGACCAAAAAATCA'
# 18S metazoan primers: 'F-GTACACACCGCCCGTC' and 'R-TGATCCTTCTGCAGGTTCACCTAC'
# coral DEEP: 17F_gomx_deepcor, 5’-CGTGAAACCGYTRRAAGGG-3’ 631R_gom_cor, 5’- TTGGTCCGTGTTTCAAGACG-3’
# MiFish-U-F: GTGTCGGTAAAACTCGTGCCAGC MiFish-U-R: CATAGTGGGGTATCTAATCCCAGTTTG

# Other settings
# MINLENGTH="100" for all primer sets except 18S, which uses MINLENGTH="50"
# MAXLENGTH="300" for all primer sets

### CHANGE THESE PARAMETERS IF NEEDED ###

# PRIMERS
# Use 5'-3' sequence; do NOT reverse-complement the reverse primer!
FPRIMER="GTGTCGGTAAAACTCGTGCCAGC"
RPRIMER="CATAGTGGGGTATCTAATCCCAGTTTG"

# LENGTH LIMIT FOR FINAL READS
# Throw out reads if the are shorter than MINLENGTH or longer than MAXLENGTH
# This will be set for R1 and R2s (futher trimming of low quality reads will be done in DADA2 later)
MINLENGTH="100"
MAXLENGTH="300"


### DON'T CHANGE ANYTHING BELOW THIS UNLESS YOU KNOW WHAT YOU ARE DOING! :-) ###

# SET STANDARD VARIABLES
# Number of threads available
NPROC=$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)

# DEPENDENCY CHECK
echo -n "Checking dependencies... "
for name in fastp seqkit cutadapt rename fastqc
do
  [[ $(which $name 2>/dev/null) ]] || { 
    echo -en "\n$name cannot be found.";deps=1; 
    }
done
    [[ $deps -ne 1 ]] && echo "OK" || { 
        echo -en "\nPlease install or activate the above dependencies and rerun this script.\n";exit 1; 
}


# Make necessary directories
mkdir -p raw stats stats/raw stats/clean stats/read_stats deduplicated trimmed clean 


# Move all files from individual directories into a single directory

find . -name "*.fastq.gz" ! -path "./raw/*" -exec cp {} ./raw/ \;

# rename files
for f in ./raw/*_R1_*.fastq.gz
    do mv $f ${f%%_*}_R1.fastq.gz
done
for f in ./raw/*_R2_*.fastq.gz
    do mv $f ${f%%_*}_R2.fastq.gz
done


# run fastqc to generate html reports of raw reads
fastqc -t $NPROC --svg --noextract -o ./stats/raw ./raw/*.fastq.gz

# Deduplicate, error correct, and detect any remaining adapters with fastp (also output quality report)
### note: fastp doesn't recognize non-ACTG characters so primers must be removed after deduplication
for f in ./raw/*_R1.fastq.gz
    do

    base=$(basename ${f%%_*})    
    R1=${base}_R1
    R2=${base}_R2

    echo -en "\nProcessing sample $base ..."

    # first fastp call to deduplicate and error-correct reads
    fastp --in1 ./raw/$R1.fastq.gz --in2 ./raw/$R2.fastq.gz --out1 ./deduplicated/$R1.dedup.fq.gz --out2 ./deduplicated/$R2.dedup.fq.gz --dedup --dup_calc_accuracy 6 --correction --thread=$NPROC -j ./stats/raw/$R1.fastp.json

    # cutadapt call to remove primers
    # note: when json is supported by multiqc, `add --json=./stats/raw/$R1.cutadapt.json` to call and remove `> ./stats/raw/$R1.cutadapt.log`
    cutadapt -e 1.5 -m 1 --cores=0 -g $FPRIMER -G $RPRIMER -o ./trimmed/$R1.trimmed.fq.gz -p ./trimmed/$R2.trimmed.fq.gz ./deduplicated/$R1.dedup.fq.gz ./deduplicated/$R2.dedup.fq.gz > ./stats/raw/$R1.cutadapt.log

    # second fastp call to clean reads
    fastp --in1 ./trimmed/$R1.trimmed.fq.gz --in2 ./trimmed/$R2.trimmed.fq.gz --out1 ./clean/$R1.clean.fq.gz --out2 ./clean/$R2.clean.fq.gz --dedup --dup_calc_accuracy 6 --trim_poly_x --cut_front --cut_tail --cut_window_size=5 --cut_mean_quality=30 --qualified_quality_phred=30 --unqualified_percent_limit=20 --length_required=$MINLENGTH --length_limit=$MAXLENGTH --n_base_limit=5 --thread=$NPROC -j ./stats/clean/$R1.fastp.json

	done

echo "Generating read statistics..."

# run fastqc to generate html reports of clean reads
fastqc -t $NPROC --svg --noextract -o ./stats/clean ./clean/*.clean.fq.gz

# compile stats from fastqc, fastp and cutadapt
multiqc --interactive -n multiqc_report_raw_reads.html -o ./stats -i "Read quality before processing" ./stats/raw/.
multiqc --interactive -n multiqc_report_clean_reads.html -o ./stats -i "Read quality after processing" ./stats/clean/.

seqkit stats -aT ./raw/*.fastq.gz > ./stats/read_stats/raw_reads.stats.tsv 
seqkit stats -aT ./deduplicated/*.dedup.fq.gz > ./stats/read_stats/dedup.stats.tsv
seqkit stats -aT ./trimmed/*.trimmed.fq.gz > ./stats/read_stats/trimmed.stats.tsv
seqkit stats -aT ./clean/*.clean.fq.gz > ./stats/read_stats/clean.stats.tsv

# cleanup

echo "Finished!"


