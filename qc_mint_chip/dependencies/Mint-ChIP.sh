#!/bin/bash
###################################################################################################
# **MINT-ChIP-seq epigenetics QA/QC Singularity pipeline**
# 
# **Authors:**  Derek Ng, Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# MINT-ChIP-seq epigenetics QA/QC Singularity pipeline
# 
# **Citation:** None
# 
# **Disclaimer:**
# DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.
# 
# This material is based upon work supported by the Defense Advanced Research 
# Projects Agency under Air Force Contract No. (FA8702- 15-D-0001). Any opinions, 
# findings and conclusions or recommendations expressed in this material are 
# those of the author(s) and do not necessarily reflect the views of the 
# Defense Advanced Research Projects Agency.
# 
# © 2023 Massachusetts Institute of Technology
# 
# The software/firmware is provided to you on an As-Is basis
# 
# Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS
# Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright notice,
# U.S. Government rights in this work are defined by DFARS 252.227-7013 or
# DFARS 252.227-7014 as detailed above. Use of this work other than as specifically
# authorized by the U.S. Government may violate any copyrights that exist in this work.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
###################################################################################################
#          FILE: Mint-ChIP.sh
#         USAGE: path-to-pipe/Mint-ChIP.sh -1 <fastq1> [options]
#
#       OPTIONS: -1 <string>    The first input read filename (required option; don't type io/<filename>)
#                -2 <string>    The second input read filename (only if paired-end)
#                -n <int>       The length of the name of a fastq read; used to allocate array (default = 128)
#                -s <int>       The length of the sequence of a fastq read; used to allocate array (default = 256)
#                -t <int>       The number of threads (default = 16)
#                -h             Print this help message
#
#   DESCRIPTION: Mint-ChIP QAQC pipeline
###################################################################################################

# Pipe Start
###################################################################################################
# Read all necessary parameters and prepare data structure

# Mint-ChIP Pipeline

# Codes to change text color
green='\033[1;32m'
red='\033[1;31m'
cyan='\033[1;36m'
nc='\033[0m'
black='\033[0;30m'
cyanb='\e[106m'

# Start time
startTime=$(date +%s%N%c)
startSec=${startTime: 0:10}
startNanSec=${startTime: 10:9}
startDate=${startTime: 19}

echo -e "\n\n\n${black}${cyanb}Mint-ChIP Pipeline${nc}\n"
echo -e "${startDate}\n"

# get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ok=true

###################################################################################################
# Command line arguments
clargs() {
    while getopts 1:2:n:s:t:h opts
    do
        case "$opts" in
            1) file1="$OPTARG";; # First Fastq file
            2) file2="$OPTARG";; # Second Fastq file
            n) nameLen="$OPTARG";; # Length of name of fastq read
            s) seqLen="$OPTARG";; #length of sequence of fastq read
            t) threads="$OPTARG";; # Number of threads
            h) usage;;
            *) exit;;
        esac
    done

    types="PE"
    f1exten=false
    f2exten=false
    f1=false
    f2=false
    n=false
    s=false
    t=false

    if [ ! -f /io/$file1 ]
    then
        f1=true
    fi
    if [ ! -f /io/$file2 ]
    then
        f2=true
    fi

    if [ -z "$file1" ]
    then
        usage
    elif [ -z "$file2" ]
    then
        types="SE"
        if $f1
        then
            echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"io/\" directory.$nc"
        fi
    else
        if $f1
        then
            echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"io/\" directory.$nc"
        fi
        if $f2
        then
            echo -e "${red}The file $file2 does not exist. Please check spelling and make sure that it is in the \"io/\" directory.$nc"
        fi
    fi

    if [ -z "$nameLen" ]
    then
        nameLen=128
    elif [ $nameLen -le 0 ]
    then
        echo -e "${red}Your argument to the -n option is ${nameLen}. However, it must be a positive integer greater than 0.$nc"
        n=true
    fi

    if [ -z "$seqLen" ]
    then
        seqLen=256
    elif [ $seqLen -le 0 ]
    then
        echo -e "${red}Your argument to the -s option is ${seqLen}. However, it must be a positive integer greater than 0.$nc"
        s=true
    fi

    if [ -z "$threads" ]
    then
        threads=8
    elif [ $threads -le 0 ]
    then
        echo -e "${red}Your argument to the -t option is ${threads}. However, it must be a positive integer greater than 0.$nc"
        t=true
    fi

    if $f1exten || $f1 || $n || $s || $t
    then
        exit
    fi
}

###################################################################################################
usage() {
    echo "
Program Name: Mint-ChIP Pipeline

Usage: path-to-pipe/Mint-ChIP.sh -1 <fastq1> [options]

Options: -1 <string>    The first input read filename (required option; don't type io/<filename>)
         -2 <string>    The second input read filename (only if paired-end)
         -n <int>       The length of the name of a fastq read; used to allocate array (default = 128)
         -s <int>       The length of the sequence of a fastq read; used to allocate array (default = 256)
         -t <int>       The number of threads (default = 16)
         -h             Print this help message
"
    exit
}

# Analysis Code
# Each step would assume the previous steps have been processed.
###################################################################################################
# Step: Preparation
s_preprocess() {
    name=`echo ${file1%.fastq}`
    name=`echo ${name%.fq}`
    name=`echo ${name%_1}`
    name=`echo ${name%_R1}`
    iof="/io/Processed_"$name

    if [ ! -d /io/'Processed_'$name/ ]
    then
      mkdir /io/'Processed_'$name/
    fi
    cd /io/'Processed_'$name/

    # start record
    echo "Mint-ChIP Pipeline" > $iof/QC_pipe_processing.log
    echo "Arguments: $@" 2>&1 | tee -a $iof/QC_pipe_processing.log
    if [[ $types == SE ]]
    then
        echo "Target file is $file1" 2>&1 | tee -a $iof/QC_pipe_processing.log
    elif [[ $types == PE ]]
    then
        echo "Target files are $file1 & $file2" 2>&1 | tee -a $iof/QC_pipe_processing.log
    fi
    echo -e "Type of reads is ${types}\n" 2>&1 | tee -a $iof/QC_pipe_processing.log
}

###################################################################################################
# Step: Cutadapt
s_cutadapt() {
    echo -e "\n**************************************************"
    echo -e "${cyan}Step 1: Cutadapt - Trimming reads$nc"

    cutadapt -q 15 -m 10 -u 16 -a AGATCGGAAGAGCACACGTCTGAAC -j $threads -o $iof/'cut_'$file1 /io/$file1 \
      && echo -e "${green}Step: Cutadapt - Trimming DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
      || { echo -e "${red}Step: Cutadapt - Trimming FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    if [[ $types == PE ]]
    then
      cutadapt -q 15 -m 10 -u 16 -a AGATCGGAAGAGCGTCGTGTAGGGA -j $threads -o $iof/'cut_'$file2 /io/$file2 \
        && echo -e "${green}Step: Cutadapt - Trimming reads DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: Cutadapt - Trimming reads FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    fi
    echo "**************************************************"
}

###################################################################################################
# Step: BWA
s_bwa() {
    echo -e "\n**************************************************"
    echo -e "${cyan}Step: BWA - Align reads$nc"

    if [[ $types == PE ]]
    then
        bwa mem -t $threads /D/GRCh38 /io/$file1 /io/$file2 > $iof/$name'_bwa.sam' \
            && samtools view $iof/$name'_bwa.sam' | grep -v -e 'XA:Z:' -e 'SA:Z:' | wc -l > $iof/unique.txt \
            && echo -e "${green}Step: BWA - Align reads DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
            || { echo -e "${red}Step: BWA - Align reads FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    elif [[ $types == SE ]]
    then
         bwa mem -t $threads /D/GRCh38 /io/$file1 > $iof/$name'_bwa.sam' \
            && samtools view $iof/$name'_bwa.sam' | grep -v -e 'XA:Z:' -e 'SA:Z:' | wc -l > $iof/unique.txt \
            && echo -e "${green}Step: BWA - Align reads DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
            || { echo -e "${red}Step: BWA - Align reads FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    fi

    echo "**************************************************"
}

###################################################################################################
# Step: Samtools sort
s_samtools_sort() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Samtools Sort - Convert from SAM to BAM$nc"

    samtools view -bS $iof/$name'_bwa.sam' > $iof/$name'_bwa.bam' \
        && echo -e "${green}Step: Samtools Sort - Convert from SAM to BAM DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: Samtools Sort - Convert from SAM to BAM FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    samtools sort -@ $threads $iof/$name'_bwa.bam' -n -o $iof/$name'_bwa_sort_o.bam' \
        && echo -e "${green}Step: Samtools Sort - sort BAM DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: Samtools Sort - sort BAM FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    echo "**************************************************"
}

###################################################################################################
# Step: Samtools markdup
s_samtools_markdup() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Samtools Markdup - Remove duplicates$nc"

    samtools fixmate -m $iof/$name'_bwa_sort_o.bam' $iof/$name'_bwa_fixout.bam' \
        && samtools sort $iof/$name'_bwa_fixout.bam' -o $iof/$name'_bwa_fixmate_sort.bam' -@ $threads \
        && samtools markdup $iof/$name'_bwa_fixmate_sort.bam' $iof/$name'_bwa_rem.bam' -r -@ $threads \
        && samtools view -F 0x400 $iof/$name'_bwa_rem.bam' | grep -v -e 'XA:Z:' -e 'SA:Z:' | wc -l > $iof/dedup_unique.txt \
        && samtools index $iof/$name'_bwa_rem.bam' \
        && echo -e "${green}Step: Samtools Markdup - Remove duplicates DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: Samtools Markdup - Remove duplicates FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    echo "**************************************************"
}

###################################################################################################
# Step: Homer tag directory
s_homer_tagdir() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: HOMER Tag Directory - Make tag directory$nc"

    makeTagDirectory $name'_tagDir' $iof/$name'_bwa_rem.bam' -sspe \
        && echo -e "${green}Step: HOMER Tag Directory - Make tag directory DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: HOMER Tag Directory - Make tag directory FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    echo "**************************************************"
}

###################################################################################################
# Step: Homer UCSC
s_homer_ucsc() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: HOMER UCSC tracks - Make UCSC tracks$nc"

    makeUCSCfile $name'_tagDir' -o $name'_bwa_homer.bedgraph' \
        && echo -e "${green}Step: HOMER UCSC tracks - Make UCSC tracks DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: HOMER UCSC tracks - Make UCSC tracks FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    echo "**************************************************"
}

###################################################################################################
# Step: Homer call peaks
s_homer_peaks() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: HOMER peaks - Find peaks$nc"

    findPeaks $name'_tagDir' -o $name'_homer_peaks.txt' \
        && echo -e "${green}Step: HOMER peaks - Find peaks DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: HOMER peaks - Find peaks FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    echo "**************************************************"
}

###################################################################################################
# Step: Homer annotate peaks
s_homer_annotate_peaks() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: HOMER annotate peaks - Annotate peaks$nc"

    annotatePeaks.pl $name'_homer_peaks.txt' hg38 > $name'_home_annotatedPeaks.txt' \
        && echo -e "${green}Step: HOMER annotate peaks - Annotate peaks DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: HOMER annotate peaks - Annotate peaks FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    echo "**************************************************"
}

###################################################################################################
# Step: Sequence Lengths
s_seq_lengths() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Count reads and calculate sequence lengths$nc"

    if [[ $types == SE ]]
    then
        /S/fqSeqLen -n $nameLen -s $seqLen -o $iof/$name'_seqLen.txt' /io/$file1 \
            && echo -e "${green}Step: Count reads and calculate read lengths DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
            || { echo -e "${red}Step: Count reads and calculate read lengths FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    elif [[ $types == PE ]]
    then
        /S/fqSeqLen -n $nameLen -s $seqLen -o $iof/$name'_seqLen.txt' /io/$file1 /io/$file2 \
            && echo -e "${green}Step: Count reads and calculate read lengths DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
            || { echo -e "${red}Step: Count reads and calculate read lengths FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    fi

    echo "**************************************************"
}

###################################################################################################
# Step: Bowtie2
s_bowtie2() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Bowtie2$nc"
    bowtie2 -p $threads -x /D/grch38 -U $iof/'cut_'$file1 -S $iof/$name'-R1_bt2.sam' 2> $iof/$name'-R1_bt2.log' \
        && echo -e "${green}Step: Bowtie2 DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: Bowtie2 FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    if [[ $types == SE ]]
    then
      ln -s $name'-R1_bt2.sam' $name_bt2.sam
    elif [[ $types == PE ]]
    then
      bowtie2 -p $threads -x /D/grch38 -U $iof/'cut_'$file2 -S $iof/$name'-R2_bt2.sam' 2> $iof/$name'-R2_bt2.log' 
      bowtie2 -p $threads -x /D/grch38 -1 /io/$file1 -2 /io/$file2 -S $name'_bt2.sam' 2> $iof/$name'_bt2.log' \
          && echo -e "${green}Step: Bowtie2 DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
          || { echo -e "${red}Step: Bowtie2 FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    fi
    echo "**************************************************"
}

###################################################################################################
# Step: Genrich peak calling
s_genrich_peaks() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Genrich$nc"
    if [[ $types == SE ]]
    then
      Genrich -t $iof/$name'_bwa.sam' -o $iof/$name'.genrichPeak' -v -y \
          && echo -e "${green}Step: Genrich DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
          || { echo -e "${red}Step: Genrich FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    elif [[ $types == PE ]]
    then
      Genrich -t $iof/$name'_bwa.sam' -o $iof/$name'.genrichPeak' -v \
          && echo -e "${green}Step: Genrich DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
          || { echo -e "${red}Step: Genrich FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    fi
    echo "**************************************************"
}

###################################################################################################
# Step: MACS2
s_macs2_peaks() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: MACS2 Peak calling$nc"
    macs2 callpeak -t $iof/$name'_bwa.sam' --broad -f SAM -g hs -n $iof/$name'_macs2_B' -B --broad-cutoff 0.1 --keep-dup 1 --nomodel --extsize 100 \
        && echo -e "${green}Step: MACS2 broad peaks DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: MACS2 broad peaks FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    macs2 callpeak -t $iof/$name'_bwa.sam' --broad -f SAM -g hs -n $iof/$name'_macs2_N' -B -q 0.01 \
        && echo -e "${green}Step: MACS2 narrow peaks DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: MACS2 narror peaks FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    echo "**************************************************"
}

###################################################################################################
# Step: Blacklist
s_blacklist() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Bedtools blacklist peaks$nc"
    bedtools intersect -v -a $iof/$name'_macs2_B_peaks.broadPeak' -b /D/hg38-blacklist.sort.bed -sorted > $iof/$name'_blacklist.bed' \
        && echo -e "${green}Step: Bedtools blacklisting peaks DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: Bedtools blacklisting peaks FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    echo "**************************************************"
}

###################################################################################################
# Step: Postprocess
s_postprocess() {
    cd /io/'Processed_'$name/

    fastqc --extract /io/$file1 --outdir=/io/'Processed_'$name/
    if [[ $types == PE ]]
    then
      fastqc --extract /io/$file2 --outdir=/io/'Processed_'$name/
    fi
    grep Deduplicated /io/'Processed_'$name/*/fastqc_data.txt > /io/'Processed_'$name/deduplicated.txt

    # multiqc .

    # Sort the Bowtie2 SAM alignment file.
    samtools view -bS $iof/$name'_bt2.sam' > $iof/$name'_bt2.bam'
    # samtools sort -@ $threads $iof/$name'_bt2.bam' -o aaa > $iof/$name'_bt2_sort.bam' 
    samtools sort -@ $threads $iof/$name'_bt2.bam' -o $iof/$name'_bt2_sort.bam' 
    samtools view $iof/$name'_bt2_sort.bam' > $iof/$name'_bt2_sort.sam'

    java -Xmx64G -Xms64G -jar /S/miniconda/share/picard-slim-2.23.3-0/picard.jar MarkDuplicates I=$iof/$name'_bt2_sort.bam' O=$iof/$name'_dups.bam' M=$iof/$name'_dup_metrics.txt'

    # Parsing results for QA/QC report.
    python3 /S/QcParser.py $name mint-chip

    # End Time
    endTime=$(date +%s%N%c)
    endSec=${endTime: 0:10}
    endNanSec=${endTime: 10:9}
    endDate=${endTime: 19}
    echo $endDate

    # Check all steps completed
    if [[ $ok == true ]]
    then
        echo -e "\n${green}All steps DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log
    else
        echo -e "\n${red}Some steps FAILED$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log
    fi
    sed -i -e 's/\x1b\[[0-9;]*m//g' $iof/QC_pipe_processing.log

    echo -e "\nMint-ChIP Pipeline"
    if [[ $types == SE ]]
    then
        echo "Target file is $file1" 2>&1 | tee -a $iof/QC_pipe_processing.log
    elif [[ $types == PE ]]
    then
        echo "Target files are $file1 & $file2" 2>&1 | tee -a $iof/QC_pipe_processing.log
    fi
    echo -e "\nStart time: $startDate" 2>&1 | tee -a $iof/QC_pipe_processing.log
    echo "End time: $endDate" 2>&1 | tee -a $iof/QC_pipe_processing.log

    # Calculate time elapsed
    let elapsedNanSec=${endNanSec#"${endNanSec%%[!0]*}"}-${startNanSec#"${startNanSec%%[!0]*}"}
    millSec=${elapsedNanSec: -3}

    let elapsedSec=$endSec-$startSec
    D=$((elapsedSec/60/60/24))
    H=$((elapsedSec/60/60%24))
    M=$((elapsedSec/60%60))
    S=$((elapsedSec%60))
    printf 'Elapsed time: ' 2>&1 | tee -a $iof/QC_pipe_processing.log
    (( $D > 0 )) && printf '%d days ' $D 2>&1 | tee -a $iof/QC_pipe_processing.log
    (( $H > 0 )) && printf '%d hours ' $H 2>&1 | tee -a $iof/QC_pipe_processing.log
    (( $M > 0 )) && printf '%d minutes ' $M 2>&1 | tee -a $iof/QC_pipe_processing.log
    (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and ' 2>&1 | tee -a $iof/QC_pipe_processing.log
    printf '%d.%d seconds\n\n' $S $millSec 2>&1 | tee -a $iof/QC_pipe_processing.log
}

# run pipe ()
###################################################################################################
# step-by-step
clargs "$@"
s_preprocess "$@"
s_cutadapt
s_bwa
s_samtools_sort
s_samtools_markdup
# s_homer_tagdir
# s_homer_ucsc
# s_homer_peaks
# s_homer_annotate_peaks
s_seq_lengths
s_bowtie2
s_genrich_peaks
s_macs2_peaks
s_blacklist
s_postprocess
