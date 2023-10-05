#!/bin/bash
#===============================================================================
# **snmC-seq epigenetics QA/QC Singularity pipeline**
# 
# **Authors:**  Derek Ng, Philip Fremont-Smith, Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# snmC-seq epigenetics QA/QC Singularity pipeline
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
#
#          FILE: snmC-seq.sh
#
#         USAGE: path-to-pipe/snmC-seq.sh -1 <read_file1> [options]
#
#       OPTIONS: -1 <string>    The first input read filename (required option; don't type io/<filename>)
#                -2 <string>    The second input read filename (only if paired-end)
#                -t <int>       The number of threads (default = 8)
#                -h             Print this help message
#
#   DESCRIPTION: snmC-seq QAQC pipeline
#
#===============================================================================


# Pipe Start
###################################################################################################
# Read all necessary parameters and prepare data structure

# snmC-seq Pipeline

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

echo -e "\n\n\n${black}${cyanb}snmC-Seq Pipeline${nc}\n"
echo -e "${startDate}\n"

# Get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ok=true

# Command line arguments
clargs() {
    while getopts 1:2:d:t:h opts
    do
        case "$opts" in
            1) file1="$OPTARG";; # First Fastq file
            2) file2="$OPTARG";; # Second Fastq file
            t) threads="$OPTARG";; # Number of threads
            h) usage;; # Print the help information
            *) exit;;
        esac
    done

    seqLen=256
    nameLen=256

    types="PE"
    f1exten=false
    f2exten=false
    f1=false
    f2=false
    d=false
    t=false

    if [[ $file1 != *.fastq ]] && [[ $file1 != *.fq ]]
    then
        f1exten=true
    fi
    if [[ $file2 != *.fastq ]] && [[ $file2 != *.fq ]]
    then
        f2exten=true
    fi
    if [ ! -f /io/fastq/$file1 ]
    then
        f1=true
    fi
    if [ ! -f /io/fastq/$file2 ]
    then
        f2=true
    fi

    if [ -z "$file1" ]
    then
        usage
    elif [ -z "$file2" ]
    then
        types="SE"
        if $f1exten
        then
            echo -e "${red}The name of your file is $file1. However, the extension must be .fastq or .fq.$nc"
        fi
        if $f1
        then
            echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"/io/fastq/\" directory.$nc"
        fi
    else
        if $f1exten
        then
            echo -e "${red}The name of your first file is $file1. However, the extension must be .fastq or .fq.$nc"
        fi
        if $f2exten
        then
            echo -e "${red}The name of your second file is $file2. However, the extension must be .fastq or .fq.$nc"
        fi

        if $f1
        then
            echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"/io/fastq/\" directory.$nc"
        fi
        if $f2
        then
            echo -e "${red}The file $file2 does not exist. Please check spelling and make sure that it is in the \"/io/fastq/\" directory.$nc"
        fi
    fi

    if [ -z "$threads" ]
    then
        threads=8
    elif [ $threads -le 0 ]
    then
        echo -e "${red}Your argument to the -t option is ${threads}. However, it must be a positive integer greater than 0.$nc"
        t=true
    fi

    if $f1exten || $f1 || $t
    then
        exit
    fi
}

usage() {
    echo "
Program Name: snmC-seq Pipeline

Usage: path-to-pipe/snmC-seq.sh -1 <read_file1> [options]

Options: -1 <string>    The first input read filename (required option; don't type io/fastq/<filename>)
         -2 <string>    The second input read filename (only if paired-end)
         -d <string>    Directory name
         -t <int>       The number of threads (default = 8)
         -h             Print this help message
"
    exit
}

# Analysis Code
# Each step would assume the previous steps have been processed.
###################################################################################################
# Step: Preparation
s_preprocess() {
    name=`echo ${file1%.raw.fq}`
    name=`echo ${name%.fq}`
    name=`echo ${name%-R1}`
    name=`echo ${name%-R1.raw.fq}`

    echo $name

    # if [ -d /io/'Processed_'$name/ ]
    # then
    #     while true
    #     do
    #         echo -e "${red}The directory /io/Processed_${name}/ exists. Do you wish to overwrite? (y/n)$nc"
    #         read yn
    #         case $yn in
    #             [Yy]*) rm -r /io/'Processed_'$name/; break;;
    #             [Nn]*) echo "Please rename or move the existing directory. You may also rename the input file(s)."; exit;;
    #             *) echo "Please answer y or n.";;
    #         esac
    #     done
    # fi
    if [ ! -d /io/'Processed_'$name/ ]
    then
      mkdir /io/'Processed_'$name/
    fi
    cd /io/'Processed_'$name/
    pathdir='/io/Processed_'$name

    # start record
    echo "snmC-seq Pipeline" > QC_pipe_processing.log
    echo "Arguments: $@" 2>&1 | tee -a QC_pipe_processing.log
    if [[ $types == SE ]]
    then
        echo "Target file is $file1" 2>&1 | tee -a QC_pipe_processing.log
    elif [[ $types == PE ]]
    then
        echo "Target files are $file1 & $file2" 2>&1 | tee -a QC_pipe_processing.log
    fi
    echo -e "Type of reads is ${types}\n" 2>&1 | tee -a QC_pipe_processing.log
}

###################################################################################################
# Step: Cutadapt
s_cutadapt() {
    echo -e "\n**************************************************"
    echo -e "${cyan}Step 1: Cutadapt - Trimming reads$nc"

    cd /io/

    cutadapt -q 20 -m 30 -u 16 -a AGATCGGAAGAGCACACGTCTGAAC -j $threads -o $pathdir/'cut_'$name'-R1.fq' /io/fastq/$file1 \
      && echo -e "${green}Step: Cutadapt - Trimming DONE$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
      || { echo -e "${red}Step: Cutadapt - Trimming FAIL$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    if [[ $types == PE ]]
    then
      cutadapt -q 20 -m 30 -u 16 -a AGATCGGAAGAGCGTCGTGTAGGGA -j $threads -o $pathdir/'cut_'$name'-R2.fq' /io/fastq/$file2 \
        && echo -e "${green}Step: Cutadapt - Trimming reads DONE$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
        || { echo -e "${red}Step: Cutadapt - Trimming reads FAIL$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    fi
    echo "**************************************************"
}

###################################################################################################
# Step: Bismark
s_bismark() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step 2: Bismark - Map bisulfite$nc"

    bismark --bowtie2 --genome /D/genome --pbat -p $threads $pathdir/'cut_'$name'-R1.fq' \
        && echo -e "${green}Step: Bismark - Map bisulfite DONE$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
        || { echo -e "${red}Step: Bismark - Map bisulfite FAIL$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    if [[ $types == PE ]]
    then
      bismark --bowtie2 --genome /D/genome -p $threads $pathdir/'cut_'$name'-R2.fq' \
        && echo -e "${green}Step: Bismark - Map bisulfite DONE$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
        || { echo -e "${red}Step: Bismark - Map bisulfite FAIL$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    fi
    echo "**************************************************"
}

###################################################################################################
# Step: Samtools sort
s_samtools_sort() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step 3: Samtools Sort - Sorting$nc"

    samtools sort -O 'bam' -o $pathdir/$name'-R1.sort.bam' -T temp_aln -@ $threads /io/'cut_'$name'-R1_bismark_bt2.bam' \
        && echo -e "${green}Step: Samtools Sort - Sorting R1 DONE$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
        || { echo -e "${red}Step: Samtools Sort - Sorting R1 FAIL$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    if [[ $types == PE ]]
    then
      samtools sort -O 'bam' -o $pathdir/$name'-R2.sort.bam' -T temp_aln -@ $threads /io/'cut_'$name'-R2_bismark_bt2.bam' \
        && echo -e "${green}Step: Samtools Sort - Sorting R2 DONE$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
        || { echo -e "${red}Step: Samtools Sort - Sorting R2 FAIL$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    fi
    echo "**************************************************"
}

###################################################################################################
# Step: Picard
s_picard() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Picard - Mark duplicates$nc"

    picard MarkDuplicates I=$pathdir/$name'-R1.sort.bam' O=$pathdir/$name'-R1.markedDups.bam' M=marked_dup_metrics.txt REMOVE_DUPLICATES=true \
        && echo -e "${green}Step: Picard - R1 Mark duplicates DONE$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
        || { echo -e "${red}Step: Picard - R1 Mark duplicates FAIL$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    if [[ $types == PE ]]
    then
      picard MarkDuplicates I=$pathdir/$name'-R2.sort.bam' O=$pathdir/$name'-R2.markedDups.bam' M=marked_dup_metrics.txt REMOVE_DUPLICATES=true \
        && echo -e "${green}Step: Picard - R2 Mark duplicates DONE$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
        || { echo -e "${red}Step: Picard - R2 Mark duplicates FAIL$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    fi
    echo "**************************************************"
}

###################################################################################################
# Step: Samtools view
s_samtools_view() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Samtools View - Filter non-clonal reads$nc"

    samtools view -b -q 10 $pathdir/$name'-R1.markedDups.bam' -o $pathdir/$name'-R1.bam' -@ $threads \
        && echo -e "${green}Step: Samtools View - R1 Filter non-clonal reads DONE$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
        || { echo -e "${red}Step: Samtools View - R1 Filter non-clonal reads FAIL$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    if [[ $types == PE ]]
    then
      samtools view -b -q 10 $pathdir/$name'-R2.markedDups.bam' -o $pathdir/$name'-R2.bam' -@ $threads \
        && echo -e "${green}Step: Samtools View - R2 Filter non-clonal reads DONE$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
        || { echo -e "${red}Step: Samtools View - R2 Filter non-clonal reads FAIL$nc" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    fi
    echo "**************************************************"
}

###################################################################################################
# Step: Sequence Lengths
s_seq_lengths() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Count reads and calculate sequence lengths$nc\n"

    cd /io/'Processed_'$name/
    /S/fqSeqLen -n $nameLen -s $seqLen -o $pathdir/$name'_seqLen.txt' 'cut_'$name'-R1.fq' \
      && echo -e "${green}Step: Count reads and calculate read lengths DONE$nc\n" 2>&1 | tee -a $pathdir/QC_pipe_processing.log \
      || { echo -e "${red}Step: Count reads and calculate read lengths FAIL$nc\n" 2>&1 | tee -a $pathdir/QC_pipe_processing.log && ok=false; }
    echo "**************************************************"
}

###################################################################################################
# Step: Postprocess
s_postprocess() {
    cd /io/'Processed_'$name/

    reportname='cut_'$name'-R1_bismark_bt2_SE_report.txt'
    grep "Sequences analyzed in total" /io/$reportname > snmC_log.txt
    grep "unique best hit from" /io/$reportname >> snmC_log.txt
    grep "Mapping efficiency" /io/$reportname >> snmC_log.txt
    grep "C methylated in CpG context" /io/$reportname >> snmC_log.txt
    grep "C methylated in CH" /io/$reportname >> snmC_log.txt

    if [[ $types == PE ]]
    then
      reportname2='cut_'$name'-R2_bismark_bt2_SE_report.txt'
      grep "Sequences analyzed in total" /io/$reportname2 > snmC_log2.txt
      grep "unique best hit from" /io/$reportname2 >> snmC_log2.txt
      grep "Mapping efficiency" /io/$reportname2 >> snmC_log2.txt
      grep "C methylated in CpG context" /io/$reportname2 >> snmC_log2.txt
      grep "C methylated in CH" /io/$reportname2 >> snmC_log2.txt
    fi

    # Calculate CCC coverage.
    ccc_total=$(zgrep -P -c '.*\t*\d+\t[+|-]\tCCC\t\d+\t\d+\t\d+' '/io/allc/'$name'.allc.tsv')
    ccc_meth=$(zgrep -P -c '.*\t*\d+\t[+|-]\tCCC\t[^0]+\t\d+\t\d+' '/io/allc/'$name'.allc.tsv')
    VAR=$(echo "scale=10; $ccc_meth/$ccc_total" | bc)
    # echo $1','$VAR >> allc_ccc.txt
    echo $VAR >> snmC_log.txt
    echo $VAR >> snmC_log2.txt

    # Calculate genome coverage
    samtools sort -@ 5 '/io/cut_'$name'-R1_bismark_bt2.bam' -o 'sort_'$name'-R1.bam'    
    bedtools genomecov -ibam 'sort_'$name'-R1.bam' | grep -P "genome\t0\t" >> snmC_log.txt
    if [[ $types == PE ]]
    then
      samtools sort -@ 5 '/io/cut_'$name'-R2_bismark_bt2.bam' -o 'sort_'$name'-R2.bam'    
      bedtools genomecov -ibam 'sort_'$name'-R2.bam' | grep -P "genome\t0\t" >> snmC_log2.txt
    fi

    echo "QcParser running"
    # Parsing results for QA/QC report.
    python3 /S/QcParser.py $name snmC-seq
    echo "QcParser done"

    # End time
    endTime=$(date +%s%N%c)
    endSec=${endTime: 0:10}
    endNanSec=${endTime: 10:9}
    endDate=${endTime: 19}
    echo $endDate

    # Check all steps completed
    if [[ $ok == true ]]
    then
        echo -e "\n${green}All steps DONE$nc" 2>&1 | tee -a QC_pipe_processing.log
    else
        echo -e "\n${red}Some steps FAILED$nc" 2>&1 | tee -a QC_pipe_processing.log
    fi
    sed -i -e 's/\x1b\[[0-9;]*m//g' QC_pipe_processing.log

    echo -e "\nsnmC-seq Pipeline"
    if [[ $types == SE ]]
    then
        echo "Target file is $file1" 2>&1 | tee -a QC_pipe_processing.log
    elif [[ $types == PE ]]
    then
        echo "Target files are $file1 & $file2" 2>&1 | tee -a QC_pipe_processing.log
    fi
    echo -e "\nStart time: $startDate" 2>&1 | tee -a QC_pipe_processing.log
    echo "End time: $endDate" 2>&1 | tee -a QC_pipe_processing.log

    # Calculate time elapsed
    let elapsedNanSec=${endNanSec#"${endNanSec%%[!0]*}"}-${startNanSec#"${startNanSec%%[!0]*}"}
    millSec=${elapsedNanSec: -3}

    let elapsedSec=$endSec-$startSec
    D=$((elapsedSec/60/60/24))
    H=$((elapsedSec/60/60%24))
    M=$((elapsedSec/60%60))
    S=$((elapsedSec%60))
    printf 'Elapsed time: ' 2>&1 | tee -a QC_pipe_processing.log
    (( $D > 0 )) && printf '%d days ' $D 2>&1 | tee -a QC_pipe_processing.log
    (( $H > 0 )) && printf '%d hours ' $H 2>&1 | tee -a QC_pipe_processing.log
    (( $M > 0 )) && printf '%d minutes ' $M 2>&1 | tee -a QC_pipe_processing.log
    (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and ' 2>&1 | tee -a QC_pipe_processing.log
    printf '%d.%d seconds\n\n' $S $millSec 2>&1 | tee -a QC_pipe_processing.log
}

# run pipe ()
###################################################################################################
# step-by-step
clargs "$@"
s_preprocess "$@"
s_cutadapt
s_bismark
s_samtools_sort
s_picard
s_samtools_view
s_seq_lengths
s_postprocess
