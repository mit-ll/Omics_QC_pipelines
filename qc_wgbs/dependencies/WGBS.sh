#!/bin/bash
###################################################################################################
# **WGBS-seq epigenetics QA/QC Singularity pipeline**
# 
# **Authors:**  Derek Ng, Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# WGBS-seq epigenetics QA/QC Singularity pipeline
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
#          FILE: WGBS.sh
#
#         USAGE: path-to-pipe/WGBS.sh -1 <read_file1> [options]
#       Options: -1 <string>    The first input read filename (don't type io/<filename>)
#                -2 <string>    The second input read file (only if paired-end)
#                -a <string>    Sequence of an adapter ligated to the 3' end of the first read (default = AGATCGGAAGAGC)
#                -b <string>    3' adapter to be removed from second read in a pair (default = AGATCGGAAGAGC)
#                -t <int>       The number of threads (default = 16)
#                -h             Print this help message
#
#   DESCRIPTION: WGBS QAQC pipeline
###################################################################################################


# Pipe Start
###################################################################################################
# Read all necessary parameters and prepare data structure

# WGBS Pipeline

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

echo -e "\n\n\n${black}${cyanb}WGBS Pipeline${nc}\n"
echo -e "${startDate}\n"

pipe_version='v1.00'
host="mitll/wgbs pipe"

# get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
md5=`md5sum $0 | awk '{print $1}'`
ok=true

###################################################################################################
# read parameters
clargs() {
    while getopts 1:2:a:b:t:h opts
    do
        case "$opts" in
            1) file1="$OPTARG";; # First Fastq file
            2) file2="$OPTARG";; # Second Fastq file
            a) adapter_1="$OPTARG";; # Add adapter1
            b) adapter_2="$OPTARG";; # Add adapter2
            t) threads="$OPTARG";; # Number of threads
            h) usage;; # Print the help information
            *) exit;;
        esac
    done

    types="PE"
    f1exten=false
    f2exten=false
    f1=false
    f2=false
    t=false
    nameLen=128
    seqLen=256

    if [[ $file1 != *.fastq ]] && [[ $file1 != *.fq ]]
    then
        f1exten=true
    fi
    if [[ $file2 != *.fastq ]] && [[ $file2 != *.fq ]]
    then
        f2exten=true
    fi
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
        if $f1exten
        then
            echo -e "${red}The name of your file is $file1. However, the extension must be .fastq or .fq.$nc"
        fi
        if $f1
        then
            echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"io/\" directory.$nc"
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
            echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"io/\" directory.$nc"
        fi
        if $f2
        then
            echo -e "${red}The file $file2 does not exist. Please check spelling and make sure that it is in the \"io/\" directory.$nc"
        fi
    fi

    if [ -z "$adapter_1" ]
    then
        adapter_1="AGATCGGAAGAGC"
        adapter_2="AGATCGGAAGAGC"
    fi

    if [ -z "$threads" ]
    then
        threads=16
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

###################################################################################################
usage() {
    echo "
        WGBS Pipeline

usage: path-to-pipe/WGBS.sh -1 <read_file1> [options]

Options: -1 <string>    The first input read filename (don't type io/<filename>)
         -2 <string>    The second input read file (only if paired-end)
         -a <string>    Sequence of an adapter ligated to the 3' end of the first read (default = AGATCGGAAGAGC)
         -b <string>    3' adapter to be removed from second read in a pair (default = AGATCGGAAGAGC)
         -t <int>       The number of threads (default = 16)
         -h             Print this help message
"
    exit
}

# Analysis Code
# Each step would assume the previous steps have been processed.
###################################################################################################
# Step 0, Preparation
s0_preprocess() {
    name=`echo ${file1%.fastq}`
    name=`echo ${name%.fq}`
    name=`echo ${name%_1}`
    name=`echo ${name%_R1}`

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

    source $pipe_path'/qc_source.sh'

    # start record
    echo "WGBS Pipeline" > QC_pipe_processing.log
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
# Step 1, Cutadapt
s1_cutadapt() {
    # TruSeq adapter
    echo -e "\n**************************************************"
    echo -e "${cyan}Step 1: Cutadapt - Trimming reads$nc"

    if [[ $types == SE ]]
    then
      $cutadapt -j $threads -a $adapter_1 -q 15,10 --minimum-length 36 -o $name'_trim_1.fq' /io/$file1 > $name'_cutadapt.trimlog' \
        && echo "step1.1, cutadapt with adapter $adapter_1 successful" >> QC_pipe_processing.log \
        || echo "step1.1, cutadapt with adapter $adapter_1 $adapter_2 fail......" >> QC_pipe_processing.log
    elif [[ $types == PE ]]
    then
      $cutadapt -j $threads -a $adapter_1 -A $adapter_2 -q 15,10 --minimum-length 36 -o $name'_trim_1.fq' -p $name'_trim_2.fq' /io/$file1 /io/$file2 > $name'_cutadapt.trimlog' \
        && echo "step1.1, cutadapt with adapter $adapter_1 $adapter_2 successful" >> QC_pipe_processing.log \
        || echo "step1.1, cutadapt with adapter $adapter_1 $adapter_2 fail......" >> QC_pipe_processing.log
    fi

    echo 'finish removing adapter'
    echo "**************************************************"
}

###################################################################################################
# Step 2, Bismark
s2_bismark(){
    head -40000000 $name'_trim_1.fq' > test1.fq
    if [[ $types == PE ]]
    then
      head -40000000 $name'_trim_2.fq' > test2.fq
    fi

    # echo -e "\n\n**************************************************"

    # 2.1
    echo 'start mapping reads to Lambda genome with bismark...'
    if [[ $types == SE ]]
    then
      $BISMARK_PATH/bismark --bowtie2 -p 4 --bam -B 'Lambda_genome_alignment' $LAMBDA_PATH test1.fq \
          && grep "Total number of C" Lambda_genome_alignment*.txt > lambda_C.txt \
          && grep "Total unmethylated C" Lambda_genome_alignment*.txt >> lambda_C.txt \
          && echo "step2.2, bowtie2 Lambda mapping successful" >> QC_pipe_processing.log \
          || echo "step2.2, bowtie2 Lambda mapping fail......" >> QC_pipe_processing.log
    elif [[ $types == PE ]]
    then
      $BISMARK_PATH/bismark --bowtie2 -p 4 --bam -B 'Lambda_genome_alignment' $LAMBDA_PATH -1 test1.fq -2 test2.fq \
          && grep "Total number of C" Lambda_genome_alignment*.txt > lambda_C.txt \
          && grep "Total unmethylated C" Lambda_genome_alignment*.txt >> lambda_C.txt \
          && echo "step2.2, bowtie2 Lambda mapping successful" >> QC_pipe_processing.log \
          || echo "step2.2, bowtie2 Lambda mapping fail......" >> QC_pipe_processing.log
    fi
    echo "**************************************************"

    # 2.3
    echo -e "\n\n**************************************************"
    echo 'start mapping reads to human genome with bismark...'
    if [[ $types == SE ]]
    then
      $BISMARK_PATH/bismark --bowtie2 --multicore 2 -p 2 --bam $GENOME_PATH $name'_trim_1.fq' \
          && echo "step2.3, bismark bowtie2 mapping successful" >> QC_pipe_processing.log \
          || echo "step2.3, bismark bowtie2 mapping fail......" >> QC_pipe_processing.log
    elif [[ $types == PE ]]
    then
      $BISMARK_PATH/bismark --bowtie2 --multicore 2 -p 2 --bam $GENOME_PATH -1 $name'_trim_1.fq' -2 $name'_trim_2.fq' \
          && echo "step2.3, bismark bowtie2 mapping successful" >> QC_pipe_processing.log \
          || echo "step2.3, bismark bowtie2 mapping fail......" >> QC_pipe_processing.log
    fi
    echo "**************************************************"
}

###################################################################################################
# Step 8, Sequence Lengths
s8_seq_lengths() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step 8: Count reads and calculate sequence lengths$nc\n"

    if [[ $types == SE ]]
    then
        /S/fqSeqLen -n $nameLen -s $seqLen -o $name'_seqLen.txt' $name'_trim_1.fq' \
            && echo -e "${green}Step 8: Count reads and calculate read lengths DONE$nc\n" 2>&1 | tee -a QC_pipe_processing.log \
            || { echo -e "${red}Step 8: Count reads and calculate read lengths FAIL$nc\n" 2>&1 | tee -a QC_pipe_processing.log && ok=false; }

    elif [[ $types == PE ]]
    then
        /S/fqSeqLen -n $nameLen -s $seqLen -o $name'_seqLen.txt' $name'_trim_1.fq' $name'_trim_2.fq' \
            && echo -e "${green}Step 8: Count reads and calculate read lengths DONE$nc" 2>&1 | tee -a QC_pipe_processing.log \
            || { echo -e "${red}Step 8: Count reads and calculate read lengths FAIL$nc" 2>&1 | tee -a QC_pipe_processing.log && ok=false; }
    fi

    echo "**************************************************"
}

###################################################################################################
# Step 9, Bowtie2
s9_bowtie2() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step 9: Bowtie2$nc"
    if [[ $types == SE ]]
    then
      bowtie2 -p 8 -x /D/hs38/Bisulfite_Genome/CT_conversion/BS_CT -U /io/$file1 -S $name'_CT_bt2.sam' --un unaligned.fq 2> $name'_ct_bt2.log' \
          && bowtie2 -p 8 -x /D/hs38/Bisulfite_Genome/GA_conversion/BS_GA -U unaligned.fq -S $name'_GA_bt2.sam' 2> $name'_ga_bt2.log' \
          && echo -e "${green}Step 9: Bowtie2 DONE$nc" 2>&1 | tee -a QC_pipe_processing.log \
          || { echo -e "${red}Step 9: Bowtie2 FAIL$nc" 2>&1 | tee -a QC_pipe_processing.log && ok=false; }
    elif [[ $types == PE ]]
    then
      bowtie2 -p 8 -x /D/hs38/Bisulfite_Genome/CT_conversion/BS_CT -1 /io/$file1 -2 /io/$file2 -S $name'_bt2.sam' 2> $name'_ct_bt2.log' \
          && bowtie2 -p 8 -x /D/hs38/Bisulfite_Genome/GA_conversion/BS_GA -U unaligned.fq -S $name'_bt2.sam' 2> $name'_ga_bt2.log' \
          && echo -e "${green}Step 9: Bowtie2 DONE$nc" 2>&1 | tee -a QC_pipe_processing.log \
          || { echo -e "${red}Step 9: Bowtie2 FAIL$nc" 2>&1 | tee -a QC_pipe_processing.log && ok=false; }
    fi
    echo "**************************************************"
}

###################################################################################################
s10_chr1() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step 10: Chromosome1$nc"
    samtools sort $name'_CT_bt2.sam' -o $name'_CT_sort.bam'
    samtools index $name'_CT_sort.bam'
    samtools view $name'_CT_sort.bam' chr01_CT_converted -o chr1_CT.sam

    samtools sort $name'_GA_bt2.sam' -o $name'_GA_sort.bam'
    samtools index $name'_GA_sort.bam'
    samtools view $name'_GA_sort.bam' chr01_GA_converted -o chr1_GA.sam
    echo "**************************************************"
}

###################################################################################################
# Step 11, MACS2
s11_macs2_peaks() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step 11: MACS2 Peak calling$nc"
    macs2 callpeak -t $name'_trim_1_bismark_bt2.bam' -f BAM -g hs -n $name'_macs2_B' -B -q 0.01 \
        && echo -e "${green}Step 11: MACS2 narrow peaks DONE$nc" 2>&1 | tee -a QC_pipe_processing.log \
        || { echo -e "${red}Step 11: MACS2 narror peaks FAIL$nc" 2>&1 | tee -a QC_pipe_processing.log && ok=false; }
    echo "**************************************************"
}


###################################################################################################
# Postprocess
s12_postprocess() {
    cd /io/'Processed_'$name/

    # Run FastQC and MultiQC
    fastqc --extract /io/$file1 --outdir=/io/'Processed_'$name/
    if [[ $types == PE ]]
    then
      fastqc --extract /io/$file2 --outdir=/io/'Processed_'$name/
    fi
    grep Deduplicated /io/'Processed_'$name/*/fastqc_data.txt > /io/'Processed_'$name/deduplicated.txt

    # multiqc .

    # Parsing results for QA/QC report.
    python3 /S/QcParser.py $name wgbs

    # End Time
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

    echo -e "\nWGBS Pipeline"
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

###################################################################################################
# run pipe ()
# step-by-step
clargs "$@"
s0_preprocess "$@"
s1_cutadapt
s2_bismark
s8_seq_lengths
s9_bowtie2
s10_chr1
s12_postprocess
