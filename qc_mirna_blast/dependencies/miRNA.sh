#!/bin/bash
###################################################################################################
# **miRNA-seq transcriptomics QA/QC Singularity pipeline (BLAST version)**
# 
# **Authors:**  Derek Ng, Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# miRNA-seq transcriptomics QA/QC Singularity pipeline (BLAST version)
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
#          FILE: miRNA.sh
#         USAGE: path-to-pipe/miRNA.sh -1 <fastq> [options]
#
#       OPTIONS: -1 <string>    The first input read filename (required option; don't type io/<filename>)
#                -h             Print this help message
#
#   DESCRIPTION: miRNA QAQC pipeline
###################################################################################################


# Pipe Start
###################################################################################################
# Read all necessary parameters and prepare data structure

# miRNA Pipeline

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

echo -e "\n\n\n${black}${cyanb}miRNA Pipeline${nc}\n"
echo -e "${startDate}\n"

# get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ok=true

################################################################################
# Command line arguments
clargs() {
    while getopts 1:h opts
    do
        case "$opts" in
            1) file1="$OPTARG";; # Fastq file
            h) usage;;
            *) exit;;
        esac
    done

    nameLen=128
    seqLen=512
    types="SE"
    f1exten=false
    f1=false

    if [[ $file1 != *.fastq ]] && [[ $file1 != *.fq ]]
    then
        f1exten=true
    fi
    if [ ! -f /io/$file1 ]
    then
        f1=true
    fi

    if [ -z "$file1" ]
    then
        usage
    else
        types="SE"
        if $f1exten
        then
            echo -e "${red}The name of your file is $file1. However, the extension must be .fastq or .fq.$nc"
        fi
        if $f1
        then
            echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"io/\" directory.$nc"
        fi
    fi

    if $f1exten || $f1 
    then
        exit
    fi
}

################################################################################
usage() {
    echo "
Program Name: miRNA Pipeline

Usage: path-to-pipe/miRNA.sh -1 <fastq> 

Options: -1 <string>    The first input read filename (required option; don't type io/<filename>)
         -h             Print this help message
"
    exit
}

# Analysis Code
# Each step would assume the previous steps have been processed.
###################################################################################################
# Step, Preparation
s_preprocess() {
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
    if [ -d /io/'Processed_'$name/ ]
    then
      echo "/io/Processed_"$name" exists already"
    else
      mkdir /io/'Processed_'$name/
    fi
    cd /io/'Processed_'$name/

    # start record
    echo "miRNA Pipeline" > QC_pipe_processing.log
    echo "Arguments: $@" 2>&1 | tee -a QC_pipe_processing.log
    if [[ $types == SE ]]
    then
        echo "Target file is $file1" 2>&1 | tee -a QC_pipe_processing.log
    fi
    echo -e "Type of reads is ${types}\n" 2>&1 | tee -a QC_pipe_processing.log
}

################################################################################
# Step, miRNA process
s_mirna() {
    echo -e "\n**************************************************"
    echo -e "${cyan}Step: MicroRnas - Align miRNA reads$nc"

    if [[ $types == SE ]]
    then
       ruby /S/fastq2a.rb /io/$file1 $name'.fa' \
       && /S/blastn -query $name'.fa' -db /D/human.fa -out $name'.out' -word_size 15 \
       && java -Xmx128G -jar /S/BlastParser.jar $name'.out' > $name'.txt' \
       && wc -l $name'.txt' > 'miRNA_counts.txt' \
          && echo -e "${green}Step: miRNAs - Align reads DONE$nc" 2>&1 | tee -a QC_pipe_processing.log \
          || { echo -e "${red}Step: miRNAs - Align reads FAIL$nc" 2>&1 | tee -a QC_pipe_processing.log && ok=false; }
    fi

    echo -e "**************************************************"
}

###################################################################################################
# Step, Sequence Lengths
s_seq_lengths() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Count reads and calculate sequence lengths$nc"

    if [[ $types == SE ]]
    then
        echo -e "/S/fqSeqLen -n ${nameLen} -s ${seqLen} -o ${name}_seqLen.txt /io/${file1}"

        /S/fqSeqLen -n $nameLen -s $seqLen -o $name'_seqLen.txt' /io/$file1 
    fi

    echo -e "**************************************************"
}

################################################################################
# Step, Postprocess
s_postprocess() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Postprocess$nc"

    cd /io/'Processed_'$name/

    echo "FASTQC processing"
    fastqc --extract /io/$file1 --outdir=/io/'Processed_'$name/
    if [[ $types == PE ]]
    then
      fastqc --extract /io/$file2 --outdir=/io/'Processed_'$name/
    fi
    # grep Deduplicated /io/'Processed_'$name/*/fastqc_data.txt > /io/'Processed_'$name/deduplicated.txt

    # echo "MultiQC processing"
    # multiqc .

    # Parsing results for QA/QC report.
    python3 /S/QcParser.py $name miRNA

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

    echo -e "\nmiRNA Pipeline"
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
    # printf '%d.%d seconds\n\n' $S $millSec 2>&1 | tee -a QC_pipe_processing.log
}

# run pipe ()
###################################################################################################
# step-by-step
clargs "$@"
s_preprocess "$@"
s_mirna
s_seq_lengths
s_postprocess
