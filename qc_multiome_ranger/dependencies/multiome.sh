#!/bin/bash
###################################################################################################
# **10X-Multipome epigenetics and transcriptomics QA/QC Singularity pipeline**
# 
# **Authors:**  Derek Ng, Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# 10X-Multipome epigenetics and transcriptomics QA/QC Singularity pipeline.
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
#          FILE: multiome.sh
#         USAGE: path-to-pipe/multiome.sh -1 <fastq1> [options]
#
#       OPTIONS: -1 <string>    The first input read filename (required option; don't type io/<filename>)
#                -d <string>    The project directory
#                -s <string>    The sample name
#                -t <int>       The number of threads (default = 1)
#                -h             Print this help message
#
#   DESCRIPTION: multiome QAQC pipeline
###################################################################################################


# Pipe Start
###################################################################################################
# Read all necessary parameters and prepare data structure

# multiome Pipeline

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

echo -e "\n\n\n${black}${cyanb}multiome Pipeline${nc}\n"
echo -e "${startDate}\n"

# get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ok=true

# Command line arguments
clargs() {
    while getopts a:d:r:s:t:h opts
    do
        case "$opts" in
            a) file_atac="$OPTARG";; # scATAC-seq fastq file
            d) nameDir="$OPTARG";; # Directory name
            r) file_gex="$OPTARG";; # scRNA-seq fastq file
            s) sample="$OPTARG";; #Sample name
            t) threads="$OPTARG";; # Number of threads
            h) usage;;
            *) exit;;
        esac
    done

## Defaults
    fa=false
    fr=false
    faexten=false
    frexten=false
    d=false
    s=false
    t=false
    nameLen=128
    seqLen=256
    io=false
    atac=false
    gex=false
    type="Gene Expression"

 ## Check directory name
    if [ -z "nameDir" ]
    then
        usage
    else
        echo -e "Your argument to the -d option is ${nameDir}."
        d=true
    fi

 ## Figure out assay type & check file extension
    if [[ $file_gex != *.fastq ]]
    then
      frexten=true
      echo -e "File r does not have the correct file extension, make sure it is fastq |$file_gex|"
    else
      gex=true
      type="Gene Expression"

      echo -e "File r is a GEX file with the correct fastq.gz format"
    fi

    if [[ $file_atac != *.fastq ]]
    then
      faexten=true
      echo -e "File a does not have the correct file extension, make sure it is fastq |$file_atac|"
    else
      atac=true
      type="Chromatin Accessibility"
      echo -e "File a is an ATAC file with the correct fastq.gz format: |${file_atac}|"
    fi

## Thread check
    if [ -z $threads ]
    then
        threads=1
    elif [ ${threads} -le 0 ]
    then
        echo -e "${red}Your argument to the -t option is ${threads}. However, it must be a positive integer greater than 0.$nc"
        t=true
    fi

## Sample check
  # if [ -z $sample ]
  # then
  #    s=true
  #    echo -e "${red}Your argument to the -s option is ${sample}. However, it must be a the sample name.$nc"
  # fi


  if $io || $faexten || $frexten || $fa || $fr || $t
  then
      exit
  fi
}


usage() {
    echo "
Program Name: multiome Pipeline

Usage: path-to-pipe/multiome.sh -1 <fastq1> [options]

Options: -a <string>    The input scATAC-seq FASTQ filename
	 -d <string>    Directory name
	 -r <string>    The input scRNA-seq FASTQ filename
         -t <int>       The number of threads (default = 1)
         -h             Print this help message
"
    exit
}

# Analysis Code
# Each step would assume the previous steps have been processed.
###################################################################################################
# Step 0a, Preparation
s_preprocess() {
    name=`echo ${file_gex%.fastq}`
    name=`echo ${name%.fq}`
    name=`echo ${name%_1}`
    name=`echo ${name%_R1}`
    iof="/io/Processed_"$nameDir

## Make output dir if it doesn't exist
    if [ ! -d /io/'Processed_'$nameDir/ ]
    then
      mkdir /io/'Processed_'$nameDir/
    fi
    cd /io/'Processed_'$nameDir/

## Start  record
    echo "multiome 10X Ranger Pipeline: " > $iof/QC_pipe_processing.log
    echo "Arguments: $@" 2>&1 | tee -a $iof/QC_pipe_processing.log
    echo "Target scRNA-seq file is $file_gex" 2>&1 | tee -a $iof/QC_pipe_processing.log
    echo "Target scATAC-seq file is $file_atac" 2>&1 | tee -a $iof/QC_pipe_processing.log
}

################################################################################################
## Step, Run cellranger

s_cellranger-arc() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}cellranger-arc:$nc"

        # /S/cellranger-arc-2.0.1/cellranger-arc count --id=$nameDir --reference=/S/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 --libraries=libraries.csv

        # cp /io/processed/$nameDir/outs/summary.csv . \
        # && echo -e "${green}CellRanger-arc Step: DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        # || { echo -e "${red}CellRanger-arc Step: FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    echo "**************************************************"
}

###################################################################################################
# Step, Unzip fastqs

s_fastq_processing() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}fastq processing:$nc"

    if [ -f "/io/fastqs_gex/$nameDir/$file_gex" ]; then
      echo "$file_gex exists"
    else
      gzip -dk /io/fastqs_gex/$nameDir/$file_gex.gz 
    fi

    if [ -f "/io/fastqs_atac/$nameDir/$file_atac" ]; then
      echo "$file_atac exists"
    else
      gzip -dk /io/fastqs_atac/$nameDir/$file_atac.gz 
    fi

    echo -e "${green}fastq processing: DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log 
    echo "**************************************************"
}

###################################################################################################
# Step, Sequence Lengths

s_seq_lengths() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Count reads and calculate sequence lengths$nc"


    /S/fqSeqLen -n $nameLen -s $seqLen -o 'atac_seqLen.txt' /io/fastqs_atac/$nameDir/$file_atac \
    && /S/fqSeqLen -n $nameLen -s $seqLen -o 'rna_seqLen.txt' /io/fastqs_gex/$nameDir/$file_gex \
        && echo -e "${green}Step: Count reads and calculate read lengths DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: Count reads and calculate read lengths FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }


    echo "**************************************************"
}

###################################################################################################
# Step, Postprocess

s_postprocess() {
    # cd /io/'Processed_'$nameDir/
    cd /io/
    find . -name "summary.csv" -print | grep $nameDir > $nameDir'_summary.list'
    find . -name "filtered_feature_bc_matrix" -print | grep $nameDir > $nameDir'_feature.list'

    # Parsing results for QA/QC report.
    python3 /S/QcParser.py $nameDir multiome

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

    echo -e "\nmultiome Pipeline"
    echo "Target scATAC-seq file is $file_atac" 2>&1 | tee -a $iof/QC_pipe_processing.log
    echo "Target scRNA-seq file is $file_gex" 2>&1 | tee -a $iof/QC_pipe_processing.log

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


###################################################################################################
## Run pipeline
clargs "$@"
s_preprocess "$@"
# s_cellranger-arc "$@"
s_fastq_processing "$@"
s_seq_lengths "$@"
s_postprocess "$@"
