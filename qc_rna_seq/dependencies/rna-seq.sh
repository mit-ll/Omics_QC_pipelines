#! /bin/bash
###################################################################################################
# **RNA-seq transcriptomics QA/QC Singularity pipeline**
# 
# **Authors:**  Derek Ng, Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# RNA-seq transcriptomics QA/QC Singularity pipeline
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
#          FILE: RNA-seq.sh
#         USAGE: path-to-pipe/RNA-seq.sh -1 <fastq1> [options]
#
#       OPTIONS: -1 <string>    The first input read filename (required option; don't type io/<filename>)
#                -2 <string>    The second input read filename (only if paired-end)
#                -d <string>    The name of the dataset directory
#                -n <int>       The length of the name of a fastq read; used to allocate array (default = 128)
#                -s <int>       The length of the sequence of a fastq read; used to allocate array (default = 256)
#                -t <int>       The number of threads (default = 16)
#                -h             Print this help message
#
#   DESCRIPTION: total RNA-seq QAQC pipeline
###################################################################################################


# Pipe Start
###################################################################################################
# Read all necessary parameters and prepare data structure

# total RNA-seq Pipeline

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

echo -e "\n\n\n${black}${cyanb}scRNA-seq Pipeline${nc}\n"
echo -e "${startDate}\n"

# get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ok=true

###################################################################################################
# Command line arguments
clargs() {
  while getopts 1:2:d:t:h opts
  do
    case "$opts" in
        1) file1="$OPTARG";; # First Fastq file
        2) file2="$OPTARG";; # Second Fastq file
	d) nameDir="$OPTARG";; # Directory name
        t) threads="$OPTARG";; # Number of threads
        h) usage;;
        *) exit;;
    esac
  done

# Defaults

  types="PE"
  f1exten=false
  f2exten=false
  f1=false
  f2=false
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
  if [ ! -f /io/fastqs/$nameDir/$file1 ]
  then
      f1=true
  fi
  if [ ! -f /io/fastqs/$nameDir/$file2 ]
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
          echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"io/fastqs/$nameDir\" directory.$nc"
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
          echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"io/fastqs/$nameDir/\" directory.$nc"
      fi
      if $f2
      then
          echo -e "${red}The file $file2 does not exist. Please check spelling and make sure that it is in the \"io/fastqs/$nameDir/\" directory.$nc"
      fi
  fi

if [ -z ${threads} ]; then
    threads=16
fi

if [ -z ${name} ]; then
    name="$(cut -d'.' -f1 <<< $file1 | rev | cut -d'/' -f1 | rev)"
fi
}

###################################################################################################
usage() {
    echo "
        RNA-seq Pipeline

usage: rna-seq.sh -1 <read_file1> [options]

Options: -1 <string>    The first input read filename (don't type io/<filename>)
         -2 <string>    The second input read file (only if paired-end)
	 -d <string>    Directory name
         -t <int>       The number of threads (default = 16)
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
  
    # if [ -z ${nameDir} ]; then nameDir=$name; fi
  
    # if [ -d /io/'Processed_'$nameDir/ ]
    # then
    #     while true
    #     do
    #         echo -e "${red}The directory /io/Processed_${name}/ exists. Do you wish to overwrite? (y/n)$nc"
    #         read yn
    #         case $yn in
    #             [Yy]*) rm -r /io/'Processed_'$nameDir/; break;;
    #             [Nn]*) echo "Please rename or move the existing directory. You may also rename the input file(s)."; exit;;
    #             *) echo "Please answer y or n.";;
    #         esac
    #     done
    # fi
    if [ ! -d /io/'Processed_'$nameDir/ ]
    then
      mkdir /io/'Processed_'$nameDir/
    fi
    cd /io/'Processed_'$nameDir/

    # start record
    echo "RNA-seq Pipeline" > QC_pipe_processing.log
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
# Step, Preparation
s_rna_seq() {

# Setup output directories
if ! [ -d featureCountsOutput ]; then
    mkdir featureCountsOutput
fi

if ! [ -d star_alignments ]; then
    mkdir star_alignments
fi


STAR_INDEX=/D/star_index
STAR_GTF=/D/GRCh38.gtf
if [ ${types} == 'PE' ]; then
    
    # Trim
    skewer /io/fastqs/$nameDir/$file1 /io/fastqs/$nameDir/$file2 -o ${name}

    # Align with STAR
    mkdir star_alignments/${name}
    STAR --runThreadN $threads --genomeDir ${STAR_INDEX} \
         --readFilesIn ${name}-trimmed-pair1.fastq \
                       ${name}-trimmed-pair2.fastq \
         --outFileNamePrefix star_alignments/${name}/ \
         --quantMode TranscriptomeSAM GeneCounts
    grep ENSG star_alignments/${name}/ReadsPerGene.out.tab | sort > gene_counts.txt
    grep "Number of input reads" star_alignments/${name}/Log.final.out > rna_info.txt
    grep "Uniquely mapped reads %" star_alignments/${name}/Log.final.out >> rna_info.txt
    grep "Average input read length" star_alignments/${name}/Log.final.out >> rna_info.txt

    # Feature counts on the the SAM files
    featureCounts -p -a ${STAR_GTF} -g gene_id \
                  -o featureCountsOutput/${name}_featureCountsOutput.out \
                  star_alignments/${name}/Aligned.out.sam
else
    # Trim
    skewer /io/fastqs/$nameDir/$file1 -o ${name}

    # Align with STAR
    mkdir star_alignments/${name}
    STAR --runThreadN $threads --genomeDir ${STAR_INDEX} \
         --readFilesIn ${name}-trimmed.fastq \
         --outFileNamePrefix star_alignments/${name}/ \
         --quantMode TranscriptomeSAM GeneCounts
    grep ENSG star_alignments/${name}/ReadsPerGene.out.tab | sort > gene_counts.txt
    grep "Number of input reads" star_alignments/${name}/Log.final.out > rna_info.txt
    grep "Uniquely mapped reads %" star_alignments/${name}/Log.final.out >> rna_info.txt
    grep "Average input read length" star_alignments/${name}/Log.final.out >> rna_info.txt

    # Feature counts on the the SAM files
    featureCounts -a ${STAR_GTF} -g gene_id \
                  -o featureCountsOutput/${name}_featureCountsOutput.out \
                  star_alignments/${name}/Aligned.out.sam

fi
}

###################################################################################################
# Step, Sequence Lengths
s_seq_lengths() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Count reads and calculate sequence lengths$nc\n"

    if [[ $types == SE ]]
    then
        /S/fqSeqLen -n $nameLen -s $seqLen -o $name'_seqLen.txt' /io/fastqs/$nameDir/$file1 \
            && echo -e "${green}Step: Count reads and calculate read lengths DONE$nc\n" 2>&1 | tee -a QC_pipe_processing.log \
            || { echo -e "${red}Step: Count reads and calculate read lengths FAIL$nc\n" 2>&1 | tee -a QC_pipe_processing.log && ok=false; }

    elif [[ $types == PE ]]
    then
        /S/fqSeqLen -n $nameLen -s $seqLen -o $name'_seqLen.txt' /io/fastqs/$nameDir/$file1 /io/fastqs/$nameDir/$file2 \
            && echo -e "${green}Step: Count reads and calculate read lengths DONE$nc" 2>&1 | tee -a QC_pipe_processing.log \
            || { echo -e "${red}Step: Count reads and calculate read lengths FAIL$nc" 2>&1 | tee -a QC_pipe_processing.log && ok=false; }
    fi

    echo "**************************************************"
}

###################################################################################################
# Step, Postprocess
s_postprocess() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Postprocess$nc"

    cd /io/'Processed_'$nameDir/

    echo "FASTQC processing"
    fastqc --extract /io/fastqs/$nameDir/$file1 --outdir=/io/'Processed_'$nameDir/
    if [[ $types == PE ]]
    then
      fastqc --extract /io/fastqs/$nameDir/$file2 --outdir=/io/'Processed_'$nameDir/
    fi
    grep Deduplicated /io/'Processed_'$nameDir/*/fastqc_data.txt > /io/'Processed_'$nameDir/deduplicated.txt
    grep "%GC" /io/'Processed_'$nameDir/*/fastqc_data.txt >> rna_info.txt

    # echo "MultiQC processing"
    # multiqc .

    echo "QcParser running"
    # Parsing results for QA/QC report.
    python3 /S/QcParser.py $name rna-seq
    echo "QcParser done"

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

    echo -e "\ntotal RNA-seq Pipeline"
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
s_seq_lengths
s_rna_seq
s_postprocess

