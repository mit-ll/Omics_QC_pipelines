#!/bin/bash

#  Author:     Derek Ng, Philip Fremont-Smith, Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2022 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
#
# DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.
#
# This material is based upon work supported by the Defense Advanced Research 
# Projects Agency under Air Force Contract No. (FA8702- 15-D-0001). Any opinions, 
# findings and conclusions or recommendations expressed in this material are 
# those of the author(s) and do not necessarily reflect the views of the 
# Defense Advanced Research Projects Agency.
#
# Copyright © 2022 Massachusetts Institute of Technology.
#
# The software/firmware is provided to you on an As-Is basis
#
# Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS Part 252.227-7013 
# or 7014 (Feb 2014). Notwithstanding any copyright notice, U.S. Government rights in this 
# work are defined by DFARS 252.227-7013 or DFARS 252.227-7014 as detailed above. Use of this 
# work other than as specifically authorized by the U.S. Government may violate any copyrights 
# that exist in this work.
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

#===============================================================================
#          FILE: ATAC-seq.sh
#         USAGE: path-to-pipe/ATAC-seq.sh -1 <fastq1> [options]
#
#       OPTIONS: -1 <string>    The first input read filename (required option; don't type io/<filename>)
#                -2 <string>    The second input read filename (only if paired-end)
#                -n <int>       The length of the name of a fastq read; used to allocate array (default = 128)
#                -s <int>       The length of the sequence of a fastq read; used to allocate array (default = 256)
#                -t <int>       The number of threads (default = 16)
#                -m <string>    Marker (default = no_annotation)
#                -h             Print this help message
#
#   DESCRIPTION: ATAC-seq QAQC pipeline
#
#===============================================================================


# Pipe Start
###################################################################################################
# Read all necessary parameters and prepare data structure

# ATAC-seq Pipeline

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

echo -e "\n\n\n${black}${cyanb}ATAC-seq Pipeline${nc}\n"
echo -e "${startDate}\n"

pipe_version="v1.00"
host="mitll/atac-seq pipe"

# get the absolute path
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ok=true

# Command line arguments
clargs() {
    while getopts 1:2:d:n:s:t:m:c:h opts
    do
        case "$opts" in
            1) file1="$OPTARG";; # First Fastq file
            2) file2="$OPTARG";; # Second Fastq file
            d) nameDir="$OPTARG";;
            n) nameLen="$OPTARG";; # Length of name of fastq read
            s) seqLen="$OPTARG";; #length of sequence of fastq read
            t) threads="$OPTARG";; # Number of threads
            m) marker="$OPTARG";; # default 'no_annotation'
            h) usage;;
            *) exit;;
        esac
    done

    types="PE"
    f1exten=false
    f2exten=false
    f1=false
    f2=false
    d=false
    n=false
    s=false
    t=false
    m=false

## dir check
    if [! -z "nameDir" ]
    then
            echo -e "${red}Please retry with directory name in -d option"
            d=true
    fi

## file check
    if [[ $file1 != *.fastq ]] && [[ $file1 != *.fq ]]
    then
        f1exten=true
    fi
    if [[ $file2 != *.fastq ]] && [[ $file2 != *.fq ]]
    then
        f2exten=true
    fi

    if [ ! -f /io/$nameDir/$file1 ]
    then
            f1=true
    fi
    if [ ! -f /io/$nameDir/$file2 ]
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
            echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"io/$nameDir\" directory.$nc"
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
            echo -e "${red}The file $file1 does not exist. Please check spelling and make sure that it is in the \"io/$nameDir\" directory.$nc"
        fi

        if $f2
        then
            echo -e "${red}The file $file2 does not exist. Please check spelling and make sure that it is in the \"io/$nameDir\" directory.$nc"
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
        threads=16
    elif [ $threads -le 0 ]
    then
        echo -e "${red}Your argument to the -t option is ${threads}. However, it must be a positive integer greater than 0.$nc"
        t=true
    fi

    if [ -z "$marker" ]
    then
        marker='no_annotation'
    fi

    if $f1exten || $f1 || $n || $s || $t || $m || $d
    then
        exit
    fi
}

usage() {
    echo "
Program Name: ATAC-seq Pipeline

Usage: path-to-pipe/ATAC-seq.sh -1 <fastq1> [options]

Options: -1 <string>    The first input read filename (required option; don't type io/<filename>)
         -2 <string>    The second input read filename (only if paired-end)
         -d <string>    The directory of input read filenames within io
         -n <int>       The length of the name of a fastq read; used to allocate array (default = 128)
         -s <int>       The length of the sequence of a fastq read; used to allocate array (default = 256)
         -t <int>       The number of threads (default = 16)
         -m <string>    Marker (default = no_annotation)
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
    iof="/io/Processed_"$nameDir

    if [ ! -d /io/'Processed_'$nameDir/ ]
    then
      mkdir /io/'Processed_'$nameDir/
    fi

    cd /io/'Processed_'$nameDir/

    # start record
    echo "ATAC-seq Pipeline" > $iof/QC_pipe_processing.log
    echo "Arguments: $@" 2>&1 | tee -a $iof/QC_pipe_processing.log
    if [[ $types == SE ]]
    then
        echo "Target file is $file1" 2>&1 | tee -a $iof/QC_pipe_processing.log
    elif [[ $types == PE ]]
    then
        echo "Target files are $file1 & $file2" 2>&1 | tee -a $iof/QC_pipe_processing.log
    fi
    echo -e "Type of reads is ${types}\n" 2>&1 | tee -a $iof/QC_pipe_processing.log

    if [[ $types == PE ]]
    then
        raw1=/io/$nameDir/$file1
        raw2=/io/$nameDir/$file2
    elif [[ $types == SE ]]
    then
        raw1=/io/$nameDir/$file1
    fi
} # s_preprocess

###################################################################################################
# Step, Trim ATAC-seq adapters and QC on seq file
s_cutadapt() {
    echo -e "\n**************************************************"
    echo -e "${cyan}Step 1.1: Cutadapt - Trimming ATAC reads$nc"

        cutadapt -q 15 -m 10 -u 16 -a AGATCGGAAGAGCACACGTCTGAAC -j $threads -o $iof/'cut_'$file1 /io/$nameDir/$file1 \
      && echo -e "${green}Step: Cutadapt - Trimming DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
      || { echo -e "${red}Step: Cutadapt - Trimming FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    if [[ $types == PE ]]
    then
      cutadapt -q 15 -m 10 -u 16 -a AGATCGGAAGAGCGTCGTGTAGGGA -j $threads -o $iof/'cut_'$file2 /io/$nameDir/$file2 \
        && echo -e "${green}Step: Cutadapt - Trimming reads DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: Cutadapt - Trimming reads FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    fi

    echo "**************************************************"
}  # s_cutadapt


###################################################################################################
# Step, BWA
s_bwa() {
    echo -e "\n**************************************************"
    echo -e "${cyan}Step: BWA - Align reads$nc"

    if [[ $types == PE ]]
    then
        bwa mem -t $threads /D/GRCh38 /io/$nameDir/$file1 /io/$nameDir/$file2 > $iof/$name'_bwa.sam' \
            && echo -e "${green}Step: BWA - Align reads DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
            || { echo -e "${red}Step: BWA - Align reads FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    elif [[ $types == SE ]]
    then
         bwa mem -t $threads /D/GRCh38 /io/$nameDir/$file1 > $iof/$name'_bwa.sam' \
            && echo -e "${green}Step: BWA - Align reads DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
            || { echo -e "${red}Step: BWA - Align reads FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    fi

    echo "**************************************************"
} # s_bwa


###################################################################################################
# Step, Samtools sort
s_samtools_sort() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Samtools Sort - Convert from SAM to BAM$nc"

    /S/samtools-1.15.1/samtools sort -@ $threads -n $iof/$name'_bwa.sam' -o $iof/$name'_bwa.bam' \
        && echo -e "${green}Step: Samtools Sort - Convert from SAM to BAM DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: Samtools Sort - Convert from SAM to BAM FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    echo "**************************************************"
} # s_samtools_sort


###################################################################################################
# Step, Samtools markdup
s_samtools_markdup() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Samtools Markdup - Remove duplicates$nc"

    /S/samtools-1.15.1/samtools fixmate -@ $threads -m $iof/$name'_bwa.bam' $iof/$name'_bwa.fixout.bam' \
        && /S/samtools-1.15.1/samtools sort $iof/$name'_bwa.fixout.bam' -o $iof/$name'_bwa.sort.bam' -@ $threads \
        && /S/samtools-1.15.1/samtools index $iof/$name'_bwa.sort.bam' \
        && /S/samtools-1.15.1/samtools view -h -o $iof/$name'_chr1.sam' $iof/$name'_bwa.sort.bam' chr1 \
	&& /S/samtools-1.15.1/samtools flagstat $iof/$name'_bwa.sort.bam' -@ $threads > $iof/$name'_flagstat.txt' \
        && /S/samtools-1.15.1/samtools markdup $iof/$name'_bwa.sort.bam' $iof/$name'_bwa.rem.bam' -r -@ $threads -s > $iof/$name'_markdup.txt' \
	&& /S/samtools-1.15.1/samtools flagstat $iof/$name'_bwa.rem.bam' -@ $threads > $iof/$name'_flagstat2.txt' \
        && echo -e "${green}Step: Samtools Markdup - Remove duplicates DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: Samtools Markdup - Remove duplicates FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    echo "**************************************************"
} # s_samtools_markdup


###################################################################################################
# Step, Sequence Lengths
s_seq_lengths() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Count reads and calculate sequence lengths$nc"

    if [[ $types == SE ]]
    then
        /S/fqSeqLen -n $nameLen -s $seqLen -o $iof/$name'_seqLen.txt' /io/$nameDir/$file1 \
            && echo -e "${green}Step: Count reads and calculate read lengths DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
            || { echo -e "${red}Step: Count reads and calculate read lengths FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }

    elif [[ $types == PE ]]
    then
        /S/fqSeqLen -n $nameLen -s $seqLen -o $iof/$name'_seqLen.txt' /io/$nameDir/$file1 /io/$nameDir/$file2 \
            && echo -e "${green}Step: Count reads and calculate read lengths DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
            || { echo -e "${red}Step: Count reads and calculate read lengths FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    fi

    echo "**************************************************"
} # s_seq_lengths


###################################################################################################
# Step, Bowtie2
s_bowtie2() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Bowtie2$nc"
    if [[ $types == SE ]]
    then
      # bowtie2 -p 8 -x /D/grch38 -U /io/$nameDir/$file1 -S $name'_bt2.sam' --very-sensitive 2> $name'_bt2.log' \
      bowtie2 -p 8 -x /D/grch38 -U /io/$nameDir/$file1 -S $iof/$name'_bt2.sam' 2> $iof/$name'_bt2.log' \
          && echo -e "${green}Step: Bowtie2 DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
          || { echo -e "${red}Step: Bowtie2 FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    elif [[ $types == PE ]]
    then
      # bowtie2 -p 8 -x /D/grch38 -1 /io/$nameDir/$file1 -2 /io/$file2 -S $iof/$name'_bt2.sam' --very-sensitive 2> $iof/$name'_bt2.log' \
      bowtie2 -p 8 -x /D/grch38 -1 /io/$nameDir/$file1 -2 /io/$nameDir/$file2 -S $iof/$name'_bt2.sam' 2> $iof/$name'_bt2.log' \
          && echo -e "${green}Step: Bowtie2 DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
          || { echo -e "${red}Step: Bowtie2 FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    fi

    # /S/samtools-1.15.1/samtools sort -o $name'_bt2.bam' -T temp_aln $name'_presort.bam' \
    #     && /S/samtools-1.15.1/samtools sort -n -o $name'_trimmedDup.bam' $name'_presort.bam' \
    #     || good=false

    # if [[ $good == true ]]
    # then
    #     echo -e "${green}Step: Bowtie2 - Alignment DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log
    # else
    #     echo -e "${red}Step: Bowtie2 - Alignment FAILED$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log
    #     ok=false
    # fi

    echo "**************************************************"
} # s_bowtie2

###################################################################################################
# Step, Genrich peak calling
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
} # s_genrich_peaks


###################################################################################################
# Step, MACS2
s_macs2_peaks() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: MACS2 Peak calling$nc"
    macs2 callpeak -t $iof/$name'_bwa.sam' --broad -f SAM -g hs -n $iof/$name'_macs2_B' -B --broad-cutoff 0.1 \
        && echo -e "${green}Step: MACS2 broad peaks DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: MACS2 broad peaks FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    macs2 callpeak -t $iof/$name'_bwa.sam' --broad -f SAM -g hs -n $iof/$name'_macs2_N' -B -q 0.01 \
        && echo -e "${green}Step: MACS2 narrow peaks DONE$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log \
        || { echo -e "${red}Step: MACS2 narror peaks FAIL$nc" 2>&1 | tee -a $iof/QC_pipe_processing.log && ok=false; }
    echo "**************************************************"
} # s_macs2_peaks


###################################################################################################
s_tss() {
    echo -e "\n\n**************************************************"
    echo -e "${cyan}Step: Calculating TSS enrichment$nc"
    cp /S/getMat.py $iof
    cp /S/calTSSscore.r $iof
    /S/getCov1.sh $iof/$name'_bwa.bam' /D/hgRef.TSS.bed
    /S/getCov2.sh $iof/$name'_bwa.bam' /D/hgRef.TSS.bed
    /S/formatTSS.sh /D/hgRef.TSS.bed $iof/$name'_bwa.bam.sel.bam.gc'
    /S/calTSS.sh $iof/$name'_bwa.bam.sel.bam.gc.TSSmat.txt' 50 100 5
    echo "**************************************************"
}  # s_tss

###################################################################################################
# Step, Postprocess
s_postprocess() {
    cd /io/'Processed_'$nameDir/

    fastqc -t $threads --extract /io/$nameDir/$file1 --outdir=/io/'Processed_'$nameDir/
    if [[ $types == PE ]]
    then
      fastqc -t $threads --extract /io/$nameDir/$file2 --outdir=/io/'Processed_'$nameDir/
    fi
    grep Deduplicated $iof/*/fastqc_data.txt > $iof/deduplicated.txt

    # multiqc .

    # bedtools intersect -a $iof/$name'_bwa.sort.bam' -b /D/hgRef.TSS.bed > $iof/intersected.bam
    # bedtools genomecov -trackline -bg -ibam $iof/intersected.bam > $iof/intersected.bedgraph
    # perl /S/Histrader.pl --bedGraph $iof/intersected.bedgraph --peaks $iof/$name'.genrickPeak' --outBG > $iof/histrader.txt

    # tssenrich --genome hg19 --log log.txt --samtools-path /S/samtools-1.15.1/samtools --processes 16 $iof/$name'_bwa.sort.bam' > $iof/tss_score.txt

    # Parsing results for QA/QC report.
    python3 /S/frip.py $iof/$name'_bwa.sort.bam' $iof/$name'.genrichPeak' > $iof/frip_score.txt
    python3 /S/QcParser.py $name atac-seq

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

    echo -e "\nATAC-seq Pipeline"
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
# s_bowtie2
s_bwa
s_samtools_sort
s_samtools_markdup
# s_tss
s_seq_lengths
s_genrich_peaks
# s_macs2_peaks
s_postprocess
