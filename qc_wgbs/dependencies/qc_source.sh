#!/bin/bash

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


# get pipe path, though readlink/realpath can do it, some version doesn't have that
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# executives
adapter_1="CTGTCTCTTATACACATCT"
adapter_2="CTGTCTCTTATACACATCT"
GENOME_PATH="/D/hs38"
LAMBDA_PATH='/D/lambda'
G_size='/D/chrom_sizes.txt'
BISMARK_PATH='/S/miniconda/bin/'
PRESEQ_PATH='/S/miniconda/bin/'
find_id='/D/find_image_ID_digest.sh'
cutadapt='/usr/local/bin/cutadapt'

BOWTIE2=/S/miniconda/bin/bowtie2
SAMTOOLS=/S/miniconda/bin/samtools
FASTQC=/S/miniconda/bin/fastqc
methylQA=/S/methylQA/methylQA
macs2=/usr/local/bin/macs2

macs2_genome='hs'
