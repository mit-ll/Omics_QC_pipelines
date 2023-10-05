

# **Mass CyTOF epigenetics QA/QC Singularity pipeline**
# 
# **Author:**  Philip Fremont-Smith (mailto: Philip.Fremont-Smith@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# ATAC-seq epigenetics QA/QC Singularity pipeline
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

# Get the command line arguments; filename passed as 2nd argument
args = commandArgs()

# Setup flowCore library; Read FCS file
library(flowCore)
samp <- read.FCS(args[2])
fcm <- as.data.frame((exprs(samp)))

# Count cells and viable cells; write to results.csv file
out <- file("results.csv", "w")
cat('"Cells",', nrow(fcm), "\n", file=out)
cat('"Viable",', nrow(fcm[fcm$Rh103Di<7.0, ]), "\n", file=out)
cat('"Pt195 1+",', nrow(fcm[fcm$Pt195Di>1.0, ]), "\n", file=out)
close(out)

# Cleanup
rm(samp, fcm, out)
q(save="no")
