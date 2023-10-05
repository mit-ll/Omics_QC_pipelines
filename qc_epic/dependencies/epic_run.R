
# **EPIC epigenetics QA/QC Singularity pipeline**
# 
# **Author:**  Philip Fremont-Smith, Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# EPIC epigenetics QA/QC Singularity pipeline
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


library(ggplot2)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#library(IlluminaHumanMethylationEPICmanifest)

#args <- "AS06-11984_5"
main <- function() {
  print(args)
  args <- commandArgs(trailingOnly = T)
  for (i in args) {
    print("in for loop")
    print(i)
    # dataset_loc <- "/io/rawdata"
    dataset_loc <- "/io"
    rgSet <- read.metharray.exp(file.path(dataset_loc,i))

    ## Get p-value fail percent ##
    fail <- as.data.frame(detectionP(rgSet))
    percent_fail <- (length(which(fail[,1] > 0.01)) / nrow(fail))*100

    ## beta_peaks plot, diff, mid peak##
    beta <- as.data.frame(getBeta(rgSet))
    colnames(beta) <- 'value'
    title <- i
    bet <- ggplot(beta, aes(x=value)) +
      geom_density() +
      ggtitle(i)

    clean <- as.data.frame(beta[!is.nan(beta$value),])
    colnames(clean) <- 'value'
    peak_1 <- max(density(clean$value)$y[density(clean$value)$x > 0.75])
    peak_2 <- max(density(clean$value)$y[density(clean$value)$x < 0.25])
    peak_3 <- max(density(clean$value)$y[density(clean$value)$x > 0.4 & density(clean$value)$x < .6])
    diff <- abs(peak_1 - peak_2)

    ## Get CpG Methylation ##
    gset <- mapToGenome(rgSet)
    gset <- addSnpInfo(gset)
    ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    keep <- !(featureNames(gset) %in% ann$Name[ann$chr %in% c("chrX","chrY")])
    gset <- gset[keep,]
    M <- getM(gset)
    cpg_per <- (length(which(M[,1] > 0)) / nrow(M))*100

    ## Write Output ##
    output <- paste(title, percent_fail, peak_1, peak_2, peak_3, diff, cpg_per, sep = ',')

    dir.create(paste0('/io/Processed_',i))
    base_dir <- paste0('/io/Processed_',i)
    summary_stats <- file.path(base_dir, paste0('summary_stats_',i,'.csv'))
    # print('Writing M_file_<i>.csv')
    # m_file <- file.path(base_dir, paste0('M_file_',i,'.csv'))
    print('Writing beta_<i>.png')
    beta_out <-file.path(base_dir, paste0('beta_',i,'.png'))

    write.table(output, summary_stats, quote = FALSE, col.names = FALSE, row.names = FALSE)
    # write.table(M, m_file, quote = FALSE, col.names = FALSE, row.names = TRUE)
    print('Tables written')
    ggsave(beta_out)
  }
}

main()
