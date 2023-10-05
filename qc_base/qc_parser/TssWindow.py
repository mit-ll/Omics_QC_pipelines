
###############################################################################
# **This class represents a TssWindow record ATAC-seq pipeine**
# 
# **Author:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# This class represents a TssWindow record ATAC-seq pipeine.
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
###############################################################################


import string

###############################################################################
class TssWindow ():

  #############################################################################
  #
  # TssWindow constructor.
  #
  def __init__ ( self, tss ):
    self.setTss( tss )

  #############################################################################
  #
  # TssWindow constructor.
  #
  def setTss( self, tss ):
    self.tss_pos = int((tss.chrom_end + tss.chrom_start) / 2)
    self.win_start = self.tss_pos - 2000
    self.win_end = self.tss_pos + 2000
    self.win_counts = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    self.core_total = 0
    self.max_count = 0
    # print("TssWindow ", self.tss_pos, " [", self.win_start, ", ", self.win_end, "]" )

  #############################################################################
  #
  # Function to add a peak to the  TssWindow read.
  #
  def addPeak( self, sam_read ):
    start_offset = int((sam_read.pos - self.tss_pos + 2000) / 100)
    end_offset = int((sam_read.pos + len(sam_read.seq) - self.tss_pos + 2000) / 100)

    if end_offset < 0 or end_offset >= 40:
      return

    # Increment the bin counts based on the start_offset
    if start_offset >= 0 and start_offset <= 39:
      self.win_counts[ start_offset ] += 1
      # print( "addPeak start", start_offset, sam_read.pos, sam_read.pos - self.tss_pos )

      if start_offset > 0 and end_offset < 39:
        self.core_total += 1

    # Increment the bin counts based on the end_offset if higher than start offset
    if end_offset > 0 and end_offset > start_offset & end_offset < 40:
      self.win_counts[ end_offset ] += 1
      # print( "addPeak end", end_offset, ", len: ", len(sam_read.seq) )

  #############################################################################
  #
  # Calculate the Transcription Enrichment for this promoter start site
  #
  def computeTss( self ):
    base_counts = (self.win_counts[0] + self.win_counts[39]) / 2.0
    self.max_counts = 0
    for i in range(1, 39):
      if self.win_counts[i] > self.max_counts:
        self.max_counts = self.win_counts[i]

    # print( "computeTss", self.max_counts, ", base: ", base_counts )
    enrichment = 0
    if base_counts > 0:
      enrichment = (self.max_counts * 1.0) / base_counts

    return enrichment

###############################################################################
