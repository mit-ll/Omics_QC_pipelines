
###############################################################################
# **This class represents a SAM record for a HTS sequence read**
# 
# **Author:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# This class represents a SAM record for a HTS sequence read.
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

#
#
###############################################################################
class SamRead ():

  #############################################################################
  #
  # SamRead constructor.
  #
  def __init__ ( self, samfile_line = "" ):
    self.dup = False
    self.parseLine( samfile_line )
    tokens = samfile_line.split()
    if ( len(tokens) >= 11 ):
      self.qname = tokens[0]
      self.flag = int(tokens[1])
      self.rname = tokens[2]
      self.pos = int(tokens[3])
      self.mapq = int(tokens[4])
      self.cigar = tokens[5]
      self.rnext = tokens[6]
      self.pnext = int(tokens[7])
      self.tlen = int(tokens[8])
      self.seq = tokens[9]
      self.qual = tokens[10]
      if ( ( self.flag & 1024) == 1024 ):
        self.dup = True
    else:
      self.qname = ""
      self.flag = 0
      self.rname = ""
      self.pos = 0
      self.mapq = 0
      self.cigar = ""
      self.rnext = ""
      self.pnext = 0
      self.tlen = 0
      self.seq = ""
      self.qual = ""

  #############################################################################
  #
  # Function to parse a SAM read line.
  #
  def parseLine( self, samfile_line = "" ):
    self.dup = False
    tokens = samfile_line.split()
    if ( len(tokens) >= 11 ):
      self.qname = tokens[0]
      self.flag = int(tokens[1])
      self.rname = tokens[2]
      self.pos = int(tokens[3])
      self.mapq = int(tokens[4])
      self.cigar = tokens[5]
      self.rnext = tokens[6]
      self.pnext = int(tokens[7])
      self.tlen = int(tokens[8])
      self.seq = tokens[9]
      self.qual = tokens[10]
      if ( ( self.flag & 1024) == 1024 ):
        self.dup = True
    else:
      self.qname = ""
      self.flag = 0
      self.rname = ""
      self.pos = 0
      self.mapq = 0
      self.cigar = ""
      self.rnext = ""
      self.pnext = 0
      self.tlen = 0
      self.seq = ""
      self.qual = ""

  #############################################################################
  #
  # Function to classify read as a PCR duplicate.
  #
  def isDuplicate( self, sam_read ):
    if ( self.seq == sam_read.seq and self.rname == sam_read.rname and self.pos == sam_read.pos ):
      self.dup = True
    return self.dup

###############################################################################
