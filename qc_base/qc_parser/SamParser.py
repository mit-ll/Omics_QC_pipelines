
###############################################################################
# **SAM file parser class**
# 
# **Author:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# SAM file parser class
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
import sys
from InputFile import InputFile;
from Peak import Peak;
from SamRead import SamRead;

###############################################################################
class SamParser ():

  #############################################################################
  #
  # SamParser constructor.
  #
  def __init__ ( self, name = "" ):
    self.name = name 

  #############################################################################
  #
  # SAM file parser.
  #
  def parseSam( self, sam_filename, peak_filename ):
    sam_file = InputFile()
    sam_file.setFileName( sam_filename )
    sam_file.openFile()

    peak_file = InputFile()
    peak_file.setFileName( peak_filename )
    peak_file.openFile()
    line = peak_file.nextLine()
    peak = Peak( line )

    duplicates = 0
    total = 0
    peak_total = 0
    peak_dups = 0
    read1 = SamRead( "" )
    while ( sam_file.isEndOfFile() == False ):
      line = sam_file.nextLine()
      sam_read = SamRead( line )
      if ( len(sam_read.seq) > 0 ):
        total += 1
        if ( sam_read.isDuplicate( read1 ) == True ):
          duplicates += 1

        # Advance Peak if past current peak
        while ( (sam_read.rname != peak.chrom or sam_read.pos > peak.chrom_end) and peak_file.isEndOfFile() == False ):
          line = peak_file.nextLine()
          peak.parseLine( line )

        # Check for read within a Peak
        if ( peak.overlaps( sam_read.rname, sam_read.pos ) ):
          peak_total += 1

          # Check for a duplicate read within a Peak
          if ( sam_read.dup == True ):
            peak_dups += 1

      read1 = sam_read

    sam_file.closeFile()
    peak_file.closeFile()

    print( "PCR duplicates: " + str(duplicates) )
    print( "# sequences: " + str(total) )
    print( "Reads in peaks: " + str(peak_total) )
    print( "Duplicate reads in peaks: " + str(peak_dups) )
    return duplicates, total, peak_total, peak_dups

  #############################################################################
  #
  #

###############################################################################
#
# This code tests the SamParser object.
#
# test = SamParser ()
# arg_count = len(sys.argv)
# print("# of arguments: " + str(arg_count) )
# print("argv[1]: " + str(sys.argv[1]) )
# if ( arg_count >= 3 ):
  # sam_name = sys.argv[1]
  # peak_name = sys.argv[2]
  # test.parseSam( sam_name, peak_name )
