
###############################################################################
# **CpG counter class**
# 
# **Author:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# CpG counter class
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

import sys
import string
from InputFile import InputFile;
from MpileupBase import MpileupBase;

###############################################################################
class CgCounter ():

  #############################################################################
  #
  # Samtools Mpileup file parser.
  #
  def countCgs( self, mpileup_name ):
    mpileup_file = InputFile()
    mpileup_file.setFileName( mpileup_name )
    mpileup_file.openFile()

    cg_total = 0
    line = mpileup_file.nextLine()
    base1 = MpileupBase()
    base2 = MpileupBase()
    base1.parseLine( line )
    while ( mpileup_file.isEndOfFile() == False ):
      line = mpileup_file.nextLine()
      base2.parseLine( line )
      if ( ( base1.base == "C" ) and ( base2.base == "G" ) ):
        if ( base1.pos +1 == base2.pos ):
          cg_total += 1
      base1.parseLine( line )
    mpileup_file.closeFile()

    print( "CG bases: " + str(cg_total) )
    return cg_total

###############################################################################
#
# This code tests the CgCounter object.
#
# test = CgCounter ()
# arg_count = len(sys.argv)
# if ( arg_count >= 2 ):
#   mpileup_name = sys.argv[1]
#   test.countCgs( mpileup_name )
