
###############################################################################
# **This class represents a Genes record for a HTS alignment**
# 
# **Author:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# This class represents a Genes record for a HTS alignment.
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
from InputFile import InputFile;

###############################################################################
class Genes():

  #############################################################################
  #
  # Genes constructor.
  #
  def __init__ ( self, filename ):
    self.genes = {}
    self.total = 0

    gene_file = InputFile()
    gene_file.setFileName( filename )
    gene_file.openFile()
    while ( gene_file.isEndOfFile () == 0 ):
      text = gene_file.nextLine ()
      if ( text != "" ):
        self.genes[ text ] = True
    gene_file.closeFile ()

  #############################################################################
  #
  # Function to check if a gene is in the gene list.
  #
  def checkGene( self, genename, count ):
    if genename in self.genes:
      self.total += count


  #############################################################################
  #
  # Function to get the count of reads for this gene list.
  #
  def getTotal( self ):
    return self.total


