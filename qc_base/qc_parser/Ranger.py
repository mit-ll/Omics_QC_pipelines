
###############################################################################
#
#  Author:     Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
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
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
# This class represents a Ranger record for a HTS alignment.
#
###############################################################################

import string
import gzip
from InputFile import InputFile;

###############################################################################
class Ranger():

  #############################################################################
  #
  # Ranger constructor.
  #
  def __init__ ( self, dataset_name, pathname ):
    self.barcodes = {}
    self.genes = {}
    self.gene_ids = {}
    self.cell_total = {}
    self.cell_mt = {}
    self.mt20_count = 0
    self.percent_mito = 0
    self.reads_total = 0

    # pathname = "/io/processed/" + dataset_name + "/outs/filtered_feature_bc_matrix/"
    # pathname = "/io/processed/" + dataset_name + "/" + dataset_name + "-<suffix>/outs/filtered_feature_bc_matrix/"

    # Read in the barcodes
    i = 1
    with gzip.open( pathname + 'barcodes.tsv.gz', 'rb') as f:
      for line in f:
        # self.barcodes.append(line)
        self.barcodes[i] = line.decode().rstrip()
        i = i + 1
      print( "Barcodes count: " + str(len(self.barcodes)) )

    # Read in the Genes
    i = 1
    with gzip.open( pathname + 'features.tsv.gz', 'rb') as f:
      for line in f:
        tokens = line.decode().rstrip().split( '\t' )
        self.gene_ids[i] = tokens[0]
        self.genes[i] = tokens[1]
        i = i + 1
      print( "Genes count: " + str(len(self.genes)) )

    # Read in the Gene Expression matrix
    i = 1
    with gzip.open( pathname + 'matrix.mtx.gz', 'rb') as f:
      for line in f:
        if i > 3:
          tokens = line.decode().rstrip().split( ' ' )
          gene_id = int(tokens[0])
          barcode_id = int(tokens[1])
          count = int(tokens[2])
          gene = self.genes[gene_id]
          if barcode_id in self.cell_total:
            self.cell_total[barcode_id] = self.cell_total[barcode_id] + count
          else:
            self.cell_total[barcode_id] = count

          # Check for a mitochondrial gene.
          if gene.startswith( "MT-" ):
            if barcode_id in self.cell_mt:
              self.cell_mt[barcode_id] = self.cell_mt[barcode_id] + count
              # print( gene + ": " + tokens[2] )
            else:
              self.cell_mt[barcode_id] = count
        i = i + 1

    # Compute the percentage of cells with >= 20% mitochondrial counts
    self.mt20_count = 0.0
    self.reads_total = 0.0
    for id, count in self.cell_total.items():
      self.reads_total = self.reads_total + 1.0
      if id in self.cell_mt:
        cell_mito = self.cell_mt[id] * 1.0 / count * 1.0
        if cell_mito >= 0.2:
          self.mt20_count = self.mt20_count + 1.0
    self.percent_mito = self.mt20_count * 100.0 / self.reads_total
    print( "Mito cells: " + str(self.mt20_count) + " /  " + str(self.reads_total) + " = " + str(self.percent_mito) )

  #############################################################################
  #
  # Function to get the count of reads for this gene list.
  #
  def getPercentMito( self ):
    return self.percent_mito

  #############################################################################
  #
  # Function returns the barcode length.
  #
  def getBarcodeLength( self ):
    if 1 in self.barcodes:
      return len(self.barcodes[1])-2
    else:
      return 16		# default barcode length

# print( "Ranger running" )
# app = Ranger( "DU19-01S0003431_1" )

