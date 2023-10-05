
# **FASTQ sequence object model**
# 
# **Authors:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# FASTQ sequence object model
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


require 'fasta_sequence.rb'

# This class provides an object model for a FASTQ sequence (header line, sequence, quality).
         
###############################################################################
class FastqSequence < FastaSequence

# quality - ascii values of quality string
attr_accessor :quality
# quality_offset - starting zero value for quality encoding
attr_accessor :quality_offset
# sequence_quality - FASTQ quality string for this sequence
attr_accessor :sequence_quality

###############################################################################
# Create a new FASTA sequence object.
def initialize( name, desc, seq, seq_type, seq_qual )
  super( name, desc, seq, seq_type )
  @quality_offset = 33
  parse_quality( seq_qual )
end  # method initialize

###############################################################################
# This method parses the sequence quality string to integer values
def parse_quality( seq_qual )
  @sequence_quality = seq_qual
  @quality = []
  seq_qual.each_byte do |byte|
    # p=exp[ (ord(c) - 33) / ( -10 * log(10) )
    @quality << byte - @quality_offset
  end  # do
end  # method parse_quality

###############################################################################
# Generates a FASTA string for this sequence.
def to_fasta()
  return "" if @sequence_name.nil?

  str = ">" + @sequence_name + " " + @sequence_description + "\n"
  seq_start = 0
  seq_end = 50

  while ( seq_end < @sequence_data.length )
    str += @sequence_data[ seq_start...seq_end ] + "\n"
    seq_start = seq_end
    seq_end += 50
  end  # while

  str += @sequence_data[ seq_start...@sequence_data.length ] + "\n"
  return str 
end  # method to_string

###############################################################################
# Generates a FASTQ string for this sequence.
def to_string()
  str = "@" + @sequence_name + " " + @sequence_description + "\n" + 
      @sequence_data + "\n+\n" + @sequence_quality + "\n"
  return str 
end  # method to_string


###############################################################################

end  # FastqSequence


###############################################################################
# Testing module.
def test_fastq_sequence
  fastq = FastqSequence.new( "AB00001", "/name=\"fun\" /gene=\"test\"", "ACGTGTCATAGCAT", "DNA", "'-/#.!)(&$,*%)" )
  fastq.parse_header( ">AB00001 /name=\"fun\" /gene=\"test\" " )
  print "name = ", fastq.sequence_name, "\n"
  print "desc = ", fastq.sequence_description, "\n"
  fastq.parse_annotation
  fastq.parse_quality( fastq.sequence_quality )
  print "seq  = ", fastq.sequence_data, "\n"
  print "qual = ", fastq.quality, "\n"
end  # method test_fastq_sequence

# test_fastq_sequence()

