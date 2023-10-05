
# **FASTQ iterator class**
# 
# **Authors:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# FASTQ iterator class
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


# This class provides an object model for accessing FASTQ sequences from
# a FASTQ sequence library input text file.
         
require 'input_file.rb'
require 'fastq_sequence.rb'

###############################################################################
class FastqIterator < InputFile

# Attributes:
# - end_of_file - the current end-of-file status (0 = false; 1 = true)
# - line - the current line of the input text file.
# - name - the name of the input text file.

# fastq sequences read from the current FASTQ file
attr_reader :fastqs


###############################################################################
# Create the FastqIterator object on named FASTQ sequence library file.
def initialize( name )
  super( name )		# initialize InputFile
  @fastq = nil
  @fastqs = {}
end  # method initialize


###############################################################################
# Return all of the FASTQ sequences in a hash keyed by sequence name.
def read_fastqs
  @fastqs = {}
  while ( self.is_end_of_file? == false )
    seq = self.next_sequence
    seq.sequence_data = SeqTools::trim_ns( seq.sequence_data ) if ! seq.nil?
    @fastqs[ seq.sequence_name ] = seq if (! seq.nil?) && (seq.sequence_data.length > 0)
  end  # while

  return @fastqs
end  # method read_fastqs


###############################################################################
# Return all of the FASTQ sequence names in a hash keyed by sequence name.
def read_names
  @names = {}
  while ( self.is_end_of_file? == false )
    seq = self.next_sequence
    @names[ seq.sequence_name ] = true if (! seq.nil?) && (seq.sequence_data.length > 0)
  end  # while

  return @names
end  # method read_names


###############################################################################
# Return all of the FASTQ sequences in to an array.
def fastqs_to_array
  @fastqs = []
  while ( self.is_end_of_file? == false )
    seq = self.next_sequence
    if (! seq.nil?) && (seq.sequence_data.length > 0)
      @fastqs << seq 
    end  # if
  end  # while

  return @fastqs
end  # method fastqs_to_array

###############################################################################
# Get the next sequence from the FASTQ library file.
def next_sequence
  # Create a new FASTQ sequence object.
  @fastq = FastqSequence.new( "", "", "", "", "" )

  # Parse the header line.
  self.next_line()
  @fastq.parse_header( @line.chomp )
  # puts "header: #{@line}"

  self.next_line()
  @fastq.sequence_data = @line.chomp()
  # puts "sequence: #{@line}"
 
  # Read in the '+' line. 
  self.next_line
  # puts "+ line: #{@line}"
  
  # Read in the quality data.
  self.next_line
  @fastq.sequence_quality = @line.chomp
  # puts "qual line: #{@line}"
  
  return @fastq
end  # method next_sequence

###############################################################################
# Get the next sequence from the FASTQ library file.
def next_pacbio_sequence
  # Create a new FASTQ sequence object.
  @fastq = FastqSequence.new( "", "", "", "", "" )

  # Parse the header line.
  @fastq.parse_header( @line.chomp )

  # puts "next_sequence: head: #{@line}"
  # Read in the sequence data.
  self.next_line
  # @fastq.sequence_data = @line.chomp
  seq = ""
  while ( ( self.is_end_of_file? == false ) &&
          ( @line[0] != '+' ) )
    seq = seq + @line.chomp()
    self.next_line()
  end  # while  
  @fastq.sequence_data = seq  
  # puts "seq : #{seq[0...20]}...#{seq[-20..-1]}"
  
  # puts "plus: #{@line}"
  
  # Read in the quality data.
  self.next_line
  @fastq.sequence_quality = @line.chomp
  qual = ""
  while ( ( self.is_end_of_file? == false ) &&
          ( qual.size < seq.size ) )
    qual = qual + @line.chomp()
    self.next_line()
  end  # while
  @fastq.sequence_quality = qual  
  # puts "qual: #{qual[0...20]}...#{qual[-20..-1]}"
  # puts "seq.size #{seq.size} : qual.size #{qual.size}"
  
  return @fastq
end  # method next_pacbio_sequence

###############################################################################

end  # class FastqIterator


###############################################################################
# Testing module.
def test_fastq_iterator
  in_fastq = FastqIterator.new( "test.fastq" )
  in_fastq.open_file
  in_fastq.next_line
  while ( in_fastq.is_end_of_file? == false )
    fastq = in_fastq.next_sequence
    if ( fastq != nil )
      print "name = ", fastq.sequence_name, " "
      print "desc = ", fastq.sequence_description, "\n"
      print "seq  = ", fastq.sequence_data, "\n"
      fastq.parse_quality( fastq.sequence_quality )
      print "qual = ", fastq.quality, "\n"
      puts
    end  # if
  end  # while
  in_fastq.close_file
end  # method test_fastq_iterator


###############################################################################
# test_fastq_iterator()

