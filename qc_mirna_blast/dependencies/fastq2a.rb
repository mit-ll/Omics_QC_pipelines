
# **FASTQ to FASTA format program**
# 
# **Authors:**  Derek Ng, Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# FASTA to FASTQ format program
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


# This program converts a FASTQ file to a FASTA file.

require 'fastq_sequence.rb'         
require 'fastq_iterator.rb'
require 'output_file.rb'

###############################################################################
class Fastq2a

###############################################################################
def self.fastq2a( fastq_filename, fasta_filename )
  puts "fastq2a called; #{fastq_filename} and #{fasta_filename}"

  # Assert: fastq_filename & fasta_filename set.
  return if ( fasta_filename.nil? ) && ( fasta_filename.size < 1 )
  return if ( fastq_filename.nil? ) && ( fastq_filename.size < 1 )

  puts "Processing sequences:" 

  fasta_file = OutputFile.new( fasta_filename )
  fasta_file.open_file

  in_fastq = FastqIterator.new( fastq_filename )
  in_fastq.open_file
  # in_fastq.next_line
  while ( in_fastq.is_end_of_file? == false )
    fastq = in_fastq.next_sequence
    if ( fastq != nil )
      fasta_file.write( fastq.to_fasta )
    end  # if
  end  # while
  in_fastq.close_file
  fasta_file.close_file if ! fasta_file.nil?
end  # fastq2a
 

###############################################################################

end  # class Fastq2a


###############################################################################
# Main program.
def fastq2a_main( fastq_filename, fasta_filename )
  puts "usage: ruby fastq2a.rb <fastq_file> <fasta_file>\n" if ARGV.length < 1
  
  fastq_filename = ARGV[0] if ARGV.length >= 1
  fasta_filename = ARGV[1] if ARGV.length >= 2

  Fastq2a::fastq2a( fastq_filename, fasta_filename )
end  # method fastq2a_main


###############################################################################
fastq2a_main( "test.fastq", "test.fa" )

