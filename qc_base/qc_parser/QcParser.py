
###############################################################################
# **This class provides QA/QC metrics scoring for epigenetics and transcriptomics results**
# 
# **Author:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
#  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
#  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
# 
# **RAMS request ID 1021178**
# 
# **Overview:**
# This class provides QA/QC metrics scoring for epigenetics and transcriptomics results.
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

import csv
import json
import os.path
import string
import sys
import numpy
from scipy.stats.stats import pearsonr

from CgCounter import CgCounter;
from FastaIterator import FastaIterator;
from FastaSequence import FastaSequence;
from Genes import Genes;
from InputFile import InputFile;
from Metrics import Metrics;
from MpileupBase import MpileupBase;
from OutputFile import OutputFile;
from Peak import Peak;
from Ranger import Ranger;
from SamParser import SamParser;
from SamRead import SamRead;
from TssWindow import TssWindow;

#
# This class parses QA/QC epigenetics and transcriptomics pipeline output results.
#
###############################################################################
class QcParser ():

  #############################################################################
  #
  # QcParser constructor.
  #
  def __init__ ( self, name = "" ):
    self.name = name 
    self.global_metrics = Metrics()
    self.ds_metrics = {}

  #############################################################################
  #
  # 
  #
  def firstInt( self, line ):
    tokens = line.split()
    return int(tokens[0])

  #############################################################################
  #
  # Function to return the last token
  #
  def lastToken( self, line ):
    tokens = line.split()
    return tokens[-1]

  #############################################################################
  #
  # Function to extract seqLength dataset name.
  #
  def seqLengthName( self, line ):
    tokens = line.split( "/" )
    return tokens[-1].strip( ".fastq" )


  #############################################################################
  #
  # samtools flagstat file parser.
  #
  def parseFlagstat( self, flagstat_name ):
    flagstat_file = InputFile()
    flagstat_file.setFileName( flagstat_name )
    flagstat_file.openFile()
    line = flagstat_file.nextLine()
    tokens = line.split()
    reads = float(tokens[0])
    self.global_metrics.setMetric( "Dataset size", reads, False )
    line = flagstat_file.nextLine()
    line = flagstat_file.nextLine()
    line = flagstat_file.nextLine()
    line = flagstat_file.nextLine()
    line = flagstat_file.nextLine()
    line = flagstat_file.nextLine()
    tokens = line.split()
    mapped = float(tokens[0])
    aligned = 0.0
    if reads > 0.0:
      aligned = (mapped * 100.0) / reads
    print( "Reads", reads, "Mapped", mapped, "Aligned %", aligned )
    return reads, aligned

  #############################################################################
  #
  # FRiP file parser.
  #
  def parseFrip( self, frip_name ):
    frip_file = InputFile()
    frip_file.setFileName( frip_name )
    frip_file.openFile()
    line = frip_file.nextLine()
    frip = 0.0
    if len(line) > 0:
      frip = float(line)
    return frip


  #############################################################################
  #
  # Transcript Start Site (TSS) enhancement file parser.
  #
  def parseTss( self, tss_name ):
    tss_file = InputFile()
    tss_file.setFileName( tss_name )
    tss_file.openFile()
    line = tss_file.nextLine()
    tss = 0.0
    if len(line) > 0:
      tss = float(line)
    return tss


  #############################################################################
  #
  # miRNA counts parser.
  #
  def parseMirnaCounts( self, counts_name ):
    in_file = InputFile()
    in_file.setFileName( counts_name )
    in_file.openFile()
    line = in_file.nextLine()
    tokens = line.split( " " )
    mirnas = int(tokens[0])
    in_file.closeFile()
    return mirnas


  #############################################################################
  #
  # seqLength file parser.
  #
  def parseSeqLength( self, seq_len_log, dataset_name, trim_length, r_name ):
    log_file = InputFile()
    filepath = "Processed_" + dataset_name + "/" + seq_len_log
    if not exists(filepath):
      return 0, 0
    log_file.setFileName( filepath )
    log_file.openFile()
    line = log_file.nextLine()
    line = log_file.nextLine()
    line = log_file.nextLine()
    ds1_name = self.seqLengthName( line )
    depth1 = log_file.nextLine()
    reads1 = int(self.lastToken(depth1))
    mean1 = log_file.nextLine()
    mean_len1 = float(self.lastToken(mean1))
    min_reads = reads1
    min_length = mean_len1
    self.ds_metrics[ r_name ] = Metrics()
    self.ds_metrics[ r_name ].setMetric( "Dataset size", reads1, False )
    self.ds_metrics[ r_name ].setMetric( "Mean trimmed sequence length", mean_len1-trim_length, False )
    # self.ds_metrics[ "R1 metrics" ].setMetric( "FastQC report", dataset_name + "_fastqc.html", True )

    line = log_file.nextLine()
    line = log_file.nextLine()
    if len(line) > 0:
      ds2_name = self.seqLengthName( line )
      depth2 = log_file.nextLine()
      reads2 = int(self.lastToken(depth2))
      mean2 = log_file.nextLine()
      mean_len2 = float(self.lastToken(mean2))
      self.ds_metrics[ "R2 metrics" ] = Metrics()
      self.ds_metrics[ "R2 metrics" ].setMetric( "Dataset size", reads2, False )
      self.ds_metrics[ "R2 metrics" ].setMetric( "Mean trimmed sequence length", mean_len2-trim_length, False )
      # self.ds_metrics[ "R2 metrics" ].setMetric( "FastQC report", ds2_name + "_fastqc.html", True )

      min_reads = min(reads1, reads2)
      min_length = min(mean_len1, mean_len2)
      self.global_metrics.setMetric( "Dataset size", min_reads, False )
      self.global_metrics.setMetric( "Mean trimmed sequence length", min_length-trim_length, False )

    print("Dataset size: " + str(min_reads))
    print("Mean trimmed sequence length: {0:.2f}".format(min_length-trim_length))
    return min_reads, min_length-trim_length


  #############################################################################
  def parseBismark( self, log_name ):
    log_file = InputFile()
    log_file.setFileName( log_name )
    log_file.openFile()

    # Number of alignments with a unique best hit from the different alignments:
    line = log_file.nextLine()
    tokens = line.split( ":" )
    aligned = int(tokens[1])

    # Mapping efficiency:
    line = log_file.nextLine()
    tokens = line.split( ":" )
    mapping_efficiency = float(tokens[1].replace('%',''))

    # C methylated in CpG context:
    line = log_file.nextLine()
    tokens = line.split( ":" )
    cg = float(tokens[1].replace('%',''))

    # C methylated in CHG context:
    line = log_file.nextLine()
    tokens = line.split( ":" )
    chg = float(tokens[1].replace('%',''))

    # C methylated in CHH context:
    line = log_file.nextLine()
    tokens = line.split( ":" )
    chh = float(tokens[1].replace('%',''))

    # CCC methylation rate
    line = log_file.nextLine()
    ccc = float(line)

    # genome coverage
    line = log_file.nextLine()
    tokens = line.split( "\t" )
    coverage = (1.0 - float(tokens[4])) * 100.0

    ch = (chg + chh) / 100.0
    return aligned, mapping_efficiency, cg, ch, ccc, coverage


  #############################################################################
  #
  # This reads the RNA-seq information file.
  #
  def readRnaInfo( self ):
    info_file = InputFile()
    info_file.setFileName( "rna_info.txt" )
    info_file.openFile()
    line = info_file.nextLine()
    tokens = line.split( "\t" )
    reads = int(tokens[1])
    line = info_file.nextLine()
    tokens = line.split( "\t" )
    uniquely_mapped = float(tokens[1].replace('%',''))
    line = info_file.nextLine()
    tokens = line.split( "\t" )
    ave_length = int(tokens[1])
    line = info_file.nextLine()
    tokens = line.split( "\t" ) 
    gc = int(tokens[1])
    info_file.closeFile()
    return reads, uniquely_mapped, ave_length, gc

  #############################################################################
  #
  # This function returns the Bismark feature count.
  #
  def getCount( self, line ):
    tokens = line.split( ":" )
    return float(tokens[1])

  #############################################################################
  #
  # Bismark Lambda Bisulfite conversion efficiency.
  #
  def bisulfite( self ):
    lambda_file = InputFile()
    lambda_file.setFileName( "lambda_C.txt" )
    lambda_file.openFile()
    line = lambda_file.nextLine()
    total_c = self.getCount( line )
    unmethylated_c = 0.0
    while lambda_file.isEndOfFile() == False:
      line = lambda_file.nextLine()
      if len(line) > 0:
        unmethylated_c += self.getCount( line ) 
    lambda_file.closeFile()
    return unmethylated_c * 100.0 / total_c


  #############################################################################
  #
  # BWA unique lines.
  #
  def parseUnique( self, unique_name ):
    unique_file = InputFile()
    unique_file.setFileName( unique_name )
    unique_file.openFile()
    line = unique_file.nextLine()
    unique_file.closeFile()
    unique_reads = int(line)
    return unique_reads


  #############################################################################
  #
  # FastQC deduplicated dataset name extractor.
  #
  def dupName( self, line ):
    tokens = line.split( "/" )
    if len(tokens) < 4:
      return ""
    return tokens[3].strip( "_fastqc" )

  #############################################################################
  #
  # FastQC deduplicated file parser.
  #
  def parseDuplicates( self, dups_log, dataset_name ):
    log_file = InputFile()
    log_file.setFileName( dups_log )
    log_file.openFile()
    line = log_file.nextLine()
    ds1_name = self.dupName( line )
    dedup1 = float(self.lastToken(line))
    line = log_file.nextLine()
    dup_rate = 100.0 - dedup1

    # Check for a second dataset FASTQ file.
    if len(line) > 0:
      ds2_name = self.dupName( line )
      dedup2 = float(self.lastToken(line))
      dup_rate2 = 100.0 - dedup2
     
      if len(ds2_name) > 0:
        if len(ds1_name) > 0:
          if (ds1_name in self.ds_metrics.keys()) == False:
            self.ds_metrics[ "R1 metrics" ] = Metrics()
          self.ds_metrics[ "R1 metrics" ].setMetric( "PCR duplicates rate", dup_rate, False )

        if (ds2_name in self.ds_metrics.keys()) == False:
          self.ds_metrics[ "R2 metrics" ] = Metrics()
        self.ds_metrics[ "R2 metrics" ].setMetric( "PCR duplicates rate", dup_rate2, False )

      dup_rate = 100.0 - min(dedup1, dedup2)

    self.global_metrics.setMetric( "PCR duplicates rate", dup_rate, False )

    log_file.closeFile()

    print("PCR duplicates rate: {0:.2f}%".format(dup_rate) )
    return dup_rate

  #############################################################################
  #
  # QcParser file parser.
  #
  def parseQc( self, bowtie2_log ):
    log_file = InputFile()
    log_file.setFileName ( bowtie2_log )
    log_file.openFile()
    reads = log_file.nextLine()
    paired = log_file.nextLine()
    unaligned = log_file.nextLine()
    once = log_file.nextLine()
    more = log_file.nextLine()
    log_file.closeFile()

    aligned_0 = self.firstInt( unaligned )
    aligned_1 = self.firstInt( once )
    aligned_n = self.firstInt( more )
    print("aligned  1: " + str(aligned_1) )
    print("aligned >1: " + str(aligned_n) )
    nrf = 0.0
    if aligned_1 + aligned_n > 0:
      nrf = float(aligned_1) / float(aligned_1 + aligned_n)
    print("NRF (non-redundant fraction): {0:.2f}".format(nrf) )
    unique = 0.0 
    aligned_percent = 0.0
    if aligned_0 + aligned_1 + aligned_n > 0:
      unique = float(aligned_1) / float(aligned_0 + aligned_1 + aligned_n)
      aligned_percent = (aligned_1 + aligned_n) * 100.0 / float(aligned_0 + aligned_1 + aligned_n)
    print("%uniquely_mapped: {0:.2f}%".format(unique*100.0) )
   
    return nrf, unique, aligned_1, aligned_percent


  #############################################################################
  #
  # Scoring function for Boolean.
  #
  def scoreBoolean( self, metric ):
    if metric:
      score = 1.0
    else:
      score = 0.0
    return score


  #############################################################################
  #
  # Scoring function for %GC.
  #
  def scoreGc( self, gc ):
    score = 0.0
    if gc < 20.0 or gc > 80.0:
      return score
    if gc >= 40.0 and gc <= 60.0:
      score = 2.0
    else:
      score = 1.0
    return score


  #############################################################################
  #
  # Scoring function for above thresholds.
  #
  def scoreHigh( self, metric, level1, level2 ):
    score = 0.0
    if ( metric+0.5 >= level1 ):  # round up for minimum score
      score = 1.0
      if ( metric >= level2 ):
        score = 2.0
    return score


  #############################################################################
  #
  # Scoring function for below thresholds.
  #
  def scoreLow( self, metric, level1, level2 ):
    score = 0.0
    if ( metric <= level1 ):
      score = 1.0
      if ( metric <= level2 ):
        score = 2.0
    return score


  #############################################################################
  #
  # JSON QC report writer.
  #
  def toJson( self, dataset_name ):
    json_file = OutputFile( dataset_name + ".json" )
    json_file.write( "{\n" )
    self.global_metrics.writeJson( json_file, "  ", "," )
    ds_names = list( self.ds_metrics.keys() )
    for i in range( len(ds_names) ):
      ds_name = ds_names[i]
      json_file.write( "  \"" + ds_name + " metrics\": {\n" )
      self.ds_metrics[ ds_name ].writeJson( json_file, "    ", "" )
      if i < len(ds_names) - 1:
        json_file.write( "  },\n" )
      else:
        json_file.write( "  }\n" )

    json_file.write( "}\n" )

  #############################################################################
  #
  # ATAC-seq QA/QC parser.
  #
  def atac_seq( self, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "ATAC-seq", True )

    trim_length = 0
    reads, length = self.parseSeqLength( ds_name + "_seqLen.txt", ds_name, trim_length, "R1" )
    flagstat, aligned = self.parseFlagstat( ds_name + "_flagstat.txt" )
    frip = self.parseFrip( "frip_score.txt" )

    sam_parser = SamParser()
    tss, nfr_peak = sam_parser.parseTss( ds_name + "_chr1.sam", "/D/refTss_chr1.bed" )

    depth_score = self.scoreHigh( reads, 25000000, 40000000 )
    aligned_score = self.scoreHigh( aligned, 50.0, 75.0 )
    flagstat_score = self.scoreHigh( flagstat, 15000000, 40000000 )
    frip_score = self.scoreHigh( frip, 0.05, 0.1 )
    trim_score = self.scoreHigh( length, 37.0, 50.0 )
    tss_score = self.scoreHigh( tss, 4.0, 6.0 )
    nfr_score = self.scoreBoolean( nfr_peak )
    mono_score = self.scoreBoolean( nfr_peak )

    min_score = min( aligned_score, depth_score, flagstat_score, frip_score, mono_score, nfr_score, trim_score, tss_score )
    ds_status = "not passed"
    ave_score = 0.0
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = ( aligned_score + depth_score + flagstat_score + frip_score + mono_score + nfr_score + trim_score + tss_score ) / 8.0

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Sequencing depth score", depth_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Aligned reads score", aligned_score, False )
    self.ds_metrics[ "QC" ].setMetric( "NoDup flagstat", flagstat_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Overlap FRiP", frip_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Mean trimmed sequence length score", trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "TSS enrichment score", tss_score, False )
    self.ds_metrics[ "QC" ].setMetric( "NFR peak score", nfr_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Monnuclear peak score", mono_score, False )

    self.global_metrics.setMetric( "Aligned reads", aligned, False )
    self.global_metrics.setMetric( "NoDuplicate flagstat", flagstat, False )
    self.global_metrics.setMetric( "Overlap FRiP", frip, False )
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )
    self.global_metrics.setMetric( "Sequencing depth", reads, False )
    self.global_metrics.setMetric( "TSS enrichment", tss, False )
    if nfr_peak:
      self.global_metrics.setMetric( "NFR peak", "Yes", True )
      self.global_metrics.setMetric( "Mononuclear peak", "Yes", True )
    else:
      self.global_metrics.setMetric( "NFR peak", "No", True )
      self.global_metrics.setMetric( "Mononuclear peak", "No", True )
    self.toJson( ds_name )


  #############################################################################
  #
  # ChIPmentation QA/QC parser.
  #
  def chipmentation( self, ds_name, antibody ):
    aligned_2 = 0
    nrf, unique, aligned_1 = self.parseQc( ds_name + "-R1_bt2.log" )
    if path.exists( ds_name + "-R2_bt2.log" ):
      nrf_r2, unique_r2, aligned_2 = self.parseQc( ds_name + "-R2_bt2.log" )
      nrf = max(nrf, nrf_r2)
      unique = max(unique, unique_r2)
    unique_aligned = aligned_1 + aligned_2
    min_reads, length, reads = self.parseSeqLength( ds_name + "_seqLen.txt" )
    unique_aligned = self.parseUnique( "unique.txt" )
    unique = float(unique_aligned) / float(reads)
    dup_rate = self.parseDuplicates( "deduplicated.txt", ds_name )

    unique_score = self.scoreHigh( unique, 0.6, 0.8 )
    if ( self.highLimit( antibody ) ):
      aligned_score = self.scoreHigh( unique_aligned, 20000000, 40000000 )
    else:
      aligned_score = self.scoreHigh( unique_aligned, 10000000, 20000000 )
    dup_score = self.scoreLow( dup_rate, 0.3, 0.15 )
    trim_score = self.scoreHigh( length, 50.0, 9999999.9 )
    nrf_score = self.scoreHigh( nrf, 0.5, 0.8 )
    min_score = min( aligned_score, unique_score, trim_score )
    ds_status = "not passed"
    ave_score = 0.0
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = statistics.mean( [aligned_score, unique_score, trim_score] )

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Uniquely Mapped Reads score", aligned_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Uniquely mapping efficiency score", unique_score, False )
    self.ds_metrics[ "QC" ].setMetric( "PCR duplicates rate score", dup_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Mean trimmed sequence length score", trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "NRF (non-redundant fraction) score", nrf_score, False )

    self.global_metrics.setMetric( "NRF (non-redundant fraction)", nrf, False )
    self.global_metrics.setMetric( "Uniquely Mapped Reads", unique_aligned, False )
    self.global_metrics.setMetric( "Uniquely mapped% (uniquely mapping efficiency)", unique*100.0, False )
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "ChIPmentation", True )
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )
    self.toJson( ds_name )


  #############################################################################
  #
  # Parse EPIC results file.
  #
  def epic( self, filename, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "EPIC Methylation", True )

    inf = InputFile()
    inf.setFileName( filename )
    inf.openFile()
    txt_line = inf.nextLine()
    inf.closeFile()

    tokens = txt_line.split( "," )
    percent_fail = float(tokens[1])
    peak_1 = float(tokens[2])
    peak_2 = float(tokens[3])
    peak_3 = float(tokens[4])
    diff = float(tokens[5])
    cpg = float(tokens[6])
    peaks_score = 1.0
    beta_value = "Fine"
    if peak_3 * 2.0 > peak_1 or peak_3 * 2.0 > peak_2:
      peaks_score = 0.0
      beta_value = "Multiple peaks"
    if diff * 2.0 > peak_1 or diff * 2.0 > peak_2:
      peaks_score = 0.0
      beta_value = "Peaks imbalanced"

    probe_score = self.scoreLow( percent_fail, 10.0, 1.0 )
    cpg_score = 0.0
    if cpg >= 20.0 and cpg <= 80.0:
      cpg_score = 1.0
    min_score = min( probe_score, peaks_score, cpg_score )
    ds_status = "not passed"
    ave_score = 0.0
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = numpy.mean( [probe_score, peaks_score, cpg_score] )

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "%failed probes score", probe_score, False )
    self.ds_metrics[ "QC" ].setMetric( "beta value distribution score", peaks_score, False )
    self.ds_metrics[ "QC" ].setMetric( "%CpG methylation score", cpg_score, False )

    self.global_metrics.setMetric( "% Failed probes", percent_fail, False )
    self.global_metrics.setMetric( "Beta value distribution", beta_value, True )
    self.global_metrics.setMetric( "%CpG methylation", cpg, False )
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )
    self.toJson( ds_name )

    return


  #############################################################################
  #
  # Parse massCyTOF results file.
  #
  def parse_cells( self, filename ):
    fcs = InputFile()
    fcs.setFileName ( filename )
    fcs.openFile()
    cells_line = fcs.nextLine()
    viable_line = fcs.nextLine()
    fcs.closeFile()

    tokens = cells_line.split( "," )
    cells = float( tokens[1] )
    tokens = viable_line.split( "," )
    viable = float( tokens[1] )
    return cells, viable

  #############################################################################
  #
  # massCyTOF QA/QC parser.
  #
  def mass_cytof( self, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "massCyTOF", True )

    cells = 0.0
    viable = 0.0
    cells, viable = self.parse_cells( "results.csv" )

    fraction_viable = viable / cells
    cell_score = self.scoreHigh( viable, 50000.0, 99999999.9 )
    viable_score = self.scoreHigh( fraction_viable, 0.6, 1.1 )

    min_score = min( cell_score, viable_score )
    ave_score = 0.0
    ds_status = "not passed"
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = numpy.mean( [cell_score, viable_score] )

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Cells score", cell_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Viable cells score", viable_score, False )

    self.global_metrics.setMetric( "Number of cells", cells, False )
    self.global_metrics.setMetric( "Viable cells", viable, False )
    self.global_metrics.setMetric( "Viable fraction%", fraction_viable*100.0, False )      
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )
    self.toJson( ds_name )


  #############################################################################
  #
  # Mint-ChIP QA/QC parser.
  #
  def mintchip( self, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "Mint-ChIP", True )

    nrf, unique, uniquely_mapped_reads_forward, aligned_percent, reads = self.parseQc( ds_name + "_bt2.log" )
    bad_count, length, bad_seqs = self.parseSeqLength( ds_name + "_seqLen.txt", ds_name, 0, "R1" )
    uniquely_mapped_reads = self.parseUnique( "unique.txt" )
    dedup_uniquely_mapped_reads = self.parseUnique( "dedup_unique.txt" )
    unique_reads_percent = (float(uniquely_mapped_reads) * 100.0) / float(reads)
    useful_reads_efficiency = (float(dedup_uniquely_mapped_reads) * 100.0) / float(uniquely_mapped_reads)
    dup_rate = self.parseDuplicates( "deduplicated.txt", ds_name )
    sam_parser = SamParser()
    sam_name = ds_name + "_bt2_sort.sam"
    peak_name = ds_name + "_macs2_B_peaks.broadPeak"
    duplicates, reads_total, peak_total, peak_dups = sam_parser.parseSam( sam_name, peak_name )
    reads_in_peaks = 0.0
    if reads_total > 0:
      reads_in_peaks = peak_total * 100.0 / (reads_total * 1.0)

    unique_score = self.scoreHigh( unique, 60.0, 85.0 )
    dup_score = self.scoreLow( dup_rate, 50.0, -0.15 )
    trim_score = self.scoreHigh( length, 50.0, 9999999.9 )

    mapped_score = self.scoreHigh( uniquely_mapped_reads_forward, 2000000, 3000000 )
    reads_score = self.scoreHigh( reads_in_peaks, 35.0, 50.0 )
    useful_score = self.scoreHigh( useful_reads_efficiency, 70.0, 999.9 )
    ds_status = "not passed"
    ave_score = 0.0

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Uniquely mapped% score", unique_score, False )
    self.ds_metrics[ "QC" ].setMetric( "PCR duplicates rate score", dup_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Mean trimmed sequence length score", trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Uniquely mapped reads score", mapped_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Percentage of reads within called peaks score", reads_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Useful reads efficiency score", useful_score, False )

    min_score = min( unique_score, trim_score, mapped_score, reads_score, useful_score )
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = statistics.mean( [unique_score, trim_score, mapped_score, reads_score, useful_score] )
    self.global_metrics.setMetric( "PCR duplicates rate", dup_rate, False )
    self.global_metrics.setMetric( "Percentage of reads within called peaks", reads_in_peaks, False )
    self.global_metrics.setMetric( "Uniquely mapped%", unique, False )
    self.global_metrics.setMetric( "Useful reads efficiency", useful_reads_efficiency, False )
    self.global_metrics.setMetric( "Uniquely mapped reads", uniquely_mapped_reads )
    self.global_metrics.setMetric( "Uniquely deduplicated mapped reads", dedup_uniquely_mapped_reads )
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )

    self.toJson( ds_name )


  #############################################################################
  #
  # Method to calculate the relative frequency of CG residues.
  #
  def relativeCg( self, c_total, g_total, cg_total, bp_total ):
    c_freq = c_total / bp_total
    g_freq = g_total / bp_total
    cg_freq = cg_total / bp_total
    cg_enrichment = cg_freq / c_freq * g_freq
    print( "CG enrichment", cg_enrichment, "CG freq", cg_freq, "CG total", cg_total, "Len:", bp_total, "C", c_total, "G", g_total)
    return cg_enrichment

  #############################################################################
  #
  # Method to count CpG site coverage in 100 base pair bins.
  #
  def count_chr1_cg( self ):
    chr1_file = FastaIterator( "/data/chr1" )
    chr1_file.openFile()

    # Initialize bins
    bins = []
    for i in range(101):
      bins.append(0)

    chr1 = chr1_file.nextSequence()
    chr1 = chr1_file.nextSequence()
    print("Name:", chr1.name)
    chr1_len = len(chr1.sequence)
    print("chr1 length", chr1_len)
    for i in xrange(0, chr1_len, 100):
      segment = chr1.sequence[i:i+100]
      cg_count = segment.count( "CG" )
      if ("NNNNNNNNNN" in segment) == False:
        bins[cg_count] += 1

    chr1_file.closeFile()

    print( "Bin\tCount" )
    for i in range(101):
      print(str(i)+"\t"+str(bins[i]))
    return bins

  #############################################################################
  #
  # Method to count CpG site coverage in peaks.
  #
  def count_cg( self, peak_name, ds_name ):
    genome_file = FastaIterator( "/data/GRCh38" )
    genome_file.openFile()

    peak_file = InputFile()
    peak_file.setFileName( peak_name )
    peak_file.openFile()
    line = peak_file.nextLine()
    peak = Peak( line )

    # Initialize bins
    bins = []
    for i in range(101):
      bins.append(0)

    peak_c_total = 0.0
    peak_g_total = 0.0
    peak_cg_total = 0.0
    peak_len_total = 0.0
    c_total = 0.0
    g_total = 0.0
    cg_total = 0.0
    len_total = 0.0

    # Read in the human genome sequences.
    while ( genome_file.isEndOfFile() == False ):
      fasta = genome_file.nextSequence()
      if fasta.sequence != "":
        c_total += fasta.sequence.count( "C" )
        g_total += fasta.sequence.count( "G" )
        cg_total += fasta.sequence.count( "CG" )
        len_total += len(fasta.sequence)

        # Advance Peak if past current genome sequence
        if ( peak.chrom != fasta.name and peak_file.isEndOfFile() == False ):
          line = peak_file.nextLine()
          peak.parseLine( line )

        # Bin chromosome 1 segments by 100 base pair bins.
        if ( peak.chrom == "chr1" ):
          chr1_len = len(fasta.sequence)
          for i in xrange(0, chr1_len, 100):
            bins[ fasta.sequence.count( "CG", i, i+100 ) ] += 1

        # Add peaks for this chromosome
        while ( peak.chrom == fasta.name and peak_file.isEndOfFile() == False ):
          peak_c_total += fasta.sequence.count( "C", peak.chrom_start-1, peak.chrom_end )
          peak_g_total += fasta.sequence.count( "G", peak.chrom_start-1, peak.chrom_end )
          peak_cg_total += fasta.sequence.count( "CG", peak.chrom_start-1, peak.chrom_end )
          peak_len_total += peak.chrom_end - peak.chrom_start + 1
          line = peak_file.nextLine()
          peak.parseLine( line )
          # print("Adding peak", peak.chrom, peak.chrom_start, peak.chrom_end)

    peak_file.closeFile()

    peak_ratio = peak_cg_total * 1.0 / (peak_len_total * 1.0)
    cg_ratio = (cg_total * 1.0) / (len_total * 1.0)
    print( "Peak ratio", peak_ratio, "genome CG ratio", cg_ratio)
    peak_enrichment = peak_ratio / cg_ratio
    print( "Peak CG enrichment", peak_enrichment )

    genome_file.closeFile()

    print( "Peak CG", peak_cg_total, "CG total", cg_total )
    cg_coverage = peak_cg_total * 1.0 / (cg_total * 1.0)
    peak_cg_enrichment = self.relativeCg( peak_c_total, peak_g_total, peak_cg_total, peak_len_total )
    genome_cg_enrichment = self.relativeCg( c_total, g_total, cg_total, len_total )
    cg_ratio = peak_cg_enrichment / genome_cg_enrichment
    print( "CG ratio", cg_ratio, "Peak CG", peak_cg_enrichment, "CG", genome_cg_enrichment)

    print( "Peak Bin\tCount" )
    for i in range(101):
      print(str(i)+"\t"+str(bins[i]))
    return cg_ratio, bins


  #############################################################################
  #
  # MeDIP-seq QA/QC parser.
  #
  def genome_cg( self ):
    genome_file = FastaIterator( "/D/GRCh38" )
    genome_file.openFile()

    # Read in the human genome sequences.
    cg_total = 0.0
    while ( genome_file.isEndOfFile() == False ):
      fasta = genome_file.nextSequence()
      if fasta.sequence != "":
        cg_total += fasta.sequence.count( "CG" )

    genome_file.closeFile()
    return cg_total


  #############################################################################
  #
  # MeDIP-seq QA/QC parser.
  #
  def medip( self, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "MeDIP", True )

    reads, length = self.parseSeqLength( ds_name + "_seqLen.txt", ds_name, 0, "R1" )

    mpileup_name = ds_name + "_mpileup.tsv"
    cg_counter = CgCounter()
    cg_found = cg_counter.countCgs( mpileup_name )
    cg_total = self.genome_cg()
    cg_coverage = 0.0
    if (cg_total > 0):
      cg_coverage = cg_found * 100.0 / cg_total

    trim_score = self.scoreHigh( length, 50.0, 9999999.9 )
    depth_score = self.scoreHigh( reads, 30000000, 999900000 )
    cg_score = self.scoreHigh( cg_coverage, 40.0, 60.0 )
    min_score = min( depth_score, trim_score, cg_score )
    ds_status = "not passed"
    ave_score = 0.0
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = numpy.mean( [depth_score, trim_score, cg_score] )

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Sequencing depth score", depth_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Mean trimmed sequence length score", trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "CpG site coverage score", cg_score, False )

    self.global_metrics.setMetric( "CpG site coverage", cg_coverage, False )
    self.global_metrics.setMetric( "Sequencing depth", reads, False )
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )

    self.toJson( ds_name )

  #############################################################################
  #
  # Function to count RNAs by gene lists.
  #
  def mirna( self, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "miRNA", True )

    reads, length = self.parseSeqLength( ds_name + "_seqLen.txt", ds_name, 0, "R1" )
    mirnas = self.parseMirnaCounts( "miRNA_counts.txt")

    trim_score = self.scoreHigh( length, 17.0, 9999999.9 )
    depth_score = self.scoreHigh( reads, 4000000, 8000000 )
    mirnas_score = self.scoreHigh( mirnas, 2000000, 98000000 )
    min_score = min( depth_score, mirnas_score, trim_score )
    ds_status = "not passed"
    ave_score = 0.0
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = numpy.mean( [depth_score, mirnas_score, trim_score] )

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Aligned Reads score", mirnas_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Mean trimmed sequence length score", trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Sequencing depth score", depth_score, False )

    self.global_metrics.setMetric( "Aligned Reads", mirnas, False )
    self.global_metrics.setMetric( "Sequencing depth", reads, False )
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )

    self.toJson( ds_name )

  #############################################################################
  #
  # Read 10X-CellRanger multiomics summary file.
  #
  def multiome_summary( self, pathname ):
    with open(pathname) as csvfile:
      reader = csv.DictReader(csvfile)
      for row in reader:
        gex_reads = float(row['GEX Sequenced read pairs'].replace(',', ''))
        cells = float(row['Estimated number of cells'].replace(',', ''))
        umi_count = float(row['GEX Median UMI counts per cell'].replace(',', ''))
        exons = float(row['GEX Reads mapped confidently to exonic regions'].replace(',', ''))
        cell_reads = float(row['GEX Fraction of transcriptomic reads in cells'].replace(',', ''))

        atac_reads = float(row['ATAC Sequenced read pairs'].replace(',', ''))
        frac_mapped = float(row['ATAC Confidently mapped read pairs'].replace(',', ''))
        frac_cut_peaks = float(row['ATAC Fraction of high-quality fragments overlapping peaks'].replace(',', ''))

        return gex_reads, cells, cell_reads, umi_count, exons, atac_reads, frac_mapped, frac_cut_peaks

  #############################################################################
  #
  # Read 10X-CellRanger multiomics metrics file.
  #
  def multiome_10x( self, dataset_name ):
    mito_percent = 0
    barcode_length = 0
    fealist = InputFile()
    fealist.setFileName( dataset_name + "_feature.list" )
    fealist.openFile()
    while fealist.isEndOfFile() == False:
      line = fealist.nextLine()
      if (len(line) > 0) and (line.find('-aggr') == -1):
        ranger = Ranger( dataset_name, line + "/" )
        mito_per = ranger.getPercentMito()
        mito_percent = max(mito_per, mito_percent)
        barcode_length = ranger.getBarcodeLength()

    # Initialize counters
    gex_reads = 0
    cells = 0
    cell_reads = 0
    umi_count = 0
    exons = 0
    atac_reads = 0
    frac_mapped = 0
    frac_cut_peaks = 0

    sumlist = InputFile()
    sumlist.setFileName( dataset_name + "_summary.list" )
    sumlist.openFile()
    while sumlist.isEndOfFile() == False:
      line = sumlist.nextLine()
      if (len(line) > 0) and (line.find('-aggr') == -1):
        gex_reads_l, cells_l, cell_reads_l, umi_count_l, exons_l, atac_reads_l, frac_mapped_l, frac_cut_peaks_l = self.multiome_summary( line )
        gex_reads += gex_reads_l
        cells += cells_l
        cell_reads = max(cell_reads_l, cell_reads)
        umi_count = max(umi_count_l, umi_count)
        exons = max(exons_l, exons)
        atac_reads += atac_reads_l
        frac_mapped = max(frac_mapped_l, frac_mapped)
        frac_cut_peaks = max(frac_cut_peaks_l, frac_cut_peaks)

    return gex_reads, cells, cell_reads, umi_count, exons, atac_reads, frac_mapped, frac_cut_peaks, mito_percent


  #############################################################################
  #
  # CellRanger multiome QA/QC parser.
  #
  def multiome( self, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "Multiome", True )

    # Parse scRNA-seq multiome results
    gex_reads, cells, cell_reads, umi_count, exons, atac_reads, frac_mapped, frac_cut_peaks, mito_percent = self.multiome_10x( ds_name )

    cells_score = self.scoreHigh( cells, 500.0, 1000.0 )
    exon_score = self.scoreHigh( exons * 100.0, 40.0, 60.0 )
    frac_cells_score = self.scoreHigh( cell_reads * 100.0, 40.0, 70.0 )
    mito_score = self.scoreLow( mito_percent, 50.0, 20.0 )
    umi_score = self.scoreHigh( umi_count, 500.0, 1500.0 )

    trim_rna_length = 0
    rna_reads, rna_length = self.parseSeqLength( "rna_seqLen.txt", ds_name, trim_rna_length, "GEX" )
    rna_trim_score = self.scoreHigh( rna_length * 2.0, 50.0, 9999999.9 )

    mapped_score = self.scoreHigh( frac_mapped * 100.0, 50.0, 75.0 )
    peaks_score = self.scoreHigh( frac_cut_peaks * 100.0, 5.0, 10.0 )

    trim_atac_length = 0
    atac_reads2, atac_length = self.parseSeqLength( "atac_seqLen.txt", ds_name, trim_atac_length, "ATAC" )
    atac_trim_score = self.scoreHigh( atac_length * 2.0, 37.0, 9999999.9 )

    atac_min_score = min( atac_trim_score, frac_cells_score, peaks_score )
    gex_min_score = min( cells_score, exon_score, mapped_score, mito_score, rna_trim_score, umi_score )
    atac_ds_status = "not passed"
    gex_ds_status = "not passed"
    ave_score = 0.0
    if ( gex_min_score > 0.0 ):
      gex_ds_status = "pass pending replication"
      ave_score = ( atac_trim_score + cells_score + exon_score + frac_cells_score + mapped_score + mito_score + peaks_score + rna_trim_score + umi_score ) / 9.0

    if ( atac_min_score > 0.0 ):
      atac_ds_status = "pass pending replication"
      ave_score = ( atac_trim_score + cells_score + exon_score + frac_cells_score + mapped_score + mito_score + peaks_score + rna_trim_score + umi_score ) / 9.0

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Estimated number of cells score", cells_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Reads Mapped to Exons score", exon_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Fraction Reads in Cells score", frac_cells_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Median UMI counts per cell score", umi_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Cells with >= 20% mitochondrial reads % score", mito_score, False )
    self.ds_metrics[ "QC" ].setMetric( "GEX mean trimmed sequence length score", rna_trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "ATAC mean trimmed sequence length score", atac_trim_score, False )

    self.global_metrics.setMetric( "ATAC Sequencing reads", atac_reads, False )
    self.global_metrics.setMetric( "GEX Sequencing reads", gex_reads, False )
    self.global_metrics.setMetric( "Estimated number of cells", cells, False )
    self.global_metrics.setMetric( "Median UMI counts per cell", umi_count, False )
    self.global_metrics.setMetric( "Reads Mapped to Exons", exons * 100.0, False )
    self.global_metrics.setMetric( "Fraction Reads in Cells", cell_reads * 100.0, False )
    self.global_metrics.setMetric( "Mitochondrial reads %", mito_percent, False )

    self.ds_metrics[ "QC" ].setMetric( "Mean scATAC-seq trimmed sequence length score", atac_trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "frac_mapped_confidently score", mapped_score, False )
    self.ds_metrics[ "QC" ].setMetric( "frac_cut_fragments_in_peaks score", peaks_score, False )

    self.global_metrics.setMetric( "frac_mapped_confidently", frac_mapped * 100.0, False )
    self.global_metrics.setMetric( "frac_cut_fragments_in_peaks", frac_cut_peaks * 100.0, False )
    self.global_metrics.setMetric( "Sequencing depth", atac_reads, False )

    self.global_metrics.setMetric( "ATAC QC status", atac_ds_status, True )
    self.global_metrics.setMetric( "GEX QC status", gex_ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )

    self.toJson( ds_name )


  #############################################################################
  #
  # Function to count RNAs by gene lists.
  #
  def countGenes( self ):
    globins = Genes( "/data/globin.ids" )
    rrnas = Genes( "/data/rRNA.ids" )
    mrnas = Genes( "/data/mRNA.ids" )
    mitos = Genes( "/data/mito.ids" )
    total = 0

    gene_file = InputFile()
    gene_file.setFileName( "gene_counts.txt" )
    gene_file.openFile()
    while ( gene_file.isEndOfFile () == 0 ):
      text = gene_file.nextLine ()
      if ( text != "" ):
        tokens = text.split( "\t" )
        genename = tokens[0]
        count = int(tokens[1])
        total += count
        globins.checkGene( genename, count )
        rrnas.checkGene( genename, count )
        mrnas.checkGene( genename, count )
        mitos.checkGene( genename, count )

    gene_file.closeFile ()
    return globins.getTotal(), rrnas.getTotal(), mrnas.getTotal(), mitos.getTotal(), total

  #############################################################################
  #
  # RNA-seq QA/QC parser.
  #
  def rna_seq( self, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "RNA-seq", True )

    reads, length = self.parseSeqLength( ds_name + "_seqLen.txt", ds_name, 0, "R1" )

    reads, uniquely_mapped, ave_length, gc = self.readRnaInfo()

    globins, rrnas, mrnas, mitos, total = self.countGenes()
    print("globin counts", globins)
    print("rRNA counts", rrnas)
    print("mRNA counts", mrnas)
    print("mito counts", mitos )
    print("Total counts", total)

    depth_score = self.scoreHigh( reads, 20000000, 50000000)
    gc_score = self.scoreGc( gc )
    globin_score = 0.0
    rrna_score = 0.0
    if ( reads > 0 ):
      globin_score = self.scoreLow( globins*1.0/reads, 0.3, 0.15 )
      rrna_score = self.scoreLow( rrnas*1.0/reads, 0.2, 0.1 )
    mrna_score = self.scoreHigh( mrnas, 16000000, 30000000 )
    trim_score = self.scoreHigh( ave_length, 50.0, 9999999.9 )
    unique_score = self.scoreHigh( uniquely_mapped, 40.0, 75.0 )

    min_score = min( depth_score, gc_score, globin_score, mrna_score, rrna_score, trim_score, unique_score )
    ds_status = "not passed"
    ave_score = 0.0
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = numpy.mean( [depth_score, gc_score, globin_score, mrna_score, rrna_score, trim_score, unique_score ] )

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Sequencing depth score", depth_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Mean trimmed sequence length score", trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "%GC score", gc_score, False )
    self.ds_metrics[ "QC" ].setMetric( "%globin score", globin_score, False )
    self.ds_metrics[ "QC" ].setMetric( "%rRNA score", rrna_score, False )
    self.ds_metrics[ "QC" ].setMetric( "mRNA reads score", mrna_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Uniquely mapped% score", unique_score, False )

    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )
    self.global_metrics.setMetric( "Sequencing depth", reads, False )
    self.global_metrics.setMetric( "%GC", gc, False )
    if ( reads > 0 ):
      self.global_metrics.setMetric( "%globin", globins*100.0/reads, False )
      self.global_metrics.setMetric( "%rRNA", rrnas*100.0/reads, False )
    self.global_metrics.setMetric( "mRNA reads", mrnas, False )
    self.global_metrics.setMetric( "Uniquely mapped%", uniquely_mapped, False )
    self.global_metrics.setMetric( "Mean trimmed sequence length", ave_length, False )
    self.toJson( ds_name )

  #############################################################################
  #
  # JSON float value.
  #
  def json_float( self, json_val ):
    if json_val == None:
      return -1.0
    return float(json_val)

  #############################################################################
  #
  # Read 10X-Genomics scATAC-seq metrics file.
  #
  def atac_10x( self ):
    with open('summary.json') as json_file:
      data = json.load(json_file)
      saturation = -1.0
      if 'bulk_estimated_saturation' in data:
        saturation = self.json_float(data['bulk_estimated_saturation'])
      frac_cut_peaks = self.json_float(data['frac_cut_fragments_in_peaks'])
      frac_targets = self.json_float(data['frac_fragments_overlapping_targets'])
      frac_mapped = self.json_float(data['frac_mapped_confidently'])
      return saturation, frac_cut_peaks, frac_targets, frac_mapped

  #############################################################################
  #
  # scATAC-seq QA/QC parser.
  #
  def sc_atac_seq( self, ds_name, len_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "scATAC-seq", True )

    trim_length = 0
    reads, length = self.parseSeqLength( len_name + "_seqLen.txt", ds_name, trim_length )
    saturation, frac_cut_peaks, frac_targets, frac_mapped, reads = self.atac_10x()

    trim_score = self.scoreHigh( length, 37.0, 9999999.9 )
    mapped_score = self.scoreHigh( frac_mapped, 50.0, 75.0 )
    peaks_score = self.scoreHigh( frac_cut_peaks, 5.0, 10.0 )
    min_score = min( mapped_score, peaks_score, trim_score )
    ds_status = "not passed"
    ave_score = 0.0
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = numpy.mean( [mapped_score, peaks_score, trim_score ] )

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Mean trimmed sequence length score", trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "frac_mapped_confidently score", mapped_score, False )
    self.ds_metrics[ "QC" ].setMetric( "frac_cut_fragments_in_peaks score", peaks_score, False )

    self.global_metrics.setMetric( "frac_mapped_confidently", frac_mapped, False )
    self.global_metrics.setMetric( "frac_cut_fragments_in_peaks", frac_cut_peaks, False )
    self.global_metrics.setMetric( "Sequencing depth", reads, False )
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )
    if saturation >= 0.0:
      self.global_metrics.setMetric( "Sequencing saturation", saturation, False )

    self.toJson( ds_name )


  #############################################################################
  #
  # Read 10X-Genomics metrics file.
  #
  def rna_10x( self, dataset_name ):
    ranger = Ranger( dataset_name )
    mito_percent = ranger.getPercentMito()
    barcode_length = ranger.getBarcodeLength()

    # pathname = "io/processed/" + dataset_name + "/outs/metrics_summary.csv"
    # with open( pathname ) as csvfile:
    with open('metrics_summary.csv') as csvfile:
      reader = csv.DictReader(csvfile)
      for row in reader:
        cells = float(row['Estimated Number of Cells'].replace(',', ''))
        cell_reads = float(row['Fraction Reads in Cells'].replace('%', ''))
        umi_count = float(row['Median UMI Counts per Cell'].replace(',', ''))
        exons = float(row['Reads Mapped Confidently to Exonic Regions'].replace('%', ''))
        return cells, cell_reads, umi_count, exons, mito_percent, barcode_length


  #############################################################################
  #
  # scRNA-seq QA/QC parser.
  #
  def sc_rna_seq( self, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "scRNA-seq", True )

    reads, length = self.parseSeqLength( ds_name + "_seqLen.txt", ds_name, 0, "R1" )
    cells, cell_reads, umi_count, exons, mito_percent, barcode_length = self.rna_10x( ds_name )

    cells_score = self.scoreHigh( cells, 500.0, 1000.0 )
    exon_score = self.scoreHigh( umi_count, 40.0, 60.0 )
    frac_cells_score = self.scoreHigh( cell_reads, 40.0, 70.0 )
    mito_score = self.scoreLow( mito_percent, 50.0, 20.0 )
    trim_score = self.scoreHigh( length, 50.0, 9999999.9 )
    umi_score = self.scoreHigh( umi_count, 500.0, 1500.0 )

    min_score = min( cells_score, exon_score, frac_cells_score, mito_score, trim_score, umi_score )
    ds_status = "not passed"
    ave_score = 0.0
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = numpy.mean( [cells_score, exon_score, frac_cells_score, mito_score, trim_score, umi_score] )


    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Mean trimmed sequence length score", trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Estimated number of cells score", cells_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Reads Mapped to Exons score", exon_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Fraction Reads in Cells score", frac_cells_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Median UMI counts per cell score", umi_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Cells with >= 20% mitochondrial reads % score", mito_score, False )

    self.global_metrics.setMetric( "Sequencing depth", reads, False )
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )
    self.global_metrics.setMetric( "Estimated number of cells", cells, False )
    self.global_metrics.setMetric( "Median UMI counts per cell", umi_count, False )
    self.global_metrics.setMetric( "Reads Mapped to Exons", exons, False )
    self.global_metrics.setMetric( "Fraction Reads in Cells", cell_reads, False )
    self.global_metrics.setMetric( "Cells with >= 20% mitochondrial reads %", mito_percent, False )
    self.toJson( ds_name )


  #############################################################################
  #
  # sn5mC-seq QA/QC parser.
  #
  def snmc_seq( self, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "snmC-seq", True )

    trim_length = 0
    reads, length = self.parseSeqLength( ds_name + "_seqLen.txt", ds_name, trim_length, "R1" )

    aligned, mapping_efficiency, cg, ch, ccc, coverage = self.parseBismark( "snmC_log.txt" )

    if path.exists( "snmC_log2.txt" ):
      aligned_r2, mapping_efficiency_r2, cg_r2, ch_r2, ccc_r2, coverage_r2 = self.parseBismark( "snmC_log2.txt" )

      if ("R1" in self.ds_metrics.keys()) == False:
        self.ds_metrics[ "R1" ] = Metrics()
      self.ds_metrics[ "R1" ].setMetric( "R1 aligned reads", aligned, False )
      self.ds_metrics[ "R1" ].setMetric( "R1 mapping efficiency", mapping_efficiency, False )
      self.ds_metrics[ "R1" ].setMetric( "R1 cell level mCG", cg, False )
      self.ds_metrics[ "R1" ].setMetric( "R1 cell level mCH", ch, False )
      self.ds_metrics[ "R1" ].setMetric( "R1 genome coverage", coverage, False )

      if ("R2" in self.ds_metrics.keys()) == False:
        self.ds_metrics[ "R2" ] = Metrics()
      self.ds_metrics[ "R2" ].setMetric( "R2 aligned reads", aligned_r2, False )
      self.ds_metrics[ "R2" ].setMetric( "R2 mapping efficiency", mapping_efficiency_r2, False )
      self.ds_metrics[ "R2" ].setMetric( "R2 cell level mCG", cg_r2, False )
      self.ds_metrics[ "R2" ].setMetric( "R2 cell level mCH", ch_r2, False )
      self.ds_metrics[ "R2" ].setMetric( "R2 genome coverage", coverage_r2, False )

      aligned = aligned + aligned_r2
      mapping_efficiency = max(mapping_efficiency, mapping_efficiency_r2)
      cg = max(cg, cg_r2)
      ch = max(ch, ch_r2)
      coverage = (coverage + coverage_r2)

    ccc_score = self.scoreLow( ccc, 0.03, 0.02 )
    cg_score = self.scoreHigh( cg, 50.0, 70.0 )
    ch_score = self.scoreLow( ch, 0.08, -1.0 )
    coverage_score = self.scoreHigh( coverage, 1.0, 5.0 )
    trim_score = self.scoreHigh( length, 30.0, 9999999.9 )
    mapped_score = self.scoreHigh( mapping_efficiency, 40.0, 70.0 )
    reads_score = self.scoreHigh( aligned, 10000.0, 1000000.0 )

    min_score = min( ccc_score, cg_score, ch_score, coverage_score, mapped_score, reads_score, trim_score )
    ds_status = "not passed"
    ave_score = 0.0
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = numpy.mean( [ccc_score, cg_score, ch_score, coverage_score, mapped_score, reads_score, trim_score] )

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Mean trimmed sequence length score", trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "CCC methylation rate score", ccc_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Cell level mCG score", cg_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Cell level mCH score", ch_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Percent genome score", coverage_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Reads per cell score", reads_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Read level %uniquely_mapped score", trim_score, False )

    self.global_metrics.setMetric( "Mean trimmed sequence length", length, False )
    self.global_metrics.setMetric( "CCC methylation rate", ccc, False )
    self.global_metrics.setMetric( "Cell level mCG", cg, False )
    self.global_metrics.setMetric( "Cell level mCH", ch, False )
    self.global_metrics.setMetric( "Percent genome", coverage, False )
    self.global_metrics.setMetric( "Reads per cell", aligned, False )
    self.global_metrics.setMetric( "Read level %uniquely mapped", mapping_efficiency, False )
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )
    self.toJson( ds_name )


  #############################################################################
  #
  # WGBS QA/QC parser.
  #
  def wgbs( self, ds_name ):
    self.global_metrics.setMetric( "Dataset name", ds_name, True )
    self.global_metrics.setMetric( "Dataset assay", "WGBS", True )

    reads, length = self.parseSeqLength( ds_name + "_seqLen.txt", ds_name, 0, "R1" )
    dup_rate = self.parseDuplicates( "deduplicated.txt", ds_name )
    nrf, ct_unique, aligned_ct, ct_percent = self.parseQc( ds_name + "_ct_bt2.log" )
    nrf, ga_unique, aligned_ga, ga_percent = self.parseQc( ds_name + "_ga_bt2.log" )
    aligned = aligned_ct + aligned_ga
    unique = (aligned * 100.0) / reads
    conversion_rate = self.bisulfite()
    peak_name = ds_name + "_macs2_B_peaks.broadPeak"
    cg_coverage, cg_bins = self.count_cg( peak_name, ds_name )
    chr1_bins = app.count_chr1_cg()

    trim_score = self.scoreHigh( length, 50.0, 9999999.9 )
    unique_score = self.scoreHigh( unique, 40.0, 70.0 )
    depth_score = self.scoreHigh( reads, 30000000, 999900000 )
    bisulfite_score = self.scoreHigh( conversion_rate, 98.0, 99.0 )
    min_score = min( bisulfite_score, depth_score, trim_score )
    dup_score = self.scoreLow( dup_rate, 0.15, -1.0 )
    ds_status = "not passed"
    ave_score = 0.0
    if ( min_score > 0.0 ):
      ds_status = "pass pending replication"
      ave_score = numpy.mean( [bisulfite_score, depth_score, dup_score, trim_score, unique_score ] )

    self.ds_metrics[ "QC" ] = Metrics()
    self.ds_metrics[ "QC" ].setMetric( "Sequencing depth score", depth_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Mean trimmed sequence length score", trim_score, False )
    self.ds_metrics[ "QC" ].setMetric( "C to T conversion rate (Bisulfite conversion rate) score", bisulfite_score, False )
    self.ds_metrics[ "QC" ].setMetric( "PCR duplicates rate score", dup_score, False )
    self.ds_metrics[ "QC" ].setMetric( "Uniquely mapped% (uniquely mapping efficiency) score", unique_score, False )

    self.global_metrics.setMetric( "PCR duplicates rate", dup_rate, False )
    self.global_metrics.setMetric( "C to T conversion rate (Bisulfite conversion rate)", conversion_rate, False )
    self.global_metrics.setMetric( "Uniquely mapped% (uniquely mapping efficiency)", unique, False )
    self.global_metrics.setMetric( "Sequencing depth", reads, False )
    self.global_metrics.setMetric( "QC status", ds_status, True )
    self.global_metrics.setMetric( "QC score", ave_score, False )
    self.toJson( ds_name )

###############################################################################
#
# This QcParser main program.
#
print( "QcParser running" )
app = QcParser ()
arg_count = len(sys.argv)
if ( arg_count >= 3 ):
  ds_name = sys.argv[1]
  assay_type = sys.argv[2]
  len_name = sys.argv[3]
  if assay_type == "atac-seq" or assay_type == "ATAC-seq" or assay_type == "ATACseq":
    print( "Processing ATAC-seq assay" )
    app.atac_seq( ds_name )

  elif assay_type == "ChIPmentation" or assay_type == "chipmentation": 
    antibody = sys.argv[3]
    app.chipmentation( ds_name, antibody )

  elif assay_type == "EPIC" or assay_type == "epic":
    dir_name = sys.argv[3]
    app.epic( ds_name, dir_name )

  elif assay_type == "massCyTOF" or assay_type == "mass-CyTOF" or assay_type == "massCytof":
    print( "Processing Mass Cytometry assay" )
    app.mass_cytof( ds_name )

  elif assay_type == "mint-ChIP" or assay_type == "mint-chip" or assay_type == "mintchip":
    print( "Processing Mint-ChIP assay" )
    # factor = int(sys.argv[3])
    app.mintchip( ds_name )

  elif assay_type == "medip" or assay_type == "meDIP" or assay_type == "meDIP-seq":
    print( "Processing MeDIP-seq assay" )
    app.medip( ds_name )

  elif assay_type == "mirna" or assay_type == "miRNA" or assay_type == "miRNA-seq":
    print( "Processing miRNA assay" )
    app.mirna( ds_name )

  elif assay_type == "multiome" or assay_type == "Multiome":
    print( "Processing multiome assay" )
    app.multiome( ds_name )

  elif assay_type == "rna-seq" or assay_type == "RNA-seq" or assay_type == "rnaseq":
    print( "Processing RNA-seq assay" )
    app.rna_seq( ds_name )

  elif assay_type == "scATAC-seq" or assay_type == "sc_atac_seq" or assay_type == "scATACseq":
    print( "Processing single cell ATAC-seq assay" )
    app.sc_atac_seq( ds_name, len_name )

  elif assay_type == "scRNA-seq" or assay_type == "sc_rna_seq" or assay_type == "scRNAseq":
    print( "Processing single cell RNA-seq assay" )
    app.sc_rna_seq( ds_name )

  elif assay_type == "snmC-seq" or assay_type == "sn5mC-seq" or assay_type == "snmCseq":
    print( "Processing single nuclear 5-methyl C-seq assay" )
    app.snmc_seq( ds_name )

  elif assay_type == "wgbs" or assay_type == "wgbs-seq" or assay_type == "WGBS" or assay_type == "WGBS-seq":
    print( "Processing WGBS-seq assay" )
    app.wgbs( ds_name )

  elif assay_type == "testing":
    # e = app.bisulfite()
    # print( "Bisulfite conversion efficiency", e)
    app.count_chr1_cg()

  else:
    print( "Unknown assay type: " + assay_type )
else:
  print( "To run: python QcParser.py <dataset name> <assay name>" )
