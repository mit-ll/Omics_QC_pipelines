
// import java.io.*;
// import java.sql.*;
// import java.util.Vector;

// import InputTools;

/******************************************************************************/
/**
  @author      Darrell O. Ricke, Ph.D.  

  **This class parses a BLAST output file**
  
  **Author:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
   Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
   License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
  
  **RAMS request ID 1021178**
  
  **Overview:**
  This class parses a BLAST output file.
  
  **Citation:** None
  
  **Disclaimer:**
  DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.
  
  This material is based upon work supported by the Defense Advanced Research 
  Projects Agency under Air Force Contract No. (FA8702- 15-D-0001). Any opinions, 
  findings and conclusions or recommendations expressed in this material are 
  those of the author(s) and do not necessarily reflect the views of the 
  Defense Advanced Research Projects Agency.
  
  © 2023 Massachusetts Institute of Technology
  
  The software/firmware is provided to you on an As-Is basis
  
  Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS
  Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright notice,
  U.S. Government rights in this work are defined by DFARS 252.227-7013 or
  DFARS 252.227-7014 as detailed above. Use of this work other than as specifically
  authorized by the U.S. Government may violate any copyrights that exist in this work.
 
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/******************************************************************************/

public class BlastParser extends Object
{

/******************************************************************************/

  private InputTools blast_file = new InputTools ();  // Blast input file

  private String line = null;         // Current line of the file

  // Blast file header fields.
  private BlastHeader blast_header = new BlastHeader ();


/******************************************************************************/
public BlastParser ()
{
  initialize ();
}  // constructor BlastParser


/******************************************************************************/
public void initialize ()
{
  line = null;
}  // method initialize


/******************************************************************************/
  // This method parses one alignment segment (Query ||| Sbjct).
  private Hit parseAlignment( Hit hit, String line )
  {
  // System.out.println ( " ** " + line );
  int query_start = InputTools.getInteger( line.substring( 6 ) );
  int seq_start = 11;
  while ( ( seq_start < line.length() ) && ( line.charAt( seq_start ) == ' ' ) )
    seq_start++;
  int index = line.indexOf( ' ', seq_start );
  // System.out.println( " -- query last space index = " + index );
        // System.out.println( seq_start + ": " + line );

  int query_end = InputTools.getInteger( line.substring( index ) );
  String query_seq = "";
  if ( index > seq_start )  query_seq = line.substring( seq_start, index );
  // System.out.println( "  -- query_seq '" + query_seq + "'" );

  line = blast_file.nextLine ().toString ();
  // String pipes = line.trim();
  // System.out.println( "  -- pipes     '" + pipes + "'" );
  
  line = blast_file.nextLine ().toString ();
  // System.out.println( "  -- " + line );
  
    int target_start = InputTools.getInteger( line.substring( 6 ) );
    // System.out.println( "  -- target_start " + target_start );
    index = line.indexOf( ' ', seq_start );
    // System.out.println( " -- target last space index = " + index );
    int target_end = InputTools.getInteger( line.substring( index ) );
    String target_seq = "";
    if ( index > seq_start )  target_seq = line.substring( seq_start, index );
    // System.out.println( "  -- target_seq'" + target_seq + "'" );
    
    // Add this alignment segment to the alignment.
    hit.addAlignment( query_start, query_end, query_seq, target_start, target_end, target_seq );
  return hit;
  }  // method parseAlignment
  

/******************************************************************************/
  private void parseHeader ()
  {
    blast_header.initialize();
    line = blast_file.nextLine().toString ();
    String tokens[] = line.split( "\\s" );
    blast_header.setProgramName( tokens[ 0 ] );
    blast_header.setProgramVersion( tokens[ 1 ] );
    // System.out.println( "Program " + blast_header.getProgramName() + "; version " + blast_header.getProgramVersion() );

    while ( ( blast_file.isEndOfFile () != true ) &&
            ( line.startsWith ( "Query=" ) != true ) )
    {
      if ( line.length () > 0 )
      {
      if ( line.startsWith( "Database:" ) )
      {
        tokens = line.split( ":" );
        blast_header.setDatabase( tokens[ 1 ].trim() );
          line = blast_file.nextLine ().toString ();
          tokens = line.trim().split( ";" );
          blast_header.setDatabaseSequences( InputTools.getInteger( tokens[ 0 ] ) );
          blast_header.setDatabaseLetters( InputTools.getInteger( tokens[ 1 ] ) );
          // System.out.println( "Database: " + blast_header.getDatabase() + "; sequences = " + blast_header.getDatabaseSequences() );
      }  // if
      }  // if

      line = blast_file.nextLine ().toString ();
    }  // while
  }  // method parseHeader


/******************************************************************************/
// This method parses a BLAST alignment hit.
  private Hit parseHit( String query_name, int query_length, String line )
  {
  Hit hit = new Hit();
  // String tokens[] = null;
  
    // System.out.println ( " -- " + line );
    hit.setQueryName( query_name );
    hit.setQueryLength( query_length );
    
    // Parse the target header.
    StringBuffer header = new StringBuffer( 1024 );
    while ( ( blast_file.isEndOfFile () != true ) &&
            ( line.startsWith ( "Length=" ) != true ) )
    {
      header.append( line );
      line = blast_file.nextLine().toString();      
    }  // while
    int index = header.indexOf( " ", 2 );
    if ( index > 0 )
    {
      hit.setTargetName( header.substring( 1, index ) );
      hit.setTargetDescription( header.substring( index+1 ) );
    }  // if
    else
      hit.setTargetName( header.substring( 1 ) );
    
    // Parse the target length.
    if ( line.startsWith( "Length=" ) )
    {
      hit.setTargetLength( InputTools.getInteger( line.substring( 7 ) ) );
      // System.out.println( "  -- target length = " + query_length );
      line = blast_file.nextLine().toString();    // blank line 
      line = blast_file.nextLine().toString();    // Score =
      line = blast_file.nextLine().toString();    // Identities =
    }  // if
    
    // Parse the Identities line.
    if ( line.startsWith( " Identities =" ) )
    {
      hit.setIdentities( InputTools.getInteger( line.substring( 13 ) ) );
      index = line.indexOf( '(' );
      if ( index > 13 )  hit.setPercent( InputTools.getInteger( line.substring( index + 1 ) ) );
      index = line.indexOf( "Gaps =" );
      if ( index > 13 )  hit.setGapCharacters( InputTools.getInteger( line.substring( index + 6 ) ) );
    }  // if
    
    line = blast_file.nextLine().toString();    // Strand=
    // Parse the alignment strands.
    if ( line.startsWith( " Strand=") )
    {
      index = line.indexOf( '/', 8 );
      if ( index > 8 )
      {
        hit.setQueryStrand( line.substring( 8, index ) );
        hit.setTargetStrand( line.substring( index + 1 ) );
      }  // if
      line = blast_file.nextLine().toString();    // Strand=      
    }  // if
    
    // Parse the alignment.
    while ( ( blast_file.isEndOfFile () != true ) &&
            ( line.startsWith( "Query=" ) != true ) &&
            ( line.startsWith( ">" ) != true ) )
    {
      // Check for a new alignment segment.
      if ( line.startsWith( "Query " ) )
      hit = parseAlignment( hit, line );

      line = blast_file.nextLine().toString();      
    }  // while
 
    // System.out.println( hit.toString() );
    
    return hit;
  }  // method parseHit
  

/******************************************************************************/
// This method parses the BLAST results for query sequences.
  private void parseQueries()
  {
  Hits hits = new Hits( 250 );  // Alignment hits
  String tokens[] = null;
  String query_name = "";      // query sequence name
  int query_length = 0;      // query sequence length
  
    while ( blast_file.isEndOfFile () != true )
    {
      if ( line.length () > 0 )
      {
      if ( line.startsWith( "Query=" ) )
      {
        // Summarize the hits for the last query sequence.
        if ( hits.getTotal() > 0 )
        {
        hits.summarize2();
        hits.reset();
        }  // if
        
        tokens = line.split( "\\s" );
        query_name = tokens [ 1 ];
          line = blast_file.nextLine().toString();
        // System.out.println( "Name: " + query_name );        
      }  // if
        
      if ( line.startsWith( "Length=" ) )
      {
        query_length = InputTools.getInteger( line.substring( 7 ) );
          line = blast_file.nextLine().toString();
        // System.out.println( "  -- query length = " + query_length );
      }  // if
      
      if ( line.startsWith( ">" ) )
        hits.addHit( parseHit( query_name, query_length, line ) );
      }  // if
      
      line = blast_file.getLine();
      if ( ( line.startsWith( ">" ) != true ) && ( line.startsWith( "Query=" ) != true ) )
        line = blast_file.nextLine().toString();
    }  // while 

    // Summarize the hits for the last query sequence.
    if ( hits.getTotal() > 0 )
  {
    hits.summarize2();
    hits.reset();
  }  // if
  }  // method parseQueries
  

/******************************************************************************/
// This method processes a BLAST output file.
  public void processFile ( String file_name )
  {
    // Setup Blast entry processor.

    // Set the input file name.
    blast_file.initialize();
    blast_file.setFileName( file_name );

    // Open the input file.
    blast_file.openFile();
    parseHeader();

    // Process the input file.
    parseQueries();

    // Close input file.
    blast_file.closeFile();
  }  // method processFile


/******************************************************************************/
// Command line usage information.
  private void usage ()
  {
    System.out.println ( "The command line syntax for this program is:" );
    System.out.println ();
    System.out.println ( "java -jar BlastParser.jar <Blast_file>" );
    System.out.println ();
    System.out.print   ( "where <Blast_file> is the file name of a " );
    System.out.println ( "Blast output file." );
  }  // method usage


/******************************************************************************/
// Main program.
  public static void main ( String [] args )
  {
    BlastParser app = new BlastParser ();

    if ( args.length != 1 )
      app.usage ();
    else
    {
      app.processFile ( args[ 0 ] );
    }  // else
  }  // method main

  
/******************************************************************************/
}  // class BlastParser

