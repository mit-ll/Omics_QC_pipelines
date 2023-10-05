

/******************************************************************************/
/**
  @author      Darrell O. Ricke, Ph.D.

  **This class models an alignment hit**
  
  **Author:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
   Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
   License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
  
  **RAMS request ID 1021178**
  
  **Overview:**
  This class models an alignment hit.
  
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

public class Hit {

  
/******************************************************************************/

  // Query sequence details:
  private String  query_name = "";
  
  private String  query_description = "";
  
  private int query_length = 0;
  
  private StringBuffer query_sequence = new StringBuffer( 1024 );
  
  private char query_strand = ' ';
  
  private int query_start = 0;
  
  private int query_end = 0;
  
  private String query_gene = "";
  
  private String query_organism = "";
  
  private int query_organism_hits = 0;
  
  private String query_product = "";
  
  private String query_taxonomy = "";
 
  
  // Target sequence details:
  private String  target_name = "";
  
  private String  target_description = "";
  
  private int  target_length = 0;
  
  private StringBuffer target_sequence = new StringBuffer( 1024 );

  private char target_strand = ' ';
  
  private int target_start = 0;
  
  private int target_end = 0;
  
  private String target_gene = "";
  
  private String target_organism = "";
  
  private String target_product = "";
  
  private String target_taxonomy = "";

  
  // Alignment attributes:
  private int gap_characters = 0;

  private int identities = 0;
  
  private int number_gaps = 0;
  
  private int percent = 0;    // percent identity
  
  
/******************************************************************************/
  public Hit()
  {
  initialize();
  }  // constructor Hit
  
  
/******************************************************************************/
  public void initialize()
  {
    // Query sequence details:
    query_name = "";
    query_description = "";
    query_length = 0;
    query_sequence = new StringBuffer( 1024 );
    query_strand = ' ';
    query_start = 0;
    query_end = 0;
    query_gene = "";
    query_organism = "";
    query_organism_hits = 0;
    query_product = "";
    query_taxonomy = "";
    
    // Target sequence details:
    target_name = "";
    target_description = "";
    target_length = 0;
    target_sequence = new StringBuffer( 1024 );
    target_strand = ' ';
    target_start = 0;
    target_end = 0;
    target_gene = "";
    target_organism = "";
    target_product = "";
    target_taxonomy = "";
    
    // Alignment attributes:
    identities = 0;  
    number_gaps = 0;
    gap_characters = 0; 
    percent = 0;
  }  // method initialize
  
  
/******************************************************************************/
  // This method adds one alignment segment to the alignment hit.
  public void addAlignment
      ( int qry_start, 
        int qry_end, 
        String qry_seq, 
        int subject_start, 
        int subject_end, 
        String subject_seq )
  {
    query_sequence.append( qry_seq );
    target_sequence.append( subject_seq );
  
    if ( qry_start < qry_end )
    {
      if ( query_start == 0 )  query_start = qry_start;
      if ( qry_end > query_end )  query_end = qry_end;
    }
    else
    {
      if ( query_start == 0 )  query_start = qry_end;
      if ( ( qry_end < query_start ) && ( qry_end > 0 ) )  query_start = qry_end;
      if ( qry_start > query_end )  query_end = qry_start;
    }  // if
    
    if ( subject_start < subject_end )
    {
      if ( target_start == 0 )  target_start = subject_start;
      if ( subject_end > target_end )  target_end = subject_end;
    }
    else
    {
      if ( target_start == 0 )  target_start = subject_end;
      if ( ( subject_end < target_start ) && ( subject_end > 0 ) )  target_start = subject_end;
      if ( subject_start > target_end )  target_end = subject_start;
    }  // if
  }  // addAlignment
  
  
/******************************************************************************/
  public void addQuerySequence( String sequence ) { query_sequence.append( sequence ); }
  
  
/******************************************************************************/
  public void addTargetSequence( String sequence ) { target_sequence.append( sequence ); }
  
  
/******************************************************************************/
  public String getQueryName() { return query_name;  } 
  
  
/******************************************************************************/
  public String getQueryDescription() { return query_description;  } 
  
  
/******************************************************************************/
  public int getQueryLength() { return query_length;  } 
  
  
/******************************************************************************/
  public String getQuerySequence() { return query_sequence.toString();  } 
  
  
/******************************************************************************/
  public int getQuerySize() { return query_end - query_start + 1;  } 
  
  
/******************************************************************************/
  public char getQueryStrand() { return query_strand;  } 
  
  
/******************************************************************************/
  public int getQueryStart() { return query_start;  } 
  
  
/******************************************************************************/
  public int getQueryEnd() { return query_end;  } 
  
  
/******************************************************************************/
public String getQueryGene() { return query_gene;  } 
    
    
/******************************************************************************/
public String getQueryOrganism() { return query_organism;  } 


/******************************************************************************/
public String getQueryProduct() { return query_product;  } 


/******************************************************************************/
public String getQueryTaxonomy() { return query_taxonomy;  } 



  
/******************************************************************************/
  public String getTargetName() { return target_name;  } 
  
  
/******************************************************************************/
  public String getTargetDescription() { return target_description;  } 
  
  
/******************************************************************************/
  public int getTargetLength() { return target_length;  } 
  
  
/******************************************************************************/
  public String getTargetProduct() { return target_product;  } 
  
  
/******************************************************************************/
  public String getTargetSequence() { return target_sequence.toString();  } 
  
  
/******************************************************************************/
  public int getTargetSize() { return target_end - target_start + 1;  } 
  
  
/******************************************************************************/
  public char getTargetStrand() { return target_strand;  } 
  
  
/******************************************************************************/
  public int getTargetStart() { return target_start;  } 
  
  
/******************************************************************************/
  public int getTargetEnd() { return target_end;  } 
    
    
/******************************************************************************/
  public String getTargetGene() { return target_gene;  } 
      
      
/******************************************************************************/
  public String getTargetOrganism() { return target_organism;  } 
  
  
/******************************************************************************/
  public String getTargetTaxonomy() { return target_taxonomy;  } 
  
  
  
  
/******************************************************************************/
  public int getGapCharacters() { return gap_characters;  } 
  
  
/******************************************************************************/
  public int getIdentities() { return identities;  } 
  
  
/******************************************************************************/
  public int getNumberGaps() { return number_gaps;  } 
  
  
/******************************************************************************/
  public int getPercent() { return percent;  } 


  
  
/******************************************************************************/
  public void setIdentities( int value ) { identities = value; }


/******************************************************************************/
  public void setGapCharacters( int value ) { gap_characters = value; }


/******************************************************************************/
  public void setPercent( int value ) { percent = value; }
  

/******************************************************************************/
  public void setQueryName( String value ) { query_name = value;  } 


/******************************************************************************/
  public void setQueryDescription( String value ) { query_description = value;  } 


/******************************************************************************/
  public void setQuerySequence( String value ) { query_sequence = new StringBuffer( value );  } 


/******************************************************************************/
  public void setQueryLength( int value ) { query_length = value;  } 


/******************************************************************************/
  public void setQueryStrand( char value ) { query_strand = value;  } 


/******************************************************************************/
  public void setQueryStrand( String value ) { query_strand = parseStrand( value ); } 


/******************************************************************************/
  public void setQueryStart( int value ) { query_start = value;  } 


/******************************************************************************/
  public void setQueryEnd( int value ) { query_end = value;  } 


/******************************************************************************/
  public void setQueryGene( String value ) { query_gene = value;  } 


/******************************************************************************/
  public void setQueryOrganism( String value ) { query_organism = value;  } 


/******************************************************************************/
  public void setQueryTaxonomy( String value ) { query_taxonomy = value;  } 


/******************************************************************************/
  public void setQueryProduct( String value ) { query_product = value;  } 


  
  
/******************************************************************************/
  public void setTargetName( String value ) { target_name = value;  } 


/******************************************************************************/
  public void setTargetDescription( String value ) { target_description = value;  } 


/******************************************************************************/
  public void setTargetSequence( String value ) { target_sequence = new StringBuffer( value ); } 


/******************************************************************************/
  public void setTargetLength( int value ) { target_length = value;  } 


/******************************************************************************/
  public void setTargetStrand( char value ) { target_strand = value;  } 


/******************************************************************************/
  public void setTargetStrand( String value ) { target_strand = parseStrand( value );  } 


/******************************************************************************/
  public void setTargetStart( int value ) { target_start = value;  } 


/******************************************************************************/
  public void setTargetEnd( int value ) { target_end = value;  } 


/******************************************************************************/
  public void setTargetGene( String value ) { target_gene = value;  } 


/******************************************************************************/
  public void setTargetOrganism( String value ) { target_organism = value;  } 


/******************************************************************************/
  public void setTargetTaxonomy( String value ) { target_taxonomy = value;  } 


/******************************************************************************/
  public void setTargetProduct( String value ) { target_product = value;  } 

  
/******************************************************************************/
  private String getValue( String pair )
  {
  int index = pair.indexOf( '=' );
  if ( index > 0 )
  {
    String value = pair.substring( index+1, pair.length() - 1);
    return value.replaceAll( "\"", "" );
  }  // if
  return "";
  }  // method getValue
  
  
/******************************************************************************/
  public void parseAnnotation()
  {
  if ( target_description.length() < 0 )  return;
  
  // System.out.println( "parseAnnotation: " + target_description );
  String tokens[] = target_description.split( " /" );
  for (int i = 0; i < tokens.length; i++ )
  {
    if ( tokens[ i ].startsWith( "/organism=" ) )  target_organism = getValue( tokens[ i ] );
    if ( tokens[ i ].startsWith( "gene=") )        target_gene     = getValue( tokens[ i ] );
    if ( tokens[ i ].startsWith( "product=" ) )    target_product  = getValue( tokens[ i ] );
    if ( tokens[ i ].startsWith( "taxonomy="))     target_taxonomy = getValue( tokens[ i ] );
    // System.out.println( " -- " + tokens[ i ] );
  }  // for
  
  /*
    tokens = @sequence_description.split(" /")
    tokens.each do |pair|
      if ( pair.length > 0 )
        tuple = pair.strip.split('=')
        @annotation[tuple[0].delete('/')] = tuple[1].delete('"') if tuple[1] != nil
      end  # if
    end  # do
  */  
  }  // method parseAnnotation
  
  
/******************************************************************************/
  // This method coverts "Plus" into '+' and "Minus" into '-'.
  private char parseStrand( String value )
  {
    if ( value.startsWith( "Plus" ) )  return '+';
    if ( value.startsWith( "Minus") )  return '-';
    return ' ';
  }  // method parseStrand
  
  
/******************************************************************************/
public String toString() { return toString( "," ); }
  
  
/******************************************************************************/
  public String toString( String delimiter )
  {
  return "\"" + query_name + "\"" + delimiter +
    query_strand + delimiter +
    query_start + delimiter +
    query_end + delimiter +
    "\"" + query_product + "\"" + delimiter +
    "\"" + query_organism + "\"" + delimiter +
    "\"" + query_taxonomy + "\"" + delimiter +
    identities + delimiter +
    percent + delimiter +
    number_gaps + delimiter +
    gap_characters + delimiter +
    "\"" + target_name + "\"" + delimiter +
    target_strand + delimiter +
    target_start + delimiter +
    target_end + delimiter +
    "\"" + target_organism + "\"" + delimiter +
    "\"" + target_product + "\"" + delimiter +
    "\"" + target_taxonomy + "\"" +
    "\"" + target_taxonomy + "\"" + delimiter +
    "\"" + query_sequence + "\"" + delimiter +
    "\"" + target_sequence + "\"";
  }  // method toString
  
  
/******************************************************************************/

}  // class Hit
