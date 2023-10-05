
import java.util.*;

/******************************************************************************/
/**
  @author	    Darrell O. Ricke, Ph.D.

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

public class Hits {

  
/******************************************************************************/

  private Hit[] hits = null;
  
  private int total = 0;
  
  
/******************************************************************************/
  public Hits() { initialize(); }
  
  
/******************************************************************************/
  public Hits( int size ) { initialize(); setSize( size ); }
  
  
/******************************************************************************/
  public void initialize()
  {
  hits = null;
  total = 0;
  }  // method initialize
  
  
/******************************************************************************/
  // This method adds a Hit to the array.
  public void addHit( Hit hit )
  {
  if ( hits == null )  setSize( 250 );
  if ( total < hits.length )  hits[ total++ ] = hit;
  }  // method addHit
  
  
/******************************************************************************/
  public Hit[] getHits() { return hits; }
  
  
/******************************************************************************/
  public int getTotal() { return total; }
  
  
/******************************************************************************/
  // This method resets the hits.
  public void reset()
  {
  total = 0;
  if ( hits == null )  return;
  
    for ( int i = 0; i < hits.length; i++ )
      hits[ i ] = null;
  }  // method reset
  
  
/******************************************************************************/
  public void setSize( int size )
  {
  hits = new Hit[ size ];
  for (int i = 0; i < size; i++ )
    hits[ i ] = null;
  }  // method setSize
  
  
/******************************************************************************/
  private String determineQuerySpecies( HashMap<String, Integer> species )
  {  
  String genus = "";
  String organism = "";
  String tokens[] = null;
  
  Set<String> keys = species.keySet();
  for (String key : keys)
  {
      if ( species.get( key ) >= 10 )
      {
      // System.out.println( "Species: " + key + " - " + species.get( key ) );
        tokens = key.split( " " );
        if ( genus.length() < 1 )  genus = tokens[ 0 ];
        if ( organism.length() < 1 )  organism = key.trim();
        
        if ( genus.equals( tokens[ 0 ] ) == false )  return "Mixed";
      
        // Use sp. ending as a wild card to match other species.
        if ( organism.endsWith( "sp." ) )  organism = key.trim();
      }  // if
  }  // for
  
  if ( organism.length() < 1 )  return "Mixed";
  return organism;
  }  // method determineQuerySpecies
  
  
/******************************************************************************/
  private String determineQueryProduct( HashMap<String, Integer> products )
  {
  int count = 0;
  String product = "";
  int value = 0;

  Set<String> keys = products.keySet();
  for (String key : keys)
  {
    value = products.get( key );
    if ( value >= count )
    {
    product = key;
    count = value;
    }
  }  // for
  
  return product;
  }  // method determineQueryProduct
  
  
/******************************************************************************/
  private String determineQueryTaxonomy( HashMap<String, Integer> taxonomy )
  {
  String taxon[] = null;
  String tokens[] = null;
  boolean clipped = false;
  
  // Start with the longest taxonomy.
    Set<String> keys = taxonomy.keySet();
    for (String key : keys)
      if ( ( key.length() > 0 ) && ( taxonomy.get( key ) >= 10 ) )
      {
        tokens = key.split( ";" );
        if ( ( taxon == null ) || ( tokens.length > taxon.length ) )  taxon = tokens;
      }  // for
    if ( taxon == null )  return "";
    
    for (String key : keys)
      if ( ( key.length() > 0 ) && ( taxonomy.get( key ) >= 5 ) )
      {
          // Compare the taxonomy strings.
        tokens = key.split( ";" );
        // System.out.println( "--n tax: " + key + " - " + taxonomy.get( key ) );
      
        for ( int i = 0; i < taxon.length; i++ )
        {
          // Trim the taxonomy tree to matching names.
          if ( ( i < tokens.length ) &&
             ( taxon[ i ] != null ) &&
             ( taxon[ i ].trim().equals( tokens[ i ].trim() ) == false ) )
          {
            for ( int j = i; j < taxon.length; j++ )
            taxon[ j ] = null;
          clipped = true;
          }  // if
        }  // for
      }  // if
    
    // Check for no common taxonomy names.
    if ( taxon == null )  return "";
    
    StringBuffer tax = new StringBuffer();
    for ( int i = 0; i < taxon.length; i++ )
    {
      if ( taxon[ i ] != null )
      {
        if ( tax.length() > 0 )  tax.append( ";" );
        tax.append( taxon[ i ] );
      }  // if
    }  // for
    if ( clipped == false )  tax.append( "." );
    
    return tax.toString();
  }  // method determineQueryTaxonomy
  
  
/******************************************************************************/
  public void summarize()
  {
  HashMap<String, Integer> products = new HashMap<String, Integer>( 64 );
  HashMap<String, Integer> species = new HashMap<String, Integer>( 64 );
  HashMap<String, Integer> taxonomy = new HashMap<String, Integer>( 128 );
  String product_name = "";
  String organism = "";
  String taxon = "";
  // String tokens[] = null;
 
  if ( hits.length <= 0 )  return;
 
  Integer one = 1;
  int target_size = hits[ 0 ].getTargetSize();
  int identities = hits[ 0 ].getIdentities();
  int index1 = 0;
  int index2 = 0;
  
  for ( int i = 0; i < hits.length; i++ )
  {
    // System.out.println( hits[i].toString() );

    if ( (hits[i] != null) && ( hits[ i ].getTargetSize() >= (target_size * 9) / 10 ) && ( hits[ i ].getIdentities() >= (identities * 9) / 10 ) )
      {
      hits[ i ].parseAnnotation();
    
      // Tally the observed gene names.
      product_name = hits[ i ].getTargetProduct();  
      if ( products.containsKey( product_name ) )
        products.put( product_name, products.get( product_name ) + one );
      else
      products.put( product_name, one );
    
      // Tally the observed target organism names.
      organism = hits[ i ].getTargetOrganism();
      index1 = organism.indexOf( ' ' );
      index2 = 0;
      if ( index1 > 0 )  index2 = organism.indexOf( ' ', index1+1 );
      if ( index2 > 0 )  organism = organism.substring( 0, index2 );
      if ( species.containsKey( organism ) )
      species.put( organism, species.get( organism ) + one );
      else
      species.put( organism, one );
    
      // Tally the taxonomies observed.
      taxon = hits[ i ].getTargetTaxonomy();
      index2 = taxon.indexOf( '.' );
      if ( index2 > 0 )  taxon = taxon.substring( 0, index2 );
      if ( taxonomy.containsKey( taxon ) )
      taxonomy.put( taxon, taxonomy.get( taxon ) + one );
      else
      taxonomy.put( taxon, one );
    }  // if
  }  // for
 
  hits[ 0 ].setQueryProduct( determineQueryProduct( products ) );
  hits[ 0 ].setQueryOrganism( determineQuerySpecies( species ) );
  taxon = determineQueryTaxonomy( taxonomy );
  hits[ 0 ].setQueryTaxonomy( taxon );
  if ( taxon.endsWith( "." ) == false )  hits[ 0 ].setQueryOrganism( "Mixed" );
  
  if ( total > 0 )
    System.out.println( hits[ 0 ].toString() );  
  }  // method summarize
  
  
/******************************************************************************/
  public void summarize2()
  {
    if ( hits.length <= 0 )  return;
 
    for ( int i = 0; i < hits.length; i++ )
    {
      if ( (hits[i] != null) && ( hits[i].getIdentities() >= 40 ) )
        System.out.println( hits[i].toString() );
    }  // for
  }  // method summarize2
  
  
/******************************************************************************/

}  // class Hits
