

/******************************************************************************/
/**
  **This class models the header information from a BLAST output file**
  
  **Author:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
   Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
   License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
  
  **RAMS request ID 1021178**
  
  **Overview:**
  This class models the header information from a BLAST output file.
  
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
public class BlastHeader {

  private String program_name = "";		// program name
  
  private String program_version = "";		// program version
  
  private String database = "";			// database name
  
  private int database_sequences = 0;		// number of database sequences
  
  private int database_letters = 0;		// number of database letters


/******************************************************************************/
  // Constructor BlastHeader
  public BlastHeader()
  {
    initialize();
  }  // constructor BlastHeader


/******************************************************************************/
  // Initialize class variables.
  public void initialize ()
  {
    program_name = "";
    program_version = "";
    database = "";
    database_sequences = 0;
    database_letters = 0;
  }  // method initialize
  

/******************************************************************************/
  public String getDatabase() { return database; }
  

/******************************************************************************/
  public int getDatabaseLetters() { return database_letters; }
  

/******************************************************************************/
  public int getDatabaseSequences() { return database_sequences; }

  
  /******************************************************************************/
  public String getProgramName() { return program_name; } 
  

/******************************************************************************/
  public String getProgramVersion() { return program_version; }

  
/******************************************************************************/
  public void setDatabase( String value ) { database = value; }

  
/******************************************************************************/
  public void setDatabaseLetters( int value ) { database_letters = value; }

  
/******************************************************************************/
  public void setDatabaseSequences( int value ) { database_sequences = value; }
  

/******************************************************************************/
  public void setProgramName( String value ) { program_name = value; }

  
/******************************************************************************/
  public void setProgramVersion( String value ) { program_version = value; }
  
/******************************************************************************/

}  // class BlastHeader
