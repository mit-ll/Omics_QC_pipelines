

/******************************************************************************/
/**
  @author	    Darrell O. Ricke, Ph.D.

  **This class models a sequence feature**
  
  **Author:**  Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
   Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
   License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
  
  **RAMS request ID 1021178**
  
  **Overview:**
  This class models a sequence feature.
  
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
public class Feature {

	
/******************************************************************************/

  private String  feature_name = "";		// /name=value
  
  private String  feature_value = "";


/******************************************************************************/
  public String getFeatureName () { return feature_name;  } 


/******************************************************************************/
  public String getFeatureValue () { return feature_value;  } 


/******************************************************************************/
  public void setFeatureName ( String value ) { feature_name = value;  } 


/******************************************************************************/
  public void setFeatureValue ( String value ) { feature_value = value;  } 


/******************************************************************************/
  private void parseFeature ( String value )
  {
    // Assert: value
    if ( value.length () <= 3 )  return;
    int index = value.indexOf ( "=" );
    if ( index < 0 )  return;

    feature_name = value.substring( 1, index - 2 );
    
    // Remove the quote after the equal sign.
    if ( value.length () > index + 1 )
      if ( value.charAt ( index + 1 ) == '"' )  index++;
    
    if ( value.length () > index + 1 )
      feature_value = value.substring ( index+1 );

    // Remove trailing double quote.
    if ( feature_value.charAt ( feature_value.length () - 1 ) == '"' )
      feature_value = feature_value.substring ( 0, feature_value.length () - 1 );
  }  // method getValue	
  
  
/******************************************************************************/
}  /* Feature */
