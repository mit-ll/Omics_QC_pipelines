# qc_base

**Epigenetics and Transcriptomics QA/QC Singularity pipelines base**

**Authors:** Derek Ng and Darrell O. Ricke, Ph.D.  (Darrell.Ricke@ll.mit.edu )
  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  

**RAMS request ID 1021178**

**Overview:**
Epigenetics and Transcriptomics QA/QC Singularity pipelines base.

**Citation:** None

**Disclaimer:**
DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.

This material is based upon work supported by the Defense Advanced Research 
Projects Agency under Air Force Contract No. (FA8702- 15-D-0001). Any opinions, 
findings and conclusions or recommendations expressed in this material are 
those of the author(s) and do not necessarily reflect the views of the 
Defense Advanced Research Projects Agency.

© 2023 Massachusetts Institute of Technology.

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


## Overview

This recipe is the base recipe for all the ECHO pipelines to be built on top of. This container should be built first, as it contains the basic software for each pipeline.

Once built, this container will have:

 * Java
 * Miniconda
 * R
 * Other supporting software, such as GCC, Git, wget, make, etc.
 * Several R packages

This recipe does not need to be built before building every single pipeline; it only has to be built once.

To build: 

    singularity build qc_base.sif qc_base.def
