# fqSeqLen

**This program counts the number of sequences in a fastq file and calculates the average length of the sequences, excluding any N's located at the two ends of a sequence.**

**Epigenetics and Transcriptomics QA/QC Singularity pipelines**

**Authors:** Derek Ng and Darrell O. Ricke, Ph.D.  (Darrell.Ricke@ll.mit.edu )
  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  

**RAMS request ID 1021178**

**Overview:**
This program counts the number of sequences in a fastq file and calculates the average length of the sequences.

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


## Instructions
1. To build this program, run the following commands:

    ```
    make
    ```

2. To run the program, type the following:

    ```
    fqSeqLen [options] <fastq1> [fastq2]
    ```

## Run Options

	    -n <int>       Length of the name of a fastq read (default = 80)

	    -s <int>       Length of the sequence of a fastq read (default = 160)

	    -o <string>    Name of the file to write results (write to console if not specified)

	    -h             Print this help message

