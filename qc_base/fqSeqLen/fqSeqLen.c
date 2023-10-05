/**
 * **FASTQ sequence length program**
 * 
 * **Authors:**  Derek Ng, Darrell O. Ricke, Ph.D.  (mailto: Darrell.Ricke@ll.mit.edu)
 *  Copyright:  Copyright (c) 2023 Massachusetts Institute of Technology 
 *  License:    GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
 * 
 * **RAMS request ID 1021178**
 * 
 * **Overview:**
 * FASTQ sequence length program
 * 
 * **Citation:** None
 * 
 * **Disclaimer:**
 * DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.
 * 
 * This material is based upon work supported by the Defense Advanced Research 
 * Projects Agency under Air Force Contract No. (FA8702- 15-D-0001). Any opinions, 
 * findings and conclusions or recommendations expressed in this material are 
 * those of the author(s) and do not necessarily reflect the views of the 
 * Defense Advanced Research Projects Agency.
 * 
 * © 2023 Massachusetts Institute of Technology
 * 
 * The software/firmware is provided to you on an As-Is basis
 * 
 * Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS
 * Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright notice,
 * U.S. Government rights in this work are defined by DFARS 252.227-7013 or
 * DFARS 252.227-7014 as detailed above. Use of this work other than as specifically
 * authorized by the U.S. Government may violate any copyrights that exist in this work.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>

#define PLUS_SIZE 4

struct Read {
        char* name; //Name
        char* seq; //Sequence
        char* qual; //Quality
};

/**
 * Print the usage information
 * @param status The exit status to be used
 */
void usage(int status) {
        printf("\nProgram Name: fqSeqLen\n\n");
        printf("Usage: path-to-program/fqLineLength [options] <fastq1> [fastq2]\n\n");
        printf("Options: -n <int>       Length of the name of a fastq read (default = 128)\n");
        printf("         -s <int>       Length of the sequence of a fastq read (default = 256)\n");
	printf("         -o <string>    Name of the file to write results (write to console if not specified)\n");
	printf("         -h             Print this help message\n\n");
        exit(status);
}

/**
 * Get command line options
 * @param argc The number of arguments passed into the C program
 * @param argv The array of arguments
 * @param nameLen Pointer to the length of a name in a Read
 * @param seqLen Pointer to the length of a sequence in a Read
 * @param outputFile Pointer to the name of the output file
 * @param writeToFile Pointer to a boolean where true means the output is written to a file
 * 	and false means the output is written to the console
 * @return the number of command line arguments used to store the options
 * 	Ex. "-s 180" is two args; "-h" is one arg; "-o example.txt -s 164 -n 78" is six args
 */
int getOptions(int argc, char* argv[], int* nameLen, int* seqLen, char** outputFile, bool* writeToFile) {
	char* ptr = 0;
	long argl = 0;
	int opt = 0;
	int optCount = 0;
	
	//get options: name length and sequence length
	while((opt = getopt(argc, argv, "n:s:o:h")) != -1) {
		switch(opt) {
			case 'n': //name length
				argl = strtol(optarg, &ptr, 10);
				*nameLen = (int) argl + 1; //account for \0
				if(*nameLen - 1 < 1) {
					printf("\033[1;31mError: %d is not a valid name length. Must be an integer greater than 0.\n\n\033[0m", *nameLen - 1);
					usage(EXIT_FAILURE);
				}
				optCount += 2;
				break;
			case 's': //sequence length
				argl = strtol(optarg, &ptr, 10);
				*seqLen = (int) argl + 1; //account for \0
				if(*seqLen - 1 < 1) {
					printf("\033[1;31mError: %d is not a valid sequence length. Must be an integer greater than 0.\n\n\033[0m", *seqLen - 1);
					usage(EXIT_FAILURE);
				}
				optCount += 2;
				break;
			case 'o': //filename to write to
				*outputFile = (char*) malloc(strlen(optarg) + 1);
				strcpy(*outputFile, optarg);
				*writeToFile = true;
				optCount += 2;
				break;
			case 'h':
				usage(EXIT_SUCCESS);
				optCount++;
			case '?':
				printf("\033[1;31mError: Unknown option: %c\n\033[0m", optopt);
				usage(EXIT_FAILURE);
		}
	}
	
	//defaults if options not provided
	if(*nameLen == 0) {
		*nameLen = 128;
	}
	if(*seqLen == 0) {
		*seqLen = 256;
	}

	return optCount;
}

/**
 * Get file extension
 * @param filename The name of the fastq file
 */
char* getFileExt(char* filename) {
    char* dot = strrchr(filename, '.');
    if(!dot) {
           return "";
    }
    return dot + 1;
}

/**
 * Determines the length of a fastq sequence excluding any \'N\'s at the ends of the sequence
 * @param seq A fastq sequence
 * @return the length of the fastq sequence
 */
int seq_length(char* seq) {
        int startIdx = 0;
        int endIdx = strlen(seq);

        //Determines the position of the first character that is not an \'N\'
        for(int i = 0; i < endIdx; i++) {
                if(seq[i] != 'N') { //Character not an \'N\'
                        startIdx = i;
                        break;
                } else if(i + 1 == endIdx) { //Sequence contains all \'N\'s
                        startIdx = i + 1;
                }
        }

        // Determines the position of the last character that is not an \'N\'
        if(startIdx != endIdx) {
                for(int j = endIdx - 1; j > 0; j--) {
                        if(seq[j] != 'N') { //Character not an \'N\'
                                endIdx = j + 1;
                                break;
                        }
                }
        }
        return (double)(endIdx - startIdx);
}

/**
 * Count the number of Reads in a file and calculates the average trimmed length of the Reads
 * @param filename The name of the file containing the Reads
 * @param nameLen The length of a name in a Read
 * @param seqLen The length of a sequence in a Read
 * @param numReads A pointer to the number of Reads in a file
 * @param avg A pointer to the average trimmed length of the Reads
 */
void countReadsAndCalcAvgLen(char* filename, int nameLen, int seqLen, int* numReads, double* avg) {
        //get file
        FILE* fp = fopen(filename, "r");
        if(!fp) {
                printf("\033[1;31mError: Cannot find the file %s\n\033[0m", filename);
                usage(EXIT_FAILURE);
        } else if(!(strcmp("fastq", getFileExt(filename)) != 0 ^ strcmp("fq", getFileExt(filename)) != 0)) {
		printf("\033[1;31mError: Input file must have .fastq or .fq as its file extension.\n\033[0m");
                usage(EXIT_FAILURE);
        }

        char n[nameLen];
        char s[seqLen];
        char p[PLUS_SIZE];
        char q[seqLen];
        double sum = 0.0;
        *numReads = 0;

        while(fgets(n, nameLen, fp) != NULL) {
                //gets read
                strtok(n, "\n"); //removes \n character at end of line
                fgets(s, seqLen, fp);
                strtok(s, "\n");
                fgets(p, PLUS_SIZE, fp);
                fgets(q, seqLen, fp);
                strtok(q, "\n");

                /*//put elements in Read
                struct Read* read = (struct Read*) malloc(sizeof(struct Read));
                read->name = malloc(nameLen);
                strncpy(read->name, n, nameLen);
                read->seq = malloc(seqLen);
                strncpy(read->seq, s, seqLen);
                read->qual = malloc(seqLen);
                strncpy(read->qual, q, seqLen);*/

                sum += seq_length(s); //calculates sum of trimmed sequences
                (*numReads)++; //calculates number of reads
        }
        fclose(fp);
        *avg = sum / *numReads; //average length of the trimmed reads
}

/**
 * Prints the output, either to the console or to a file
 * @param writeToFile True if the output is to be written to a file, false if the output is to be written to the console
 * @param twoFiles True if there are two input files, false if there is one input file
 * @param outputFile The name of the file to write the output
 * @param file1 The first input file
 * @param file2 The second input file
 * @param nameLen The length of a name of a Read
 * @param seqLen The length of a sequence of a Read
 */
void printOutput(bool writeToFile, bool twoFiles, char* outputFile, char* file1, char* file2, int nameLen, int seqLen) {
	int numReads = 0;
	double avg = 0;

	countReadsAndCalcAvgLen(file1, nameLen, seqLen, &numReads, &avg); //get data for first input file
	if(writeToFile) { //writes to file
		FILE* fp = fopen(outputFile, "w");
		fprintf(fp, "Current file: %s\n\n", outputFile);
		
		fprintf(fp, "Fastq file: %s\n", file1);
		fprintf(fp, "Sequencing Depth: %d\n", numReads);
		fprintf(fp, "Mean Trimmed Sequence Length: %f\n\n", avg);

		if(twoFiles) { //two files
			countReadsAndCalcAvgLen(file2, nameLen, seqLen, &numReads, &avg);
			fprintf(fp, "Fastq file: %s\n", file2);
			fprintf(fp, "Sequencing Depth: %d\n", numReads);
			fprintf(fp, "Mean Trimmed Sequence Length: %f\n\n", avg);
		}
		fclose(fp);
	} else { //write to console
		printf("Fastq file: %s\n", file1);
		printf("Sequencing Depth: %d\n", numReads);
		printf("Mean Trimmed Sequence Length: %f\n\n", avg);
		if(twoFiles) { //two files
			countReadsAndCalcAvgLen(file2, nameLen, seqLen, &numReads, &avg);
			printf("Fastq file: %s\n", file2);
			printf("Sequencing Depth: %d\n", numReads);
			printf("Mean Trimmed Sequence Length: %f\n\n", avg);
		}
	}
}

int main(int argc, char* argv[]) {
	int nameLen = 0;
	int seqLen = 0;
	char* outputFile;
	bool writeToFile = false;
	int optCount = 0;
	bool twoFiles = false;

	if(argc < 2) { //no arguments provided
		printf("\033[1;31mError: Filename required\033[0m\n");
		usage(EXIT_FAILURE);
	} else if(argc == 2) { //one filename and no options provided
		if(strcmp("-h", argv[1]) == 0) { //argument could be -h
			usage(EXIT_SUCCESS);
		}
		nameLen = 128;
		seqLen = 256;
	} else {
		optCount = getOptions(argc, argv, &nameLen, &seqLen, &outputFile, &writeToFile);
		//if there are two input files, then difference between number of args and option args is three
		if(argc - optCount == 3) {
			twoFiles = true;
		}
	}

	if(twoFiles) { //first file stored in argv[argc - 2]; second file stored in argv[argc - 1]
		printOutput(writeToFile, twoFiles, outputFile, argv[argc - 2], argv[argc - 1], nameLen, seqLen);
	} else { //file stored in argv[argc - 1]; second argv[argc - 1] is useless - just a placeholder
		printOutput(writeToFile, twoFiles, outputFile, argv[argc - 1], argv[argc - 1], nameLen, seqLen);
	}

	return EXIT_SUCCESS;
}
