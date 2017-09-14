#!/usr/bin/env python

#reads_stats.py
#created 12 September 2017

# All progress reports are flushed to the .out file

############### SET UP #########################################################

print("HERE WE GO!", flush = True)

import argparse
print("I've got argparse!", flush = True)

import numpy as np
print("Numpy, check!", flush = True)

parser = argparse.ArgumentParser(description="Assess the overall quality of a lane of Illumina sequencing data by calculating various statistics.")
parser.add_argument('-r','--read_file', help='Add absolute path to your FASTQ file.', required=True, type=str)
parser.add_argument('-l','--read_length', help='The length of your reads in bp.', required=True, type=int)
parser.add_argument('-obp','--output_bpqual', help='Add absolute path to output file for per-base-pair mean quality across your dataset.', required=True, type=str)
parser.add_argument('-orq','--output_readqual', help='Add absolute path to output file for per-read mean quality across your dataset.', required=True, type=str)

args = parser.parse_args()

############### OPEN FASTQ FILE ################################################

# Open your FASTQ files
r1 = open(args.read_file, "r")

############### MAKE NUMPY ARRAY ###############################################

# Create empty array to store quality scores
mean_scores1 = np.zeros((args.read_length, 1), dtype = np.int64)
# Progress report
print("Created mean_scores array.", flush = True)

############### MAKE DICTIONARY FOR READ QUALITY SCORES ########################

# Make empty dictionary to hold counts of the per-read mean quality scores
read_quality_counts1 = {}
# Progress report
print("Created read quality counts dictionary.", flush = True)

############### LOOP THROUGH FASTQ FILE ########################################

# Start reading through lines in the FASTQ file
linecount = 0
for line1 in r1:
    line1 = line1.strip("\n")
    linecount += 1
    if linecount % 1000000 == 0:     # Porgress report
        print("I am on line " + str(linecount) + ".", flush = True)
    if linecount % 4 == 0:     # Pull quality lines
        quality1 = line1

############### MAKE ARRAY AND DICTIOANRY ######################################

        index = 0
        read1 = []
        for char in quality1:
            phred = ord(char) - 33     # Convert ASCII quality scores to numeric qualty scores
            mean_scores1[index] += phred     # Add score to running sum in per-base quality array
            read1.append(phred)    # Appent score to list of scores for bases in this read
            index += 1
        mean_read1 = int(np.mean(read1))     # Calculate mean score for this read and cast as an integer for "automatic" binning
        if mean_read1 in read_quality_counts1:    # Count occurences of the read means in this loop
            read_quality_counts1[mean_read1] += 1     # Increment
        else:
            read_quality_counts1[mean_read1] = 1     # Increment
# Progress report
print("Created initial array of quality scores and dictionary of mean read quality counts.", flush = True)

############### CALCULATE PER-BASE MEAN SCORE ##################################

# Calulate mean quality score per base position
means_a = mean_scores1 / (linecount/4)

############### WRITE TO OUTPUT FILES ###########@@@@@##########################

# Convert numpy array to a list of lists for more aesthetic writing to output
means = means_a.tolist()

# Write results to output: per-base-pair mean quality score
with open(args.output_bpqual, 'w') as obp:
    obp.write("Base_Pair" + "\t" + "Mean_Score"  + "\n")
    for i in range(1, args.read_length + 1):
         obp.write(str(i) + "\t" + str(means[i-1][0]) + "\n")

# Write results to output: per read mean quality score
with open(args.output_readqual, 'w') as orq:
    orq.write("Quality_Score" + "\t" + "Read_Count" + "\n")
    for k,v in read_quality_counts1.items():
        orq.write(str(k) + "\t" + str(v) + "\n")

############### CLOSING STATEMENT ##############################################
print("\n", flush = True)
print("Guess what? I finished!", flush = True)
