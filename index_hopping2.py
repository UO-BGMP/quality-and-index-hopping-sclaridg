#!/usr/bin/env python

#index_hopping.py
#created 12 July 2017

# All progress reports are flushed to the .out file

############### SET UP #########################################################

print("I'm starting up!", flush = True)

import argparse
print("I've imported argparse.", flush = True)
print("\n", flush = True)

parser = argparse.ArgumentParser(description="Quality filter sequencing data and split files by index.")
parser.add_argument('-r2','--read2', help='Add absolute path to R2.', required=True, type=str)
parser.add_argument('-r3','--read3', help='Add absolute path to R3.', required=True, type=str)
parser.add_argument('-q','--min_qual', help='Provide integer minimum quality score cutoff.', required=True, type=int)
parser.add_argument('-i','--indices', help='Add absolute path to index TSV file.', required=True, type=str)
parser.add_argument('-o','--output_filepath', help='Add absolute path to output file for expected and unexpected index pairs.', required=True, type=str)
args = parser.parse_args()

# Set min_qual variable
min_qual = args.min_qual

############### OPEN FASTQ FILES ###############################################

# Open your sequencing files
r2 = open(args.read2, "r")
r3 = open(args.read3, "r")
print("Ooooo-eeee! I've opened your four FASTQ files and am ready for business!", flush = True)
print("\n", flush = True)

############### MAKE INDEX LIST ################################################

# Make list of indices used in sequencing (will use to make the unexpected_ and expected_index_pairs dictionaries)
index_list = []
with open(args.indices) as indices:
    for line in indices:
        line = line.strip("\n").split("\t")
        index_list.append(line[1])
print("I just added all of your indices to the index_list.", flush = True)
print("\n", flush = True)

############### MAKE EXPECTED INDEX PAIRS DICTIONARY ###########################

# Make dictionary with keys as all expected index pairs and values as their counts, initialized at 0
# Make "same_same" pairs of indexes, which are the index combinations that were used
expected_index_pairs = {}
for index in index_list:
    pair = index + "_" + index
    expected_index_pairs[pair] = 0     # Initialized at 0
    print("Added index pair " + pair + " to EXPECTED PAIRS dictionary.", flush = True)
print("\n", flush = True)
print("There are " + str(len(expected_index_pairs)) + " index pairs in the EXPECTED PAIRS dictionary.", flush = True)
print("\n", flush = True)

############### MAKE UNEXPECTED INDEX PAIRS DICTIONARY #########################

# Make list of unexpected index combinations using itertools
# From the documentation on itertools.combinations():
    # Elements are treated as unique based on their position, not on their value
    # So if the input elements are unique, there will be no repeat values in each combination.
    # This produces a list of tuples
import itertools
unexpected = list(itertools.combinations(index_list, 2))

# Make dictionary with keys as all unexpected index pairs and values as their counts
unexpected_index_pairs = {}
for pair in unexpected:
    pair = '_'.join(pair)     # This converts the tuple format to an "index1_index2" format
    unexpected_index_pairs[pair] = 0     # Initialized at 0
print("There are " + str(len(unexpected_index_pairs)) + " index pairs in the UNEXPECTED PAIRS dictionary.", flush = True)
print("\n", flush = True)

############### MAKE N-CONTAINING INDEX PAIRS DICTIONARY #######################

# Make dictionary with keys as all N-containing index pairs and values as their counts
n_index_pairs = {}

############### REVERSE COMPLEMENT FUNCTION ####################################

# Define function that takes a DNA sequece and returns the reverse complement
# Needed for orienting R3 reads according to how the unexpected_ and expected_index_pairs dictionaries were constructed
def rev_complement(seq):
    complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    rev = "".join([complements[base] for base in reversed(seq)])
    return rev

############### LOOP THROUGH FASTQ FILES #######################################

# Access the 2 FASTQ files in a loop simultaneously
# Run the quality filter on the index reads to determine if the indices pass
linecount = 0
for line2, line3 in zip(r2, r3):
    line2 = line2.strip("\n")
    line3 = line3.strip("\n")
    linecount += 1
    if linecount % 100000000 == 0:     # Progress report
        print("I am on line " + str(linecount) + ".", flush = True)
    elif linecount % 4 == 2:
        sequence2 = line2
        sequence3 = rev_complement(line3)     # Save sequence3 as the reverse complement of the actual read
    elif linecount % 4 == 0:
        quality2 = line2
        quality3 = line3

############### FILTER BY MIN_QUAL #############################################

# Filter index reads by quality filter with minimum threshold of min_qual
# This is a very stringent filter!
# If the quality of one base in either index read falls below the min_qual threshold,
# then the entire read is discarded and not counted in any dictionary
        scores2 = 0     # Set up a count for quality2 scores that pass the filter
        for char in quality2:
            phred2 = ord(char) - 33 # Convert ASCII score to phred score
            if phred2 >= min_qual:     # If the phred score is at or above the min_qual, increment scores counter
                scores2 +=1
            elif phred2 < min_qual:     # If the phred score is below the min_qual, reset scores to 0 and break
                scores2 = 0
                print("Uh-oh, line " + str(linecount) + " failed the quality2 filter.", flush = True)     # Progress report
                break
        scores3 = 0     # Set up a count for quality3 scores that pass the filter
        for char in quality3:
            phred3 = ord(char) - 33 # Convert ASCII score to phred score
            if phred3 >= min_qual:     # If the phred score is at or above the min_qual, increment scores counter
                scores3 += 1
            elif phred3 < min_qual:     # If the phred score is below the min_qual, reset scores to 0 and break
                scores3 = 0
                print("Uh-oh, line " + str(linecount) + " failed the quality3 filter.", flush = True)     # Progress report
                break

############### INCREMENT INDEX PAIRS DICTIONARIES #############################

# If all 8 nucleotides of both indices pass the the min_qual threshold, then oth score counters should at 8
# Test if the index pair is in the unexpected, expected, or N-contaiing dictionaries.
# If the index pair isn't in either any of these dictionaries, pass because it contains an odd, unexpected character charater
        if scores2 == 8 and scores3 == 8:
            new_index_pair = sequence2 + "_" + sequence3     # Make index pair that matches dictionary format
            if new_index_pair in expected_index_pairs:     # Is this new index pair expected?
                expected_index_pairs[new_index_pair] += 1     # Increment value
            elif new_index_pair in unexpected_index_pairs:     # Is this new index pair unexpected?
                unexpected_index_pairs[new_index_pair] += 1     # Increment value
            else:
                if "N" in new_index_pair:     # Search for Ns in the new index pair
                    if new_index_pair in n_index_pairs:     # If this index pair contains an N, count it or add it
                        n_index_pairs[new_index_pair] +=1     # Increment value
                    else:
                        n_index_pairs[new_index_pair] = 1     # Add it
                else:
                    pass
        else:
            pass

############### CLOSE FASTQ FILES ##############################################

# Close the FASTQ files
r2.close()
r3.close()
print("\n", flush = True)
print("I have closed your FASTQ files.", flush = True)     # Progress report
print("\n", flusqh = True)

############### WRITE TO OUTPUT FILE ###########################################

# Write dictionary counts to output with tab delimiter
# First column contains "expected" or "unexpected" to make grepping on the TSV file easier
print("I'm writing to your output file now!", flush = True)
with open(args.output_filepath, 'w') as output:
    output.write("@@@EXP_INDEX_PAIRS___PAIRS_ARE_ASSOCIATED_WITH_SAMPLES" + "\n")
    for k, v in expected_index_pairs.items():
        output.write("expected" + "\t" + str(k) + "\t" + str(v) + "\n")
    output.write("@@@UNEXP_INDEX_PAIRS___POTENTIAL_INDEX_SWAPPING" + "\n")
    for k, v in unexpected_index_pairs.items():
        output.write("unexpected" + "\t" + str(k) + "\t" + str(v) + "\n")
    output.write("@@@N_CONTAINING_INDEX_PAIRS" + "\n")
    for k, v in n_index_pairs.items():
        output.write("n_containing" + "\t" + str(k) + "\t" + str(v) + "\n")

############### CLOSING STATEMENT ##############################################

# Print ending text
print("\n", flush = True)
print("Guess what? I finished!", flush = True)
