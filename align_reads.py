'''
align_reads.py
short read mapping analysis over test data
'''

import sys
import random
from search_bwt import *

# Get each read's identifier, position, and string value.
# myfile is the file of reads in FASTA format.
def parse_reads(myfile):
    file_name = myfile
    read_file = open(file_name, 'r')

    # read_dict is a dictionary with: 
    # key = unique read identifier
    # value = (true position, read string) tuple
    read_dict = {}
    identifier = ''
    position = 0
    read = ''

    line_num = 0

    for line in read_file:

        line_num += 1

        # Odd numbered lines are read metadata
        # Metadata contains unique read identifer and true position
        if line_num % 2 == 1:
            word_arr = line.split()
            identifier = word_arr[0]
            position = int((word_arr[3])[4:])

        # Even numbered lines are read strings (actual nucleotide sequence)
        if line_num % 2 == 0:
            # Pick random bases to replace unknown nucleotides (rare)
            base = ''
            rand = random.randint(1,4)
            if rand == 1:
                base = 'A'
            elif rand == 2:
                base = 'T'
            elif rand == 3:
                base = 'G'
            else:
                base = 'C'

            # Add the (key, value) pair to the dictionary of reads
            # Remove newline characters and replace unknown nucleotides in sequences
            read_dict[identifier] = (position, line.rstrip().replace('N',base))

    read_file.close()
    return read_dict


# Calculates the reverse compliment of a nucleotide sequence.
# Needed because reads are pair-ended.
def reverse_complement(base_string):
    base_dict = { 'A':'T', 'T':'A', 'C':'G', 'G':'C' }
    return "".join([base_dict[base] for base in reversed(base_string)])  


# Aligns each read to the reference genome.
def align_reads(genome_file, reads_file, threshold):
    print "Calculating BWT, BWT reverse, suffix array, and aligning reads..."

    # Get read dictionary
    read_dict = parse_reads(reads_file)

    # Get reference genome from file
    fref = open(genome_file)
    ref = ''.join(fref.readlines()).replace('\n','')

    # Calculate suffix array, BWT, and reverse BWT of genome
    sa = suffix_array(ref)
    bw = bwt(ref)
    bwr = bwt(ref[::-1])

    # Search over all reads

    # New dictionary of reads with: 
    # key = unique read identifier
    # value = predicted best position in genome
    new_read_dict = {}

    # Iterate through every read
    for identifier in read_dict:
        # s is the read string
        s = read_dict[identifier][1]

        # Determine the best position of the read, and the reverse complement of the read.
        # Need both because data is pair-ended, and you don't know which direction is correct.
        best_position1, score1 = best_match_position(bw, bwr, s, threshold, sa)
        best_position2, score2 = best_match_position(bw, bwr, reverse_complement(s), threshold, sa)

        # Determine which position is better (higher score)
        if score1 >= score2 and score1 != -1:
            new_read_dict[identifier] = best_position1
        elif score1 < score2:
            new_read_dict[identifier] = best_position2
        elif score1 == score2 and score1 == -1:
            new_read_dict[identifier] = -1              
        else:
            new_read_dict[identifier] = -2

    num_correct = 0
    num_incorrect = 0
    num_missing = 0
    num_error = 0

    # Calculating results by comparing predicted positions with true positions
    for identifier in read_dict:
        if read_dict[identifier][0] == new_read_dict[identifier]:
            num_correct += 1
        elif new_read_dict[identifier] == -1:
            num_missing += 1
        elif new_read_dict[identifier] == -2:
            num_error += 1
        elif read_dict[identifier][0] != new_read_dict[identifier]:
            num_incorrect += 1

    print "\nResults (threshold=" + str(threshold) + ")..."
    print "Number of reads: \t\t\t" + str(len(read_dict))
    print "Number of correct alignments: \t\t" + str(num_correct)
    print "Number of 'no matches': \t\t" + str(num_missing)
    print "Number of incorrect predictions: \t" + str(num_incorrect)
    print "Number of read errors: \t\t\t" + str(num_error)


def main():
    # Default threshold is 3
    threshold = 3

    if "-t" in sys.argv and len(sys.argv) == 5:
        threshold = float(sys.argv[4])
    elif "-t" not in sys.argv and len(sys.argv) == 3:
        pass
    else:
        print '\nusage: python align_reads.py <genome file name> <read file name> [-t <threshold level>]\n'
        return

    align_reads(sys.argv[1], sys.argv[2], threshold)
    
main()
