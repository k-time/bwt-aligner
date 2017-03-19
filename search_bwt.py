"""
search_bwt.py
command-line tool for short read mapping
inexact search over burrows-wheeler genome data
"""

import sys
from enum import Enum
from operator import itemgetter
import time
from bwt import *

# bw is bwt of genome (B)
# bwr is bwt of reverse genome (B')
# s is short read to be matched (W)
# diff is max num of differences (z)

# alphabet of symbols in allowed
alphabet = {'A', 'C', 'G', 'T'}

# defined below
C = {}
O = {}
D = []

# enum represents the type of the alignment choice
Type = Enum('Type', 'START MATCH MISMATCH INSERTION DELETION')

# rewards/penalties
gap_open = 3
gap_ext = 1
mismatch = 1
match = 0

# option switches
NO_INDELS = False
sub_mat = {}

num_prunes = 0


def inexact_search(bw, bwr, s, diff):
    """find suffix array intervals with up to diff differences"""

    global C, O, D, num_prunes
    # totals, ranks
    # O is a dictionary with keys $,A,C,G,T, and values are arrays of counts
    O, tot = rank(bw)

    # reverse ranks
    Oprime, junk = rank(bwr)

    # C[a] := number of lexicographically smaller letters than a in bw/reference
    C = compute_C(tot)

    # D[i] := lower bound on number of differences in substring s[1:i]
    D = compute_D(s, C, Oprime, bw)

    # call the recursive search function and return a list of SA-range tuples
    sa_index_set = inexact_recursion(s, len(s)-1, diff, 0, len(bw)-1, Type.START)
    index_dict = {}

    for (i, j) in sa_index_set:
        # if index already exists, pick the higher diff value
        if i in index_dict:
            if index_dict[i] < j:
                index_dict[i] = j
                num_prunes += 1

        else:
            index_dict[i] = j

    # sort list by diff from highest to lowest
    return sorted(index_dict.items(), key=itemgetter(1), reverse=True) 


def best_match_position(bw, bwr, s, diff, sa):
    sa_index_list = inexact_search(bw, bwr, s, diff)
    if len(sa_index_list) != 0:
        best_index, score = sa_index_list[0]
        return sa[best_index]+1, score
    else:
        return -1, -1


def compute_C(totals):
    """compute C, the number of lexicographically greater symbols in the ref"""
    C = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for k in alphabet:
        for ref in alphabet:
            if ref < k:
                C[k] += totals[ref]

    return C


def compute_D(s, C, Oprime, bw):
    """compute estimated lower bounds of differences in substring s[0:i] for all  in [0,len(s)]"""
    k = 1
    l = len(bw)-2
    z = 0
    D = [0] * len(s)

    for i in range(0, len(s)):
        k = C[s[i]] + Oprime[s[i]][k-1] + 1
        l = C[s[i]] + Oprime[s[i]][l]
        if k > l:
            k = 1
            l = len(bw)-1
            z += 1
        D[i] = z

    return D


def get_D(i):
    """enforce condition that if D[i] is set to -1, its value will be considered as 0"""
    if i < 0:
        return 0
    else:
        return D[i]


def get_O(char, index):
    """see get_D()"""
    if index < 0:
        return 0
    else:
        return O[char][index]


def inexact_recursion(s, i, diff, k, l, prev_type):
    """search bwt recursively and tolerate errors"""
    
    global num_prunes

    # pruning based on estimated mistakes
    if diff < get_D(i):
        num_prunes += 1
        return set()

    # end of query condition
    temp = set()
    if i < 0:
        for j in range(k, l+1):
            temp.add((j, diff))
        return temp

    # search
    sa_idx = set()  # set of suffix array indices at which a match starts
    
    if not NO_INDELS:
        # Insertion
        if prev_type == Type.INSERTION:
            sa_idx = sa_idx.union(inexact_recursion(s, i-1, diff-gap_ext, k, l, Type.INSERTION))
        else:
            sa_idx = sa_idx.union(inexact_recursion(s, i-1, diff-gap_ext-gap_open, k, l, Type.INSERTION))

    for char in alphabet:
        temp_k = C[char] + get_O(char, k-1) + 1
        temp_l = C[char] + get_O(char, l)

        if temp_k <= temp_l:
            if not NO_INDELS:
                # Deletion
                if prev_type == Type.DELETION:
                    sa_idx = sa_idx.union(inexact_recursion(s, i, diff-gap_ext, temp_k, temp_l, Type.DELETION))
                else:
                    sa_idx = sa_idx.union(inexact_recursion(s, i, diff-gap_ext-gap_open, temp_k, temp_l, Type.DELETION))
            if char == s[i]:
                # Match!
                sa_idx = sa_idx.union(inexact_recursion(s, i-1, diff+match, temp_k, temp_l, Type.MATCH))
                
            else:
                # Mismatch
                if sub_mat:
                    sa_idx = sa_idx.union(inexact_recursion(s, i-1, diff-mismatch*sub_mat[(s[i], char)],
                                                            temp_k, temp_l, Type.MISMATCH))
                else:
                    sa_idx = sa_idx.union(inexact_recursion(s, i-1, diff-mismatch, temp_k, temp_l, Type.MISMATCH))

    return sa_idx


def estimate_substitution_mat(ref, r):
    """get likelihood of each substitution type over all possible alignments"""
    mismatches = {}

    for i in range(0, len(ref)):
        for j in range(0, len(r)):
            if ref[i] != r[j]:
                if (ref[i], r[j]) in mismatches:
                    mismatches[(ref[i], r[j])] += 1
                else:
                    mismatches[(ref[i], r[j])] = 1

    scale = max(mismatches.values())
    for k in mismatches:
        mismatches[k] = float(mismatches[k])/scale

    return mismatches


def print_output(sa_index_list, sa, s, read):
    """print formatted output"""
    sa_values = [(sa[i], j) for (i, j) in sa_index_list]

    print '\n-------------------------------------'
    print 'Reference: ' + s
    print 'Read: \t   ' + read + '\n'

    print str(len(sa_values)) + " match(es) found!\n"
    print "Score\tPos.\tSuffix\n"
    for v, x in sa_values:
        print str(x) + "\t" + str(v) + "\t" + s[v:v+35]

    print '-------------------------------------'


def test():
    # Another test: 'ATGCGTAATGCCGTCGATCG'
    s = 'CGATCCGCGCTGCTGATGATCGATG'
    read = 'GATGAT'
    threshold = 2

    sa = suffix_array(s)
    bw = bwt(s)
    bwr = bwt(s[::-1])

    print_output(inexact_search(bw, bwr, read, threshold), sa, s, read)


def main():
    start = time.time()

    threshold = 3  # this will be the z value
    usage = '\nusage: python search_bwt.py [--no-indels] [test|<reference file name>] [<read file name>]\n'

    if '--no-indels' in sys.argv:
        global NO_INDELS
        NO_INDELS = True

    if '--linear-gaps' in sys.argv:
        global gap_open, gap_ext
        gap_open = 0
        gap_ext = 1

    if len(sys.argv) == 1:
        print usage
        return

    elif sys.argv[1].lower() == 'test':
        test()
        return
    
    elif len(sys.argv) < 3:
        print usage
        return

    for i in range(0, len(sys.argv)):
        if sys.argv[i] == '-t' and i < len(sys.argv)-2:
            threshold = int(sys.argv[i+1])

    fread = open(sys.argv[-1])
    fref = open(sys.argv[-2])

    ref = ''.join(fref.readlines()).replace('\n', '')
    read = ''.join(fread.readlines()).replace('\n', '')

    # estimate the substitution matrix
    if '--no-sub-mat' not in sys.argv:
        global sub_mat
        sub_mat = estimate_substitution_mat(ref, read)

    sa = suffix_array(ref)

    bw = bwt(ref)
    bwr = bwt(ref[::-1])

    print read
    print_output(inexact_search(bw, bwr, read, threshold), sa, ref)

    elapsed = time.time() - start
    if '--show-time' in sys.argv:
        print 'time elapsed: '+str(elapsed)
    if '--count-prunes' in sys.argv:
        print str(num_prunes) + ' nodes pruned.'
    print "error score upper bound: " + str(threshold)
    fread.close()
    fref.close()


if __name__ == "__main__":
    main()
