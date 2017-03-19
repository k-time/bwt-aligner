"""
bwt.py
Implementation of burrows-wheeler transform for genome data
"""

alphabet = ['A', 'C', 'G', 'T']


def suffix_array(s):
    """build suffix array of string s"""
    sa = sorted([(s[i:], i) for i in xrange(0, len(s)+1)])
    return map(lambda x: x[1], sa)


def bwt(t):
    """compute burrows-wheeler transform of string t"""
    bw = []
    for i in suffix_array(t):
        if i == 0:
            bw.append('$')
        else:
            bw.append(t[i-1])
    return ''.join(bw)


def first_col(totals):
    """make dict of chars to range of occurrences in first column"""
    col = {}
    temp = 0
    for i, j in sorted(totals.iteritems()):
        col[i] = (temp, temp + j)
        temp += j
    return col


def ibwt(bw):
    """decode bwt"""
    ranks, totals = rank(bw)
    fc = first_col(totals)
    row = 0
    t = '$'

    while bw[row] != '$':
        char = bw[row]
        t = char+t
        row = fc[char][0] + ranks[row]

    return t


def rank(bw):
    """rank(char) := list of number of occurrences of a char for each substring R[:i] (reference)"""
    totals = {}
    ranks = {}

    for char in alphabet:
        if (char not in totals) and (char != '$'):
            totals[char] = 0
            ranks[char] = []

    for char in bw:
        if char != '$':
            totals[char] += 1
        for t in totals.iterkeys():
            ranks[t].append(totals[t])

    return ranks, totals


def count_matches_exact(bw, s):
    """return number of exact matches of s to bw"""
    ranks, totals = rank(bw)
    fc = first_col(totals)

    # char isn't in bwt
    if s[-1] not in fc:
        return 0

    l, r = fc[s[-1]]
    i = len(s)-2
    while i >= 0 and r > 1:
        char = s[i]
        l = fc[char][0] + ranks[char][l-1]  # R(aW)
        r = fc[char][0] + ranks[char][r-1]  # Rbar(aW)
        i -= 1
        print('l: '+str(l)+' r: '+str(r))
        print(''.join(bw[l:r]))
        print(''.join(sorted(bw)[l:r]))

    return r-l
