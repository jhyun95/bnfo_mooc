#!/usr/bin/env python

import re

def find_clumped_kmers(seq,k,L,t):
    ''' Find kmers that "clump" in the sequence, defined as kmers such that 
        there exists an interval L with at least t occurences of the kmer'''
    print 'Counting Kmers...'
    kmerMap = kmer_map(seq,k); # maps kmers to lists of indices for occurence
    print len(kmerMap),'kmers counted'
    print 'Filtering Kmers...'
    kmerMap = dict(filter(lambda x: len(x[1]) >= t, kmerMap.items())) # filter for kmers with >=t occurences
    print len(kmerMap),'kmers left'
    clumpedKmers = []
    print 'Screening Kmers...'
    for kmer in kmerMap:
        indices = kmerMap[kmer]
        for ii in xrange(len(indices)-t+1):
            indexSubset = indices[ii:ii+t]
            first = indexSubset[0]
            last = indexSubset[-1]
            intervalLength = last - first + k
            if intervalLength <= L:
                clumpedKmers.append(kmer); break
    return clumpedKmers

def reverse_complement(seq):
    ''' Returns the reverse complement of a DNA strand '''
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','a':'t','t':'a','c':'g','g':'c'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def kmer_map(text,k):
    ''' Returns a dictionary of kmer occurences 0-indexed for <text> {kmer:[indices]} '''
    kmers = {}
    for i in xrange(len(text)-k+1):
        kmer = text[i:i+k]
        if kmer in kmers:
            kmers[kmer].append(i)
        else:
            kmers[kmer] = [i]
    return kmers

def kmer_frequencies(text,k):
    ''' Returns a dictionray of kmer counts for <text> {kmer:count} '''
    kmers = kmer_map(text,k); counts = {}
    for kmer in kmers:
        counts[kmer] = len(kmers[kmer])
    return counts

def most_frequent_kmers(text,k):
    ''' Returns a list of the most frequently occuring kmers in <text> '''
    kmers = kmer_frequencies(text,k) 
    maxCount = max(kmers.values())
    maxKmers = filter(lambda x: kmers[x] >= maxCount, kmers.keys())
    return maxKmers

def pattern_matches(text,pattern):
    ''' Returns the indices at which <pattern> occurs in <text>
        Includes overlapping occurences separately, 0-indexed '''		
    rePattern = '(?=' + pattern + ')'
    match = [m.start() for m in re.finditer(rePattern,text)]
    return match

def pattern_count(text,pattern):
    ''' Returns the number of time <pattern> occurs in <text>
        Includes overlapping occurances separately '''
    return len(pattern_matches(text,pattern))
