# -*- coding: utf-8 -*-
"""
Created on Wed Aug 09 23:28:35 2017

@author: jhyun_000
"""

import numpy as np
codonTablePath = "Codons.txt"
aamwPath = "AAMWs.txt"

def main():
    ''' Refer to http://rosalind.info/problems/list-view/ '''
    problemKMF()
    
def problemKMF(filename='Data/rosalind_kmp.txt'):
    header, seq = parse_fasta(filename).items()[0]
    n = len(seq) # sequence length
    failure = [0] * n # failure array
    j = 1 # substring check length
    for k in xrange(1,n): # search index
        suffix = seq[k-j+1:k+1]
        prefix = seq[:j]
        if suffix == prefix: # if matching at this length, scan forward
            while suffix == prefix and j <= k:
                j += 1
                suffix = seq[k-j+1:k+1]
                prefix = seq[:j]
            j -= 1 # overcounted
        else: # if not matching at this length, scan downward
            while suffix != prefix:
                j -= 1
                suffix = seq[k-j+1:k+1]
                prefix = seq[:j]
        failure[k] = j
        j += 1
    print ' '.join(map(str,failure))
    return failure
    
def problemKMER(filename='Data/rosalind_kmer.txt'):
    import re
    def overlapping_count(search, query):
        return sum(1 for m in re.finditer('(?=%s)' % query, search))
    header,seq = parse_fasta(filename).items()[0]
    kmers = make_kmers(4)
    counts = map(lambda kmer: overlapping_count(seq,kmer),kmers)
    print ' '.join(map(str,counts))
    
def problemINOD(n=4):
    return n-2
    
def problemCORR(filename='Data/rosalind_corr.txt'):
    data = parse_fasta(filename)
    seqs = data.values()
    counts = {}
    for seq in seqs: # Make read counts
        rseq = revc(seq)
        if seq in counts:
            counts[seq] += 1
        elif rseq in counts:
            counts[rseq] += 1
        else:
            counts[seq] = 1
    verified = {}; errors = []
    for seq in counts: # Separate into verified and suspect reads
        if counts[seq] > 1:
            verified[seq] = counts[seq]
        else:
            errors.append(seq)
    for seq in errors: # Find closest match to correct suspect reads
        for correct in verified:
            rcorr = revc(correct)
            if hamming_dist(seq,correct) == 1:
                print seq + "->" + correct; break
            elif hamming_dist(seq,rcorr) == 1:
                print seq + "->" + rcorr; break
    
def problemCAT(filename='Data/rosalind_cat.txt'):
    header,rna = parse_fasta(filename).items()[0]
    matching = {'A':'U','U':'A','C':'G','G':'C'}
    solved = {}
    def is_balanced(seq):
        return seq.count('A') == seq.count('U') and seq.count('C') == seq.count('G')
    def noncrossing(seq):
        if len(seq) == 0:
            return 1
        elif seq in solved:
            return solved[seq]
        end = seq[0]; target = matching[end]
        matchings = 0
        for i in xrange(1,len(seq)):
            if seq[i] == target: 
                split1 = seq[1:i]; split2 = seq[i+1:]
                if is_balanced(split1) and is_balanced(split2):
                    matchings += noncrossing(split1) * noncrossing(split2)
        solved[seq] = matchings
        return matchings
    print noncrossing(rna) % 1000000                 

def problemTREE(filename='Data/rosalind_tree.txt'):
    with file(filename,'r+') as f:
        n = int(f.readline().strip())
        edges = []
        line = f.readline().strip()
        while line:
            a,b = map(int,line.split())
            edges.append((a,b))
            line = f.readline()
        f.close()
    subtrees = []; orphans = range(1,n+1)
    for edge in edges:
        a,b = edge; subtreeA = None; subtreeB = None
        for subtree in subtrees:
            if a in subtree:
                subtreeA = subtree
            if b in subtree:
                subtreeB = subtree
        if subtreeA != None and subtreeB != None:
            subtrees.remove(subtreeA)
            subtrees.remove(subtreeB)
            subtrees.append(subtreeA + subtreeB)
        elif subtreeA != None and subtreeB == None:
            subtreeA.append(b)
        elif subtreeA == None and subtreeB != None:
            subtreeB.append(a)
        else:
            subtrees.append([a,b])
        if a in orphans:
            orphans.remove(a)
        if b in orphans:
            orphans.remove(b)
    for orphan in orphans:
        subtrees.append([orphan])
    print subtrees
    print len(subtrees)-1
    
def problemTRAN(filename='Data/rosalind_tran.txt'):
    seq1,seq2 = parse_fasta(filename).values()
    transitions = {'A':'G','G':'A','C':'T','T':'C'}
    tranversions = {'A':['C','T'],'G':['C','T'],'C':['A','G'],'T':['A','G']}
    trn = 0.0; trv = 0.0
    for i in xrange(len(seq1)):
        bp1 = seq1[i]; bp2 = seq2[i]
        if bp2 == transitions[bp1]:
            trn += 1.0
        elif bp2 in tranversions[bp1]:
            trv += 1.0
    print trn/trv           
    
def problemSSEQ(filename='Data/rosalind_sseq.txt'):
    data = parse_fasta(filename)
    seq = data.values()[0]; motif = data.values()[1]
    if len(motif) > len(seq):
        temp = seq; seq = motif; motif = temp
    indices = []    
    
    n1 = len(seq); n2 = len(motif); i = 0
    for j in xrange(n2):
        query = motif[j]
        while query != seq[i] and i < n1:
            i += 1
        indices.append(i+1); i += 1

    assembled = ""    
    for i in indices:
        assembled += seq[i-1]
    print ' '.join(map(str,indices))      
    
def problemSIGN(n=2,fout='output_sign.txt'):
    def make_sign_permutations(nums):
        if len(nums) == 1:
            return [[nums[0]], [-nums[0]]]
        signPermutations = []
        for subPerm in make_sign_permutations(nums[1:]):
            signPermutations.append([nums[0]] + subPerm)
            signPermutations.append([-nums[0]] + subPerm)
        return signPermutations
    
    count = 2**n * factorial(n); output = ""
    print count; output += str(count) + '\n'
    unsigned = permutations(range(1,n+1)); signed = []
    for perm in unsigned:
        signed += make_sign_permutations(perm)
    for perm in signed:
        line = ' '.join(map(str,perm))
        print line; output += line + '\n'

    f = file(fout,'w+')
    f.write(output.strip())
    f.close()
    
def problemPROB(seq='ACGATACAA',gcContents='0.129 0.287 0.423 0.476 0.641 0.742 0.783'):
    gcs = map(float,gcContents.split())
    logProbs = []
    for gc in gcs:
        logGC = np.log10(gc/2.0); logAT = np.log10(0.5-gc/2.0)
        weights = {'A':logAT, 'T':logAT, 'G':logGC, 'C':logGC}
        logProb = sum(map(lambda bp: weights[bp], list(seq)))
        logProbs.append(logProb)
    fmted = map(lambda x: round(x,3), logProbs)
    print ' '.join(map(str,fmted))

def problemPPER(n,k):
    print reduce(lambda x,y: (x*y)%1000000, range(n-k+1,n+1))
    
def problemPMCH(seq='AGCUAGUCAU'):
    # Assumes AU and CG are balanced for perfect matchings
    AUpairs = seq.count('A')
    CGpairs = seq.count('C')
    print factorial(AUpairs)*factorial(CGpairs)
    
def problemLONG(filename='Data/rosalind_long.txt'):
    def fuse(seq1,seq2,minOverlap):
        # Attempt seq1-seq2 overlap, then seq2-seq1 overlap
        over12 = 0; over21 = 0
        for i in xrange(len(seq1)):
            if seq1[-i:] == seq2[:i]:
                over12 = i
            if seq2[-i:] == seq1[:i]:
                over21 = i
        # Attempt to fuse
        if over12 >= minOverlap: # seq1-seq2
            return seq1 + seq2[over12:]
        elif over21 >= minOverlap: # seq2-seq1
            return seq2 + seq1[over21:]
        else:
            return None
        
    def fuse_pass(frags,minOverlap):
        contigs = []
        for frag in frags:
            newContig = None
            for contig in contigs:
                merged = fuse(contig,frag,minOverlap)
                if merged:
                    newContig = merged; break
            if newContig:
                contigs.remove(contig)
                contigs.append(newContig)
            else:
                contigs.append(frag)
        return contigs
        
    data = parse_fasta(filename)
    frags = data.values()
    minOverlap = (len(frags[0])+1)/2
    passes = 0
    print "Pass",passes,"\tFragments",len(frags)
    while len(frags) > 1:
        passes += 1
        frags = fuse_pass(frags,minOverlap)
        print "Pass",passes,"\tFragments",len(frags)        
    print frags[0]
    
def problemLGIS(filename='Data/rosalind_lgis_sample.txt'):
    f = file(filename,'r+')
    n = int(f.readline().strip())
    perm = f.readline().strip()
    f.close()
    
    objs = map(int,perm.split())
    def find_longest_path(decreasing=True):
        dist = [] # distance of node at given index
        maxParents = {} # parent that goes along longest path
        for i in xrange(n):
            parents = [] # indices of all parent nodes
            for j in xrange(i):
                if decreasing and objs[j] > objs[i]: # decreasing
                    parents.append(j)
                if not decreasing and objs[j] < objs[i]: # increasing
                    parents.append(j)
            maxDist = -1; maxParent = None
            for parent in parents:
                if dist[parent] > maxDist:
                    maxParent = parent; maxDist = dist[parent]
            dist.append(maxDist+1); maxParents[i] = maxParent
        endpoint = dist.index(max(dist))
        path = [endpoint]; parent = maxParents[endpoint]
        while not parent == None:
            path.append(parent)
            parent = maxParents[parent]
        path.reverse()
        objPath = map(lambda x: objs[x],path)
        return objPath
      
    ipath = find_longest_path(False)          
    dpath = find_longest_path(True)    
    print ' '.join(map(str,ipath))
    print ' '.join(map(str,dpath))
    
def problemLEXF(n=2,alpha='A C G T'):
    syms = alpha.split()
    for kmer in make_kmers(n,syms):
        print kmer
    
def problemSPLC(filename="Data/rosalind_splc.txt"):
    f = file(filename,'r+')
    dnaHeader = f.readline().strip()[1:]
    f.close()
    
    introns = parse_fasta(filename)
    dna = introns[dnaHeader]
    del introns[dnaHeader]
    
    for intron in introns.values():
        dna = ''.join(dna.split(intron))
    print translate(dna.replace('T','U'))
    
def problemREVP(filename="Data/rosalind_revp.txt"):
    data = parse_fasta(filename)
    header,seq = data.items()[0]
    n = len(seq); rmin = 4; rmax = 12
    for i in xrange(n):
        for r in xrange(rmin,rmax+1,2):
            test = seq[i:i+r]
            s1 = test[:r/2]; s2 = test[r/2:]
            if s1 == revc(s2):
                print i+1,r
    
def problemPRTM(prot='SKADYEK'):
    aamws = get_aa_mw_table()
    print sum(map(lambda x: aamws[x], list(prot)))
    
def problemPERM(n=5):
    
    perms = permutations(range(1,n+1))
    print len(perms)
    for perm in perms:
        print ' '.join(map(str,perm))
            
    
def problemORF(filename='Data/rosalind_orf.txt'):
    def find_orf_1way(seq):
        orfs = []
        for i in xrange(len(seq)-2):
            if seq[i:i+3] == 'AUG':
                aa = 'M'; j = i+3; codon = seq[j:j+3]
                while codon in codonTable and not 'Stop' in aa:
                    aa += codonTable[codon]
                    j += 3; codon = seq[j:j+3]
                if 'Stop' in aa:
                    orfs.append(aa[:-4])
            i += 1
        return orfs
    header, seq = parse_fasta(filename).items()[0]
    seq1 = seq.replace('T','U')
    seq2 = problemREVC(seq).replace('T','U')
    codonTable = get_codon_table()
    orfs1 = find_orf_1way(seq1)
    orfs2 = find_orf_1way(seq2)
    orfs = {}
    for orf in orfs1 + orfs2:
        orfs[orf] = True
    for orf in orfs:
        print orf
    
def problemMRNA(seq='MA'):
    codonTable = get_codon_table()
    redundancy = {}
    for codon in codonTable:
        aa = codonTable[codon]
        if aa in redundancy:
            redundancy[aa] += 1
        else:
            redundancy[aa] = 1
    aas = list(seq) + ['Stop']
    seqRedundancy = map(lambda aa: redundancy[aa], aas)
    print reduce(lambda x,y: x*y % 1000000, seqRedundancy)
    
def problemMPRT(filename='Data/rosalind_mprt.txt'):
    import urllib2, re
    def parse_fna(text):
        data = text.strip().split("\n")
        header = data[0][1:]
        seq = ''.join(data[1:])
        return header,seq    
        
    pattern = '(?=N[^P][ST][^P])'
    for line in file(filename,'r+'):
        ID = line.strip()
        url = 'http://www.uniprot.org/uniprot/' + ID + '.fasta'
        try:
            response = urllib2.urlopen(url)
            header, seq = parse_fna(response.read())
        except urllib2.HTTPError:
            print 'TRYING:',url
            url = 'http://www.uniprot.org/uniprot/' + ID.split('_')[0] + '.fasta'
            response = urllib2.urlopen(url)
            header, seq = parse_fna(response.read())
#        print ID, [seq[m.start():m.end()] for m in re.finditer(pattern, seq)]
        indices = [m.start()+1 for m in re.finditer(pattern, seq)]
        if len(indices) > 0:
            print ID
            print ' '.join(map(str,indices))
    
def problemLIA(k,N):
    from scipy.stats import binom
    print 1.0 - binom.cdf(N-1,2**k,0.25)
#    print np.dot(distr,)
    
def problemLCSM(filename='Data/rosalind_lcsm.txt'):
    ''' TODO: Try Rabin-Karp algorithm, there are many faster ways to do this,
        i.e. binary search through first element and check all '''
    def find_lsms(seq1,seq2):
        query,check = (seq1,seq2) if len(seq1) < len(seq2) else (seq2,seq1)
        lsms = [""]; inCheck = {}
        i = 0; k = 1; # i = search index, k = search length
        while i + k <= len(query):
            test = query[i:i+k]
            if not test in inCheck:
                inCheck[test] = test in check
            if inCheck[test]:
                if len(test) > len(lsms[0]):
                    lsms = [test]
                elif len(test) == len(lsms[0]):
                    lsms.append(test)
                k += 1
            else:
                i += 1; k = len(lsms[0])
        return lsms
                
    seqs = parse_fasta(filename).values()
    seqs = sorted(seqs,cmp=lambda x,y: len(x)-len(y))
    
    lsms = find_lsms(seqs[0],seqs[1]); i = 2
    while i < len(seqs):
        submotifs = []
        for motif in lsms:
            submotifs += find_lsms(motif,seqs[i])
        longest = max(map(len,submotifs))
        lsms = filter(lambda x: len(x) == longest,submotifs)
        i += 1
    print lsms
    print lsms[0]

def problemIEV(counts="1 0 0 1 0 1"):
    data = map(int,counts.strip().split())
    print 2*(data[0] + data[1] + data[2] + 0.75*data[3] + 0.5*data[4])
    
def problemGRPH(filename='Data/rosalind_grph (3).txt'):
    def ends_only(seq,e=3):
        return (seq[:e],seq[-e:])
    data = parse_fasta(filename,ends_only)
    for entry in data.items():
        print entry
    counter = 0
    for ID1 in data.keys():
        for ID2 in data.keys():
            if ID1 != ID2 and data[ID1][1] == data[ID2][0]:
                print ID1+" "+ID2#,data[ID1],data[ID2]
            counter += 1 
    
def problemFIBD(n=7,m=3):
    from decimal import Decimal, getcontext
    getcontext().prec = 50
    F = [Decimal(0)] * n; A = [Decimal(0)] * m
    t = 1; A[0] = Decimal(1); F[0] = Decimal(1)
    while t < n:
        Ao_next = sum(A[1:])
        Anext = [Ao_next] + A[:-1]
        A = Anext; t += 1
        F[t-1] = sum(A)
    print long(F[-1])

def problemCONS(filename='Data/rosalind_cons.txt'):
    data = parse_fasta(filename)
    n = len(data.values()[0])
    profile = [[0]*n,[0]*n,[0]*n,[0]*n]
    baseIndex = {'A':0,'C':1,'G':2,'T':3}
    for seq in data.values():
        for i in xrange(n):
            base = seq[i]; j = baseIndex[base]
            profile[j][i] += 1
            
    cons = ""
    for i in xrange(n):
        baseProf = {'A':profile[0][i],'C':profile[1][i],'G':profile[2][i],'T':profile[3][i]}
        baseItems = baseProf.items()
        baseItems = sorted(baseItems,lambda x,y: y[1]-x[1])
        cons += baseItems[0][0]
    
    print cons
    print 'A:', ' '.join(map(str,profile[0]))
    print 'C:', ' '.join(map(str,profile[1]))
    print 'G:', ' '.join(map(str,profile[2]))
    print 'T:', ' '.join(map(str,profile[3]))             

def problemSUBS(seq,motif):
    n = len(seq); k = len(motif); hits = []
    for i in range(k,n):
        kmer = seq[i-k:i]
        if kmer == motif:
            hits.append(i-k+1)
    print ' '.join(map(str,hits))
    return hits
    
def problemPROT(rna):
    prot = ""; i = 0; codon = rna[i:i+3]; 
    codonTable = get_codon_table()
    while len(codon) == 3 and codon in codonTable and codonTable[codon] != 'Stop':
        prot += codonTable[codon]
        i += 3; codon = rna[i:i+3]
    return prot

def problemIPRB(k,m,n):
    k = float(k); m = float(m); n = float(n); s = k+m+n
    rec = (n/s)*(n-1)/(s-1) + (n/s)*m/(s-1) + (m/s)*(m-1)/(s-1)*0.25
    return 1 - rec

def problemHAMM(seq1,seq2):
    return sum(map(lambda x: seq1[x] != seq2[x], xrange(min(len(seq1),len(seq2)))))
    
def problemGC(filename='Data/rosalind_gc.txt'):
    def gc_count(seq):
        return float(seq.count('G')+seq.count('C')) / len(seq) * 100
        
    data = parse_fasta(filename)
    maxGCheader = max(data.keys(), key=lambda x: gc_count(data[x]))
    print maxGCheader
    print gc_count(data[maxGCheader])
    
def problemFIB(n,k):
    Fn_2 = 1; Fn_1 = 1
    for x in xrange(n-2):
        Fn = Fn_1 + Fn_2*k
        Fn_2 = Fn_1
        Fn_1 = Fn
    return Fn

def problemREVC(seq):
    rc = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(map(lambda x: rc[x], seq[::-1]))

def problemRNA(seq):
    return seq.replace('T','U')

def problemDNA(seq):
    seq = seq.upper()
    return seq.count('A'), seq.count('C'), seq.count('G'), seq.count('T')
    
''' Utilities '''

def revc(seq):
    ''' Calculates the reverse complement of a strand of DNA '''
    rc = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join(map(lambda x: rc[x], seq[::-1]))
    
def translate(rna):
    ''' Translates a strand of RNA to an AA string '''
    prot = ""; i = 0; codon = rna[i:i+3]; 
    codonTable = get_codon_table()
    while len(codon) == 3 and codon in codonTable and codonTable[codon] != 'Stop':
        prot += codonTable[codon]
        i += 3; codon = rna[i:i+3]
    return prot
    
def make_kmers(n,syms=['A','C','G','T']):
    ''' Construct all kmers length n using symbols in syms (default ACGT)
        Warning, memory-intensive, needs to rework to make generator '''
    def perm(n):
        if n == 1:
            return syms
        permutations = []
        substrings = perm(n-1)
        for sym in syms:
            for substring in substrings:
                permutations.append(sym + substring)
        return permutations
    return perm(n)
    
def hamming_dist(seq1,seq2):
    ''' Find hamming distance between two seqs, from problem HAMM '''
    return sum(map(lambda x: seq1[x] != seq2[x], xrange(min(len(seq1),len(seq2)))))

def get_aa_mw_table():
    f = file(aamwPath,'r+'); aamws = {}
    for line in f:
        aa,mw = line.strip().split()
        aamws[aa] = float(mw)
    return aamws    
    
def get_codon_table():
    f = file(codonTablePath,'r+'); codonTable = {}
    for line in f:
        data = line.strip().split()
        for i in range(0,len(data),2):
            codon = data[i]; aa = data[i+1]
            codonTable[codon] = aa
    return codonTable

def parse_fasta(filename,process=None):
    ''' Returns {Header:Sequence}, apply function "process" to each seq '''
    f = file(filename,'r+'); header = ""; data = {}
    for line in f: # Calculate GC content of each entry
        if line[0] == ">":
            if header != "" and process: # new line, post process the string
                data[header] = process(data[header])
            header = line[1:].strip(); data[header] = ""
        else:
            data[header] += line.strip()
    if process:
        data[header] = process(data[header])
    return data
    
def factorial(n,mod=0):
    if mod == 0:
        return reduce(lambda x,y:x*y, range(1,n+1))
    else:
        return reduce(lambda x,y:(x*y) % mod, range(1,n+1)) 
    
def permutations(objs,seed=[]):
    subs = []
    if len(objs) == 0: # end of permtation, return only result
        return [seed]
    for obj in objs: # cycle through all remaining possible objects
        nextSeed = seed + [obj]
        remainingObjs = objs+[]; remainingObjs.remove(obj)
        subs.append((remainingObjs,nextSeed))
    subSeqs = map(lambda x: permutations(x[0],x[1]), subs)
    subSeqs = reduce(lambda x,y: x+y,subSeqs) 
    return subSeqs
    
main()