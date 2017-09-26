#!/usr/bin/env python

import bnfo1

def main():
    problem1_2_10()

def problem1_3_2():
    seq = line_to_params('Data/bnfo1_week1/dataset_3_2.txt')[0]
    print bnfo1.reverse_complement(seq)

def problem1_2_10():
    text, k = lines_to_params('Data/bnfo1_week1/dataset_2_10.txt')
    print_list(bnfo1.most_frequent_kmers(text,int(k)))

def problem1_2_7():
    text, pattern = lines_to_params('Data/bnfo1_week1/dataset_2_7.txt')
    print bnfo1.pattern_count(text,pattern) 

''' Utilities '''

def print_list(listobj, delim=' '):
    ''' Prints elements in a list separated by delimitor '''
    print delim.join(map(str,listobj))

def lines_to_params(filename):
    ''' Creates an array of strings corresponding to lines in the file '''
    lines = []
    for line in file(filename,'r+'):
        lines.append(line.strip())
    return lines

main()
#seq,params = lines_to_params('Data/dataset_4_5.txt')
#seq = lines_to_params('Data/E_coli.txt')[0]
#k = 9; L = 500; t = 3

#soln = bnfo1.find_clumped_kmers(seq,k,L,t)
#print_list(soln)
#print len(soln)
