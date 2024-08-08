#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 15:24:57 2022

@author: ayang
Script for extracting noncoding coordinates of .sto alignment after 1 
round of flanked nhmmer
"""
import os
import re
#%%
def parser(infile):
    '''
    Parse the stockholm file made after 1 round of nhmmer
    using the flanked sequence

    Parameters
    ----------
    infile : string
        path to stockholm file made from 1 round of nhmmer
        on the flanked sequence

    Returns
    -------
    seq : string
        reference sequence for fugu in the .sto alignment
    start : integer
        coordinate of start of seq
    end : integer
        coordinate of end of seq

    '''
    query = "" # emtpy strings to append to
    identifier = "" # emtpy strings to append to
    seq = ""
    with open(infile, 'r') as fh: # open the .sto file for reading
        for line in fh: # go through every line of the file
            if "ID" in line: # take the reference sequence (first instance of fugu)
                this = line.split() # split on white space
                query += this[2] # take the identifier with coordinates
                break # break to avoid taking second, third, etc. instances of fugu
        
        for line in fh: # go through every line of the file
            if "#=GS NC_" in line: # take the reference sequence (first instance of fugu)
                this = line.split() # split on white space
                identifier += this[1] # take the identifier with coordinates
                break # break to avoid taking second, third, etc. instances of fugu

        #print("\n^^query", query);
        #print("\n^^IDENT", identifier);
        
        for line in fh: # go through the file
            if line.startswith(identifier): # and search for the identifier in the alignment portion
                seq += line.split()[1] # get the sequence from the alignment
        #print("\n^^seq", seq);
    
    # get the exonU and exonD start/send from the query
    m = re.search(r'\.(\d+)-(\d+)/\S+/(\d+)-(\d+)$', query)
    exonU_start = int(m.group(1)) 
    exonU_end = int(m.group(2))
    exonD_start = int(m.group(3)) 
    exonD_end = int(m.group(4))
    #print("exonU start/end", exonU_start, exonU_end,"exonD start/end", exonD_start, exonD_end)
    
    # get the start and end from the identifier
    m = re.search(r'(\d+)-(\d+)', identifier)
    start = int(m.group(1)) 
    end = int(m.group(2))
    #print("start/end", start, end)
    return seq, start, end, query, exonU_start, exonU_end, exonD_start, exonD_end    


def make_counters(start_coord, sequence):
    '''
    Matches a position in the sequence to a column position 
    in the .sto alignment

    Parameters
    ----------
    start_coord : integer
        start coordinate of sequence in the genome
    sequence : string
        sequence from the .sto alignment, with gaps

    Returns
    -------
    seq_counter : list of integers
        counts each position in the sequence with "@" placeholders
    col_counter : list of integers
        counts each column of the .sto alignment
    * seq_counter and col_counter should be paired so that the elements at each
    index correspond to each other

    '''
    seq_counter = [] # empty lists to append to
    col_counter = []
    for i in range(len(sequence)): # go through the sequence
        col_counter.append(i+1) # start counting the columns
        if sequence[i].isalpha(): # start counting along the sequence
            seq_counter.append(start_coord) # use genomic coordinate
            start_coord += 1 # increment
        else:
            seq_counter.append("@") # use placeholders so indices match
    return seq_counter, col_counter


cwd = os.getcwd() # current working directory where this script is saved
for subdirs, dirs, files in os.walk(cwd): # go through all files in the cwd
    for file in files:
        if file == "v0.sto": # only consider the files named v0.sto (the file made after first round of flanked nhmmer)
            complete_file = os.path.join(subdirs, file) # get the complete path
            outfile = os.path.join(subdirs, 'column_nums.txt') # make a file that includes the sequence coordinates' positions as columns
            f = open(outfile, 'w') # open outfile for writing

            alignment_seq, alignment_start, alignment_end, query, exonU_start, exonU_end, exonD_start, exonD_end = parser(complete_file) # apply the parser
            alen = len(alignment_seq);
            exonU_len = exonU_end - exonU_start + 1;
            exonD_len = exonD_end - exonD_start + 1;
            print(alen, "alignment ", alignment_start, alignment_end, "exonU", exonU_start, exonU_end, exonU_len, "exonD", exonD_start, exonD_end, exonD_len);
            
            bp_counter, column_counter = make_counters(alignment_start, alignment_seq) # apply the counter
            print(bp_counter)
            
            col_num_start = 1
            if (exonU_end > 0):
                for i in range(alen):
                    if bp_counter[i] == exonU_end:
                        col_num_start = i+1 + 1
                        break
                    
            col_num_end = alen
            if (exonD_start > 0):
                for i in range(alen):
                    if bp_counter[i] == exonD_start:
                        col_num_end = i+1 - 1
                        break
            print (" start", col_num_start, "end", col_num_end);
                
            f.write("{0:50}{1:20}{2:20}".format(query, col_num_start, col_num_end)) # write column info to a file
            f.close()
