#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 24 09:58:47 2021

@author: ayang

Script for taking top scoring hit per genome from last iteration of nhmmer
Determines which sequences to keep and writes accessions to a file (accessions_to_keep.txt)
Can then run eslalimanip.sh to keep only these sequences
"""
import os
directory = os.getcwd()

#%% 
# parses last alignment built by nhmmer, v3.sto   
def parse_last_iteration(alignment):
    '''
    Parameters
    ----------
    alignment : string
        name of a stockholm alignment made by nhmmer

    Returns
    -------
    genome : string
        the genome from which each sequence in the alignment is from (as species name)
    accession : string
        full ensembl accession for each sequence, includes coordinates
    '''
    
    genome = [] # name of species in which the hit was found
    accession = [] # full accession number, includes coordinates
    
    with open(alignment, "r") as fh: # open the file (alignment) for reading
        for line in fh: # go through every line in the alignment
            if line.startswith("#=GS"): # only consider lines of the alignment portion
                this = line.split()
                
                # parse using indexes
                genus = this[5]
                species = this[6]
                
                for character in genus:
                    if character.isalpha()==False:
                        genus = genus.replace(character, '')
                for character in species:
                    if character.isalpha()==False:
                        species = species.replace(character, '')
                
                this_genome = genus + species
                this_accession = this[1]
                
                # add to list
                genome.append(this_genome)
                accession.append(this_accession)
        
    return genome, accession

#%%
# only keep the top scoring hit from each genome
encountered_genomes = [] # for keeping track of first instance of a genome
accessions_to_keep = [] # and corresponding accession

for root, dirs, files in os.walk(directory): # go through every file in the directory containing .sto files made by nhmmer
    for file in files:
        if file == "accessions_to_keep.txt": # writes over old files, if there are any
            os.remove(file)
            
        elif file == "v3.sto": # go through the .sto made by the last iteration of nhmmer
            genomes, accessions = parse_last_iteration(file) # parse
            num_seqs = len(genomes)
            
            for j in range(num_seqs): # go through all the genomes in the alignment
                this_genome = genomes[j] 
                this_accession = accessions[j]
                
                if this_genome not in encountered_genomes: # save only the first instance
                    encountered_genomes.append(this_genome)
                    accessions_to_keep.append(this_accession)
        else:
            pass
                
#%%
# write the accessions to keep to a file
num_accessions_to_keep = len(accessions_to_keep)
filehandle = "accessions_to_keep.txt"
f = open(filehandle, "w")
for j in range(num_accessions_to_keep):
    f.write("%s\n" % (accessions_to_keep[j]))
f.close()

