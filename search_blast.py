#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import argparse
from icecream import ic
import sys
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
from io import StringIO
from typing import List
from search_rules import SnakemakeSearchRules
import config

from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr), includeContext=True)
ic.configureOutput(prefix='')

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--db', help="[pdbnt, refseq, nt]", default="pdb")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="a fast file with one sequence", default="") # nargs='+')
    return parser

def process_blast_output(blast_output_file):
    # Open the BLAST output XML file
    with open(blast_output_file) as result_handle:
        # Parse the BLAST output
        blast_records = NCBIXML.parse(result_handle)

        # Iterate over each BLAST record and get the relevant information
        for blast_record in blast_records:
            
            query_id = blast_record.query_id
            query_length = blast_record.query_length
            # Iterate over each alignment and get the relevant information
            if not (blast_record.alignments):
                print('No hits found')
                return
            
            for alignment in blast_record.alignments:
                subject_id = alignment.title
                ic(subject_id)
                subject_length = alignment.length

                # Iterate over each HSP and get the relevant information
                for hsp in alignment.hsps:
                    alignment_length = hsp.align_length
                    query_alignment = hsp.query
                    subject_alignment = hsp.sbjct
                    # Process the alignment information as needed

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    import os
    search_rules:SnakemakeSearchRules = SnakemakeSearchRules()
    run_blast_search:bool = search_rules.do_blast_run()
    
    blast_result_string:str = ''
    if run_blast_search is True:        
        for f in args.file:
            ## from Bio import SeqIO
            ## seq_record = next(SeqIO.parse(open(f),'fasta'))
            ## from Bio.Blast.Applications import NcbiblastnCommandline
            ## blastn_cline = NcbiblastnCommandline(query = "search.fasta", db = "nr", 
            ##                                      outfmt = 5, out = "results.xml")         
            
            ## stdout, stderr = blastn_cline()
            fh = open(f)
            h = fh.readline().strip()
            seq = fh.readline().strip()
            print(h)
            print(seq)
            if args.db == 'nt':
                db = config.BLAST_NT
            if args.db == 'pdbnt':
                db = '/home/rnahub/rnahub/db/pdbnt/pdbnt'
            if args.db == 'refseq':
                db = '/home/rnahub/rnahub/db/refseq_rna/refseq_rna'
            #db = '/home/rnahub/mnts/harvard/RoseTTAFold2NA/RNA/nt'
            #  -outfmt=5  # xml
            #  -outfmt=6  # regular file
            blast_output_file = 'example/output_file.xml'
            try:
                os.remove(blast_output_file)
            except FileNotFoundError:
                pass
            cmd = 'time blastn -num_threads 6 -db ' + db + ' -query ' + f + ' -outfmt=5 -out ' + blast_output_file

            """
            Need to change the bellow to make the output be sent to the stdout
            vs a file in order to make work correctly with snakemake - JMP
            """
            
            print(cmd)
            os.system(cmd)
            
            #for now send blast file to stout so that snakemake will have a rule to process the output
            #as it has its own logic and without a output it does not think this should be ran
            blast_results_text:List[str] = []
            with open(blast_output_file, 'r') as file:
                blast_results_text = file.readlines()
                
            blast_result_string = ''.join(blast_results_text)        
            process_blast_output(blast_output_file)
            ic(blast_output_file)
    else:
        # setting result string to call out rfam search as not sure what other search we are doing and wanted to show an example of use¯\_(ツ)_/¯
        blast_result_string = f'BLAST not ran\nFound good results in rfam search.\n'
    sys.stdout.write(blast_result_string) 
