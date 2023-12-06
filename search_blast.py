#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
from io import StringIO

from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr), includeContext=True)
ic.configureOutput(prefix='')

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--blast', help="[pdbnt, refseq, nt]", default="pdb")
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
            ic(blast_record.query)
            # Iterate over each alignment and get the relevant information

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

    os.system('figlet -f smblock rnahub')
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
        if args.blast == 'nt':
            pass
        if args.blast == 'pdbnt':
            # db = '/home/rnahub/mnt/pdb/pdbnt'
            db = '/home/rnahub/rnahub/db/pdbnt/pdbnt'
        if args.blast == 'refseq':
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

        print(cmd)
        os.system(cmd)
        process_blast_output(blast_output_file)
        ic(blast_output_file)
