#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
https://biopython.org/docs/1.75/api/Bio.SearchIO.HmmerIO.html
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


from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr), includeContext=True)
ic.configureOutput(prefix='')
import subprocess

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--db', help="[pdbnt, refseq, nt]", default="pdbnt")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="a fast file with one sequence", default="") # nargs='+')
    return parser

import subprocess
def exe(cmd):
    o = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out = o.stdout.read().strip().decode()
    err = o.stderr.read().strip().decode()
    return out, err


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    input_file = args.file[0]
    output_file = "nhmmer_output.txt"

    if args.db == 'nt':
        db = '/home/rnahub/mnt/nt/nt'
    if args.db == 'pdbnt':
        db = '/home/rnahub/rnahub/db/pdbnt/pdbnt.fsa'
    if args.db == 'refseq':
        db = '/home/rnahub/rnahub/db/refseq_rna/refseq_rna'

    # Construct the nhmmer command
    command = f'/usr/local/bin/nhmmer "--tblout" {output_file} -E 0.005 --cpu 5 {input_file} {db} ' # -E 0.001 

    # Run the command
    ic(command)
    stdout, stderr = exe(command)

    # Check for errors
    if stderr:
        print("Error running nhmmer:")
        print(stderr)
    else:
        print("nhmmer ran successfully. Output saved in", output_file)

    from Bio import SearchIO

    # Replace 'hmmer_output_file.txt' with the path to your HMMER output file
    hmmer_output_file = 'nhmmer_output.txt'

    # Parse the HMMER output file
    query_results = SearchIO.parse(hmmer_output_file, 'hmmer3-text')

    # Iterate over each QueryResult
    ic([x for x in query_results]) 
    for query in query_results:
        print(f"Query: {query.id}")
        for hit in query.hits:
            print(f"  Hit: {hit.id}, E-value: {hit.evalue}")

            # Iterate over each HSP (high-scoring pair)
            for hsp in hit.hsps:
                print(f"    HSP, E-value: {hsp.evalue}")
        
