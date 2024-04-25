#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')
import shutil
import subprocess
import sys
import os
from config import *
# SLURM directives are not directly used in Python scripts.
# Instead, configure your job submission script or environment accordingly.


RSCAPE_PATH = '/n/home06/mmagnus//m/opt/rscape_v2.0.4.b/bin/R-scape'

def exe(command, dry=False):
    """Execute a shell command."""
    print(command)
    if not dry: os.system(command)
    return
    try:
        subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e.cmd}")
        print(e.output.decode())
        sys.exit(1)

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('db', help="", default="")
    parser.add_argument('--job-name', default="", help="by default is input file name (wihout extension)")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--evalue", default="1e-5", help="e-value threshold")
    parser.add_argument("--rscape", help="rscape only",
                        action="store_true")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


def clean():
        for pattern in ['v1*', 'v2*', 'v3*', 'rm_v3.sto']:#, './*.txt', './*/']:
            exe(f'rm -f {pattern}')

def search():
    for i in range(1, 4):  # you can play with this one, starting from 1
            sto_file = f'v{i}.sto'
            output_file = f'v{i}.out'
            input_file = query if i == 1 else f'{j}/v{i-1}.sto'
            command = f"{nhmmer} --cpu 2 --incE {evalue} -A {j}/{sto_file} {input_file} {db} > {j}/{output_file}"
            exe(command)

def remove_multicopies():
        scriptsdir = './'
        exe(f'python {scriptsdir}/remove_multicopies5.py')
        
def esl_aliminip():
        # Run esl-alimanip
        eslalimanip = os.path.join(rscapedir, "lib/hmmer/easel/miniapps/esl-alimanip")
        exe(f"{eslalimanip} --seq-k accessions_to_keep.txt v3.sto > rm_v3.sto")

def rscape():
    # Set up for R-scape analysis
    #exe('rm -f rscape_results.txt')
    #exe('rm -rf rscape_output')
    try:
        os.makedirs(f'{j}/rscape_output')
    except FileExistsError:
       pass
    #exe(f"{RSCAPE_PATH} --outdir {job_folder}/rscape_output --cacofold --outtree rm_v3.sto > rscape_results.txt")
    exe(f"{RSCAPE_PATH} --outdir {j}/rscape_output --cacofold --outtree {j}/v3.sto > {j}/rscape_results.txt")

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    db = args.db
    evalue = args.evalue
    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        query = f
        if args.job_name:
            j = 'jobs/' + args.job_name
        else:
            j = 'jobs/' + os.path.basename(f).replace('.fa', '')

        try:
            os.mkdir(f'{j}')
        except FileExistsError:
            pass
        shutil.copy(f, j)

        #scriptsdir = "/n/eddy_lab/users/erivas/projects/SKennedy/2024_conserved_introns/shscripts/unflanked_scripts"
        # Clean up previous output files
        #clean()
        # Perform nhmmer iterations
        if not args.rscape:
            search()
        # Remove duplicate copies of genomes
        #remove_multicopies()
        # Module load (if necessary in your environment; may require a system-specific approach)
        #exe('module load Anaconda2/2019.10-fasrc01')
        #esl_aliminip()
        rscape()
