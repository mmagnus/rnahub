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
from config import RSCAPE_PATH, nhmmer, SCRIPTS_DIR, CPUs, EASEL_PATH
import logging

# SLURM directives are not directly used in Python scripts.
# Instead, configure your job submission script or environment accordingly.

def now():
    import datetime
    print(datetime.datetime.now())

def exe(command, dry=False, logging=False):
    """Execute a shell command."""
    print(f'cmd: {command}', flush=True)
    if logging:
        now()
        logger.info(command)
        logger.info(command)
    if not dry:
        os.system(command)
    return
    # try:
    #     subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # except subprocess.CalledProcessError as e:
    #     print(f"Error executing command: {e.cmd}")
    #     print(e.output.decode())
    #     sys.exit(1)

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('db', help="", default="")
    parser.add_argument('--job-name', default="", help="by default is input file name (wihout extension)")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--slurm",  action="store_true", help="send it to slumrm")
    parser.add_argument("-f", "--flanked",  action="store_true", help="run flanked mode (create extra v0)")
    parser.add_argument("--evalue", default="1e-5", help="e-value threshold")
    parser.add_argument("--iteractions", default=3, help="number of iterations", type=int)
    parser.add_argument("--dry", help="show all cmds, dont run them", action="store_true")
    parser.add_argument("--dev-skip-nhmmer0", help="show all cmds, dont run them", action="store_true")
    parser.add_argument("--rscape", help="rscape only",
                        action="store_true")
    parser.add_argument("file", help=".fa for now, don't use .fasta", default="", nargs='+')
    
    return parser

def clean():
        for pattern in ['v0*', 'v1*', 'v2*', 'v3*', 'rm_v3.sto']:#, './*.txt', './*/']:
            exe(f'rm -f {pattern}', dry)

def search():
    """Return the last sto file generated"""
    if args.flanked:
        def bp_col(alignment, seq_spec_string):
            ic(seq_spec_string)
            # directory of alignment file that you provide, the output will go into there
            DIR = os.path.dirname(alignment)

            # Read the alignment file to get the top sequence that contains reference genome
            with open(alignment, 'r') as f:
                for line in f:
                    if seq_spec_string in line:
                        seq_name = line.split()[1]
                        break
            ic(seq_name)
            if not seq_name:
                print(f"No sequence found matching {seq_spec_string}")
                sys.exit(1)

            # Extract only the sequence (gaps included) of that sequence from the MSA
            sequence = []
            with open(alignment, 'r') as f:
                for line in f:
                    if seq_name in line and not '/' in line.split()[1]:
                        sequence.append(line.split()[1])

            ic(sequence)           
            f = ''.join(sequence)
            ic(f, len(f))

            # Now create a counter that loops through and increments as it reads each character.
            bp = []
            array = []
            counter = 0
            for c in f:
                    counter += 1
                    if c in ['U', 'A', 'T', 'G', 'C']:
                        if c == 'U':
                            c = 'T'
                        bp.append(c)
                        array.append(counter)

            # Write header and results to bp_col.txt
            bp_col_txt_path = os.path.join(DIR, 'bp_col.txt')
            with open(bp_col_txt_path, 'w') as f:
                f.write("bp \t Col\n")
                for i, each in enumerate(array, start=1):
                    f.write(f"{i} \t {each}\n")

        def extractor(query):
            # gkab355_supplemental_files/tutorial/method_scripts/tutorial/YAR014C_plus_IGR/
            # Extract the relevant information from the query file
            ic(query)
            with open(query, 'r') as fn:
                s1 = None
                s2 = None
                for line in fn:
                    if line.startswith('>'):
                        print(line)
                        parts = line.split()
                        s1 = parts[3]
                        s2 = parts[4]
                        break

            if s1 is None or s2 is None:
                print("Failed to extract s1 or s2 from the query file.")
                sys.exit(1)

            ic(s1, s2)

            # Locate the first and second site in the bp_col.txt file
            bp_col_path = os.path.join(f'{j}', 'bp_col.txt')
            ic(bp_col_path)

            first_site = ''
            second_site = ''

            with open(bp_col_path, 'r') as fl:
                for line in fl:
                    parts = line.split()
                    print(parts, s1, s2)
                    if parts[0] == s1:
                        first_site = parts[1]
                    if parts[0] == s2:
                        second_site = parts[1]

            ic(s1, first_site, second_site)
            #if first_site is None or second_site is None:
            #        print("Failed to locate first_site or second_site in the bp_col.txt file.")
            #        sys.exit(1)
            cmd = ''.join([f'{EASEL_PATH}/esl-alimask -t ', j, '/flanked.sto ', first_site, '..', second_site, ' > ', j, '/noncoding.sto'])
            print(cmd)
            exe(cmd)

            cmd = ''.join([f'{EASEL_PATH}/esl-alimanip --lmin 50  ', j, '/noncoding.sto > ', j, '/trim_noncoding.sto'])
            print(cmd)
            exe(cmd)
            # esl-alimanip --lmin 50 tutorial/YAR014C_plus_IGR/noncoding.sto

        # v0 is flanked
        #cmd = nhmmer -E 1e-10 -A $DIR/$filename/flanked.sto --tblout $DIR/$filename/flanked.hmmout $query $DB > out.txt
        # nhmmer -E 1e-10 --cpu 64 -A tutorial/gly1_igr/flanked.sto --tblout tutorial/gly1_igr/flanked.hmmout tutorial/gly1_igr.fa ../../../db/1409_Acomycota_genomes-may19.fa
        cmd = f"{nhmmer} --cpu {CPUs} --incE 1e-10 -A {j}/flanked.sto {j}/{fbase}.fa {db} > {j}/flanked.out" #v0.sto is flanked # fa vs fasta #TODO
        if not args.dev_skip_nhmmer0:
            exe(cmd, dry = False)
        # subscripts/bp_col.sh tutorial/YAR014C_plus_IGR/flanked.sto S288C
        cmd = f'{SCRIPTS_DIR}/bp_col.py {j}/flanked.sto S288C' #
        dry = False
        #exe(cmd, dry)
        bp_col(f'{j}/flanked.sto', 'S288C')

        #+ subscripts/noncoding_extractor.sh tutorial/YAR014C_plus_IGR.fasta esl-alimask
        #+ esl-alimask -t tutorial/YAR014C_plus_IGR/flanked.sto 2305..
        extractor(f'{j}/{fbase}.fa') # noncoding.sto and trim_noncoding.sto

        query = f'{j}/trim_noncoding.sto' # query is now the noncoding.sto, not a single sequence
    for i in range(1, nofinteractions + 1):  # you can play with this one, starting from 1
            sto_file = f'v{i}.sto'
            output_file = f'v{i}.out'
            input_file = query if i == 1 else f'{j}/v{i-1}.sto'
            command = f"{nhmmer} --cpu 2 --incE {evalue} -A {j}/{sto_file} {input_file} {db} > {j}/{output_file}"
            exe(command, dry)
def remove_multicopies():
        scriptsdir = './'
        exe(f'python {scriptsdir}/remove_multicopies5.py', dry)
        
def esl_aliminip():
        # Run esl-alimanip
        eslalimanip = os.path.join(rscapedir, "lib/hmmer/easel/miniapps/esl-alimanip")
        exe(f"{eslalimanip} --seq-k accessions_to_keep.txt v3.sto > rm_v3.sto", dry)

def rscape():
    # Set up for R-scape analysis
    #exe('rm -f rscape_results.txt')
    #exe('rm -rf rscape_output')
    try:
        os.makedirs(f'{j}/rscape_output')
    except FileExistsError:
       pass
    #exe(f"{RSCAPE_PATH} --outdir {job_folder}/rscape_output --cacofold --outtree rm_v3.sto > rscape_results.txt")
    exe(f"{RSCAPE_PATH} --outdir {j}/rscape_output --cacofold --outtree {j}/v{nofinteractions}.sto | tee {j}/rscape_results.txt", dry)

def save_to_slurm():
        name = f'{dbbase}X{fbase}'
        t = f"""#!/bin/bash

#SBATCH -n 3 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0 # Runtime in minutes
#SBATCH -p eddy # Partition to submit to
#SBATCH --mem=50000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o {j}/%j.out # Standard out goes to this file
#SBATCH -e {j}/%j.err # Standard err goes to this filehostname
#SBATCH --job-name={j} # Name of the job
# Runs a command on all FASTA files in current directory
python {' '.join(sys.argv[:]).replace('--slurm', '')}
        """

        with open(f'{j}/run.slurm', 'w') as fi:
            fi.write(t)
        print(f'sbatch {j}/run.slurm')
        os.chmod(f'{j}/run.slurm', 0o755) #mode=stat.S_IXUSR)
        os.system(f'sbatch {j}/run.slurm')
        
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    db = args.db
    evalue = args.evalue
    nofinteractions = args.iteractions
    dry = args.dry
    
    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:
        query = f
        fbase = os.path.basename(f).replace('.fa', '')
        dbbase = os.path.basename(db).replace('.fa', '')
        if args.job_name:
            j = 'jobs/' + args.job_name
        else:
            j = 'jobs/' + fbase
        try:
            os.makedirs(f'{j}', exist_ok=True) 
        except FileExistsError:
            pass
        shutil.copy(f, j)

        now()

        print(args, flush=True)

        with open(f'{j}/cmd.sh', 'w') as fi:
            fi.write(' '.join(sys.argv))
        
        if args.slurm:
            save_to_slurm()
            exit()

        import logging
        # ic(f'{j}/log.log')
        # logging.basicConfig(filename=f'{j}/log.log', filemode='w',
        #                     encoding='utf-8', level=logging.INFO,
        #                     format='%(asctime)s - %(levelname)s - %(name)s: %(message)s')
        # logging.info(str(args))
        # logger = logging.getLogger(__name__)
        # with open(f'{j}/log2.log', 'w') as f:
        #     f.write(str(args))
        # create logger with 'spam_application'
        logger = logging.getLogger('rnahub')
        logger.setLevel(logging.DEBUG)
        # create file handler which logs even debug messages
        fh = logging.FileHandler(f'{j}/log.log')
        fh.setLevel(logging.DEBUG)

        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)  # You can adjust this as needed

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        
        logger.addHandler(fh)
        logger.addHandler(ch)
        
        logger.info('start')
        logger.info(str(args))
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

        logging.info('done')
        logger.info('done')
        print('done', flush=True)
        now()
