#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
job_path directory (j before)

.. warning ::

   R-scape version, searching for .helixcov files here!

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
from config import RSCAPE_PATH, nhmmer, EASEL_PATH, RFAM_DB_PATH, REPEAT_MASKER_PATH
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

    parser.add_argument('--db', help="", default="", nargs='+')
    parser.add_argument('--job-name', default="", help="by default is input file name (wihout extension)")
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("--slurm",  action="store_true", help="send it to slumrm")
    parser.add_argument("--evalue", default="1e-10", help="e-value threshold for all the runs but the final one")
    parser.add_argument("--evalue-final", default="1e-5", help="e-value threshold for the final run")
    parser.add_argument("--iteractions", default=3, help="number of iterations", type=int)
    parser.add_argument("--cpus", default=2, help="number of cpus for nhmmer", type=int)
    parser.add_argument("--dry", help="show all cmds, dont run them", action="store_true")
    parser.add_argument("--repeatmasker", help="", action="store_true")
    parser.add_argument("--create-job-folder", help="create a job folder based on the path to the input fasta sequence, e.g. example/seq.fa, jobs/seq/seq.fa", action="store_true")
    parser.add_argument("--dev-skip-search", help="skip v0..v3 all nhmmer searches", action="store_true")
    parser.add_argument("--dev-skip-nhmmer0", help="show all cmds, dont run them", action="store_true")
    parser.add_argument("--dev-skip-nhmmer123", help="show all cmds, dont run them", action="store_true")
    parser.add_argument("--dev-skip-cmcalibrate", help="show all cmds, dont run them", action="store_true")
    parser.add_argument("--dev-skip-rscape", help="show all cmds, dont run them", action="store_true")
    parser.add_argument("--dev-skip-infernal", help="show all cmds, dont run them", action="store_true")

    parser.add_argument("--rscape", help="rscape only",
                        action="store_true")
    parser.add_argument("fasta", help=".fa for now, don't use .fasta")#, default="", nargs='+')
    parser.add_argument("--flanked", help=".fa for now, don't use .fasta, flank the sequence including the query sequence")
    parser.add_argument("--flanks-in-header", action="store_true", help="run flanked mode (create extra v0 files), syntax in the fasta header '><seq_name> <start>-<end>, use this or --flanked fasta file")
    
    return parser

def clean():
        for pattern in ['v0*', 'v1*', 'v2*', 'v3*', 'rm_v3.sto']:#, './*.txt', './*/']:
            exe(f'rm -f {pattern}', dry)

def search(seq_path, seq_flanked_path = ''):
    """Return the last sto file generated"""
    query = seq_path
    print(f'query: {seq_path}')
    
    if args.flanked or args.flanks_in_header:
        def bp_col(alignment):
            # directory of alignment file that you provide, the output will go into there
            DIR = os.path.dirname(alignment)
            # Read the alignment file to get the top sequence that contains reference genome
            with open(alignment, 'r') as f:
                for line in f:
                    if '#=GS' in line:
                        seq_name = line.split()[1]
                        ic(seq_name)
                        break

            if not seq_name:
                print(f"No sequence found matching")
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

        def extractor(query, verbose=False):
            # gkab355_supplemental_files/tutorial/method_scripts/tutorial/YAR014C_plus_IGR/
            # Extract the relevant information from the query file
            # get full file and find substring from the long one
            if not args.flanks_in_header:
                ic(query)
                with open(query, 'r') as fn:
                    lines = fn.readlines()
                    sequence_only = "".join(line.strip().upper() for line in lines if not line.startswith('>'))
                ic(sequence_only)


                with open(args.flanked, 'r') as fn:
                    lines = fn.readlines()
                    flanked_sequence_only = "".join(line.strip().upper() for line in lines if not line.startswith('>'))
                ic(flanked_sequence_only)

                # Find the beginning index
                start_index = flanked_sequence_only.find(sequence_only)

                # Calculate the end index if substring is found
                if start_index != -1:
                    end_index = start_index + len(sequence_only) - 1
                    print(f"Substring found from index {start_index + 1} to {end_index}.")
                else:
                    raise Exception("Substring not found.")

                s1 = str(start_index + 1)
                s2 = str(end_index)

            # parse the header of the fasta file
            else:
                with open(query, 'r') as fn:
                    s1 = None
                    s2 = None
                    for line in fn:
                        if line.startswith('>'):
                            print(line)
                            parts = line.split()
                            ic(parts)
                            #>target 5-20
                            s1, s2 = parts[1].split('-') ##3]
                            ic(s1, s2)
                            break

            if s1 is None or s2 is None:
                print("Failed to extract s1 or s2 from the query file.")
                sys.exit(1)

            ic(s1, s2)
            
            # Locate the first and second site in the bp_col.txt file
            bp_col_path = os.path.join(f'{job_path}', 'bp_col.txt')

            first_site = ''
            second_site = ''

            with open(bp_col_path, 'r') as fl:
                for line in fl:
                    parts = line.split() #> parts: ['16', '16']
                    if parts[0] == s1: 
                        first_site = parts[1]
                    if parts[0] == s2:
                        second_site = parts[1]

            ic(s1, first_site, second_site) # take the last number if you can't find s2 ;-)

            #if first_site is None or second_site is None:
            #        print("Failed to locate first_site or second_site in the bp_col.txt file.")
            #        sys.exit(1)
            cmd = ''.join([f'{EASEL_PATH}/esl-alimask -t ', job_path, '/v0.sto ', first_site, '..', second_site, ' > ', job_path, '/v0_targetRegionOnly.sto'])
            print(cmd)
            exe(cmd)

            cmd = ''.join([f'{EASEL_PATH}/esl-alimanip --lmin 50  ', job_path, '/v0_targetRegionOnly.sto > ', job_path, '/v0_targetRegionOnly_trim.sto'])
            print(cmd)
            exe(cmd)
            # esl-alimanip --lmin 50 tutorial/YAR014C_plus_IGR/noncoding.sto

        # v0 is flanked
        #cmd = nhmmer -E 1e-10 -A $DIR/$filename/flanked.sto --tblout $DIR/$filename/flanked.hmmout $query $DB > out.txt
        # nhmmer -E 1e-10 --cpu 64 -A tutorial/gly1_igr/flanked.sto --tblout tutorial/gly1_igr/flanked.hmmout tutorial/gly1_igr.fa ../../../db/1409_Acomycota_genomes-may19.fa#

        # {j}/{fbase}.fa is causing missing organism problem
        if args.flanks_in_header:
            cmd = f"cat {db} | {nhmmer}  --noali --cpu {CPUs} --incE {args.evalue} -A {job_path}/v0.sto {job_path}/{fbase}.fa - "#| tee {j}/v0.out"
        else:
            cmd = f"cat {db} | {nhmmer}  --noali --cpu {CPUs} --incE {args.evalue} -A {job_path}/v0.sto {seq_flanked_masked_path} - "#| tee {j}/v0.out"            
        #v0.sto is flanked.sto # fa vs fasta #TODO
        if not args.dev_skip_nhmmer0:
            exe(cmd, dry)
        # subscripts/bp_col.sh tutorial/YAR014C_plus_IGR/flanked.sto S288C
        #cmd = f'{SCRIPTS_DIR}/bp_col.py {j}/flanked.sto S288C' #
        #dry = False
        #exe(cmd, dry)
        bp_col(f'{job_path}/v0.sto') # !!!!!!!!!!!!!!!!!!!!!!!!! #TODO
        #+ subscripts/noncoding_extractor.sh tutorial/YAR014C_plus_IGR.fasta esl-alimask
        #+ esl-alimask -t tutorial/YAR014C_plus_IGR/flanked.sto 2305..
        extractor(f'{job_path}/{fbase}.fa') # noncoding.sto and trim_noncoding.sto
        query = f'{job_path}/v0_targetRegionOnly_trim.sto'

    if not args.dev_skip_nhmmer123:
        for i in range(1, nofinteractions + 1):  # you can play with this one, starting from 1
            sto_file = f'v{i}.sto'
            output_file = f'v{i}.out'
            input_file = query if i == 1 else f'{job_path}/v{i-1}.sto'
            evalue = args.evalue
            if i == args.iteractions:
                evalue =  args.evalue_final
            #  {j}/{fbase}.fa
            command = f"time cat {db} | {nhmmer} --noali --cpu {CPUs} --incE {evalue} -A {job_path}/{sto_file} {input_file} - "#| tee {j}/{output_file}"
            exe(command, dry)

    # statistics for v3
    cmd = ''.join([f'{EASEL_PATH}/esl-alistat ', job_path, '/v3_rm.sto > ', job_path, '/v3_rm_stats.txt'])
    print(cmd)
    exe(cmd, dry)
    
            
def find_top_scoring_hits(directory=None, output_file="accessions_to_keep.txt"):
    """
    Created on Mon May 24 09:58:47 2021

    @author: ayang

    Script for taking top scoring hit per genome from last iteration of nhmmer
    Determines which sequences to keep and writes accessions to a file (accessions_to_keep.txt)
    Can then run eslalimanip.sh to keep only these sequences

    Finds the top scoring hit per genome from the last iteration of nhmmer and writes 
    the accessions to a file.
    
    Parameters
    ----------
    directory : str, optional
        The directory to search for .sto files. Defaults to the current working directory.
    output_file : str, optional
        The output file where accessions to keep will be written. Defaults to 'accessions_to_keep.txt'.
    """
    def parse_last_iteration(alignment):
        '''
        Parses the last alignment built by nhmmer.

        Parameters
        ----------
        alignment : string
            Name of a stockholm alignment made by nhmmer.

        Returns
        -------
        genome : list of strings
            The genomes from which each sequence in the alignment is from (as species name).
        accession : list of strings
            Full ensembl accession for each sequence, includes coordinates.
        '''
        
        genome = [] # name of species in which the hit was found
        accession = [] # full accession number, includes coordinates
        print(alignment)
        with open(alignment, "r") as fh: # open the file (alignment) for reading
            for line in fh: # go through every line in the alignment
                if line.startswith("#=GS"): # only consider lines of the alignment portion
                    this = line.split()
                    
                    # parse using indexes
                    genus = this[5]
                    species = this[6]
                    
                    for character in genus:
                        if character.isalpha() == False:
                            genus = genus.replace(character, '')
                    for character in species:
                        if character.isalpha() == False:
                            species = species.replace(character, '')
                    
                    this_genome = genus + species
                    this_accession = this[1]
                    
                    # add to list2
                    genome.append(this_genome)
                    accession.append(this_accession)
            
        return genome, accession

    # Initialize lists to keep track of genomes and accessions
    encountered_genomes = [] # for keeping track of first instance of a genome
    accessions_to_keep = [] # and corresponding accession

    # Go through every file in the directory containing .sto files made by nhmmer
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file == output_file: # writes over old files, if there are any
                os.remove(os.path.join(root, file))
                
            elif file.endswith(".sto"): # go through the .sto files
                #v0.sto", "v1.sto", "v2.sto", 
                if file in [f"v{nofinteractions}.sto"]: # only consider the last iteration
                    genomes, accessions = parse_last_iteration(os.path.join(root, file)) # parse
                    num_seqs = len(genomes)

                    for j in range(num_seqs): # go through all the genomes in the alignment
                        this_genome = genomes[j] 
                        this_accession = accessions[j]

                        if this_genome not in encountered_genomes: # save only the first instance
                            encountered_genomes.append(this_genome)
                            accessions_to_keep.append(this_accession)
                    print(this_genome, this_accession)
            
    # Write the accessions to keep to a file
    with open(os.path.join(directory, output_file), "w") as f:
        for accession in accessions_to_keep:
            f.write(f"{accession}\n")
    
    print(f"Accessions to keep have been written to {output_file}.")
    cmd = f'{EASEL_PATH}/esl-alimanip --seq-k {directory}/accessions_to_keep.txt {directory}/v{nofinteractions}.sto > {directory}/v{nofinteractions}_rm.sto'
    print(cmd)
    exe(cmd)

def rscape():
    # Set up for R-scape analysis
    #exe('rm -f rscape_results.txt')
    #exe('rm -rf rscape_output')
    try:
        os.makedirs(f'{job_path}/rscape_output')
    except FileExistsError:
       pass
    exe(f"{RSCAPE_PATH} --outdir {job_path}/rscape_output --cacofold --outtree {job_path}/v{nofinteractions}_rm.sto | tee {job_path}/rscape_results.txt", dry)

def rscape_infernal():
    # Set up for R-scape analysis
    #exe('rm -f rscape_results.txt')
    #exe('rm -rf rscape_output')
    try:
        os.makedirs(f'{job_path}/rscape_infernal')
    except FileExistsError:
       pass
    exe(f"{RSCAPE_PATH} --outdir {job_path}/rscape_infernal --cacofold --outtree {job_path}/infernal.sto | tee {job_path}/rscape_infernal_results.txt", dry)

def is_hit():
    """
    Get the number of base pairs covered by the helices in the .helixcov file.
If the total number of base pairs covered is greater than or equal to 3 and there are at least 2 helices that cover 2 or more base pairs, return True.
    
    # RM_HELIX 27-45 148-162, nbp = 13 nbp_cov = 2
    # RM_HELIX 172-175 183-186, nbp = 4 nbp_cov = 0
    # RM_HELIX 187-192 373-377, nbp = 5 nbp_cov = 1
    # RM_HELIX 8-14 378-384, nbp = 7 nbp_cov = 1
    # RM_HELIX 396-406 469-479, nbp = 11 nbp_cov = 0
    # RM_HELIX 385-394 727-761, nbp = 10 nbp_cov = 0
    """
    import re
    import glob

    def analyze_nbp_cov(nbp_cov):
        """
        """
        ic(results)
        total_nbp_covs = sum(nbp_cov)
        print(f"Total nbp_covs: {total_nbp_covs}")
        if total_nbp_covs >= 3 and [i for i in nbp_cov if i >= 2]:
            print("Hit")
            return True
        return False
    
    def find_helixcov_file(folder_path):
        # Use glob to find all .helixcov files in the folder
        helixcov_files = glob.glob(os.path.join(folder_path + '/rscape_output', "*.helixcov"))

        if helixcov_files:
            # Return the path of the first .helixcov file found
            return helixcov_files[0]
        else:
            return None

    def parse_nbp_cov(text):
        pattern = r'# RM_HELIX .+, nbp = \d+ nbp_cov = (\d+)'
        matches = re.findall(pattern, text)

        results = []
        for nbp_cov in matches:
            results.append(int(nbp_cov))
        return results

    helixcov_file = find_helixcov_file(job_path)
    ic(helixcov_file)
    if helixcov_file:
        print(f"Found .helixcov file: {helixcov_file}")
        with open(helixcov_file, 'r') as file:
            content = file.read()
            results = parse_nbp_cov(content)
            return analyze_nbp_cov(results)

def infernal():
    import os
    import glob

    def find_cocofold_file(folder_path):
        # Use glob to find all .cocofold files in the folder
        cocofold_files = glob.glob(os.path.join(folder_path, 'rscape_output', "*.cacofold.sto")) # cocofold

        if cocofold_files:
            # Return the path of the first .cocofold file found
            return cocofold_files[0]
        else:
            return None

    # Example usage
    cocofold_file = find_cocofold_file(job_path)
    if cocofold_file:
        print(f"Found .cocofold file: {cocofold_file}")
    else:
        print("No .cocofold file found in the specified folder.")

    if not args.dev_skip_cmcalibrate:
        # Run cmbuild
        cm = os.path.splitext(cocofold_file)[0] + '.cm'
        cmd = f'cmbuild -F {cm} {cocofold_file}'
        print(cmd)
        exe(cmd, dry)
        # Run cmcalibrate
        cmd = f'cmcalibrate {cm}'
        print(cmd)
        exe(cmd, dry)
        pass
    # Search the Rfam database with the covariance model to eliminate known case
    #cmd = f'cmscan -o {j}/cmscan_res.out {RFAM_DB_PATH} {cm}'
    #print(cmd)
    #exe(cmd)
    cmd = f'cmsearch -A {job_path}/infernal.sto -o {job_path}/cmsearch.out {cm} {db}'
    print(cmd)
    exe(cmd, dry)
    rscape_infernal()
    
def save_to_slurm():
        name = f'{dbbase}X{fbase}'
        t = f"""#!/bin/bash

#SBATCH -n 3 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0 # Runtime in minutes
#SBATCH -p eddy # Partition to submit to
#SBATCH --mem=50000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append
#SBATCH -o {job_path}/%j.out # Standard out goes to this file
#SBATCH -e {job_path}/%j.err # Standard err goes to this filehostname
#SBATCH --job-name={job_path} # Name of the job
# Runs a command on all FASTA files in current directory
python {' '.join(sys.argv[:]).replace('--slurm', '')}
        """
        with open(f'{job_path}/run.slurm', 'w') as fi:
            fi.write(t)
        print(f'sbatch {job_path}/run.slurm')
        os.chmod(f'{job_path}/run.slurm', 0o755) #mode=stat.S_IXUSR)
        os.system(f'sbatch {job_path}/run.slurm')
        
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()
    db = args.db
    evalue = args.evalue
    nofinteractions = args.iteractions
    dry = args.dry
    CPUs = args.cpus

    if list != type(args.fasta):
        args.fasta = [args.fasta]

    for f in args.fasta:
        fbase = os.path.basename(f).replace('.fa', '')
        seq_fn = os.path.basename(f)
        dbbase = os.path.basename(args.db[0]).replace('.fa', '')
        db = ' '.join(args.db)

        if args.create_job_folder:
            if args.job_name:
                job_path = 'jobs/' + args.job_name
            else:
                job_path = 'jobs/' + fbase
            try:
                os.makedirs(f'{job_path}', exist_ok=True) 
            except FileExistsError:
                pass
            if f:
                shutil.copy(f, job_path)
                shutil.copy(args.flanked, job_path + '/' + fbase + '_flanked.fa')
        else:
            job_path = os.path.dirname(f)  # the folder where the fasta file is

        now()

        seq_path = job_path + '/' + seq_fn
        seq_flanked_path = job_path + '/' +  seq_fn.replace('.fa', '') + '_flanked.fa'

        print(args, flush=True)

        with open(f'{job_path}/cmd.sh', 'w') as fi:
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
        fh = logging.FileHandler(f'{job_path}/log.log')
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
        #if not args.rscape:
        # test for seq now! ;))))

        if args.flanked:
            with open(seq_path, 'r') as fn:
                lines = fn.readlines()
                sequence_only = "".join(line.strip().upper() for line in lines if not line.startswith('>'))
            ic(sequence_only)

            with open(args.flanked, 'r') as fn:
                lines = fn.readlines()
                flanked_sequence_only = "".join(line.strip().upper() for line in lines if not line.startswith('>'))
            ic(flanked_sequence_only)

            # Find the beginning index
            start_index = flanked_sequence_only.find(sequence_only)

            # Calculate the end index if substring is found
            if start_index != -1:
                end_index = start_index + len(sequence_only) - 1
                print(f"Substring found from index {start_index + 1} to {end_index}.")
            else:
                raise Exception("Substring not found.")

        if args.repeatmasker:
            # cmd REPEAT_MASKER_PATH
            cmd = f'{REPEAT_MASKER_PATH}/RepeatMasker {seq_path}'
            print(cmd)
            seq_masked_path = seq_path + '.masked'
            print(seq_masked_path)
            if not os.path.exists(seq_masked_path):
               print('RepeatMasker: No repetitive sequences were detected in seq')
            else:
               print('RepeatMasker: Repetitive sequences were detected in seq')
               query = seq_masked_path  
               print(open(seq_masked_path).read())

            # cmd REPEAT_MASKER_PATH
            cmd = f'{REPEAT_MASKER_PATH}/RepeatMasker {seq_flanked_path}'
            print(cmd)
            seq_flanked_masked_path = seq_flanked_path + '.masked'
            if not os.path.exists(seq_flanked_masked_path):
               print('RepeatMasker: No repetitive sequences were detected in seq flanked')
            else:
               print('RepeatMasker: Repetitive sequences were detected in seq flanked')
               seq_flanked_path = seq_masked_path  
               print(open(seq_flanked_masked_path).read())

        if not args.dev_skip_search:
            search(seq_path, seq_flanked_path)
        # Remove duplicate copies of genomes
        find_top_scoring_hits(job_path)
        if not args.dev_skip_rscape:
            rscape()
        is_hit = is_hit()
        if is_hit:
            if not args.dev_skip_infernal:
                infernal()
        else:   
            logging.info('No Infernal')
        logging.info('Normal termination of rnahub')
        logger.info('Normal termination of rnahub')
        print('Normal termination of rnahub', flush=True)
        now()
        
