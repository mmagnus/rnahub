#!/home/rnahub/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
from Bio import AlignIO
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')
import os
from config import RFAM_DB_PATH, RSCAPE_PATH, EASEL_PATH, RFAM_FILES
import re

def rscape(alignment, dry=False):
    try:
        os.makedirs(f'{job_path}/rscape_rfam')
    except FileExistsError:
       pass
    exe(f"{RSCAPE_PATH} --outdir {job_path}/rscape_rfam --cacofold --outtree {alignment} | tee {job_path}/rscape_rfam_results.txt", dry)

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

    
def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--job-path', default='.', help='Path to the job directory')
    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


def coverge(alignment_fn, verbose=False):
    # Load the alignment
    alignment = AlignIO.read(alignment_fn, "stockholm")
    alignment_name = alignment_fn.split('/')[-1]
    target_id = alignment[0].id
    
    # Find the target sequence
    target_seq = None
    for record in alignment:
        if record.id == target_id:
            target_seq = record.seq
            break

    if target_seq is None:
        print(f"Target sequence '{target_id}' not found in the alignment.")
    else:
        seq_len = 0
        cols_no_coverage = 0
        # Collect columns where the target sequence has no gaps
        filtered_columns = []
        for col_index in range(alignment.get_alignment_length()):
            if target_seq[col_index] != "-":
                column = [record.seq[col_index] for record in alignment]
                filtered_columns.append(column)
                seq_len += 1

        # Display filtered columns
        if verbose: print(f"Columns where the target sequence '{target_id}' has no gaps:")
        for col in filtered_columns:
            if verbose: print(col, len(col))
            # Count the number of '-'
            gap_count = col.count('-')
            if verbose: print(f"Number of '-': {gap_count}")
            ratio_per_col = gap_count/(len(col) - 1)
            if verbose: print(f'ratio of gaps: {ratio_per_col}') # -1 target_sequence
            if ratio_per_col > 0.5:
                if verbose: print('WARNING: more than 50% gaps in column', cols_no_coverage)
                cols_no_coverage +=1
            if verbose: print()

        ratio = round((1 - (cols_no_coverage / seq_len)) * 100)
        print(f"Query {target_id} length {seq_len} in the alignment {alignment_name} number of columns with low coverage (<50%) for the query sequence {cols_no_coverage}; coverge for the query: {ratio}%") 

        
def process_rfam(family_id, seq):
    """seq ../rnahub-web/media/jobs/xrRNA-925d3bfd/seq.fa """
    # Concatenate files
    # ../rnahub-web/media/jobs/xrRNA-925d3bfd/seq(.fa) -> seq+
    seq_seed = seq.replace('.fa', f'+{family_id}_seed.fa')

    with open(seq_seed, "w") as outfile:
        content = open(seq).read()
        content += '\n' + open(f'{RFAM_FILES}/{family_id}.seed.fa').read()  # get seed fa
        #print(content)
        outfile.write(content)
        print('saved to', seq_seed)

    c = 0
    for l in open(seq_seed):
        if l.startswith('>'):
            c += 1
            if c > 5:
                break
        print(l.strip())
    print('[...]')        
    # Run cmalign
    alignment = seq.replace('.fa', f'+{family_id}_seed.sto')
    cmd = ' '.join(["cmalign", f"{RFAM_FILES}/{family_id}.cm", seq_seed, ">", alignment])
    print(cmd)
    #subprocess.run(cmd)
    os.system(cmd)
    return alignment

    
def get_rfam_accession(model_name):
    """
    Fetch the Rfam accession ID for a given model name.
    """
    import requests
    url = f"https://rfam.org/family/{model_name}/acc"
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text.strip()
    else:
        return f"Error: Unable to fetch data (status code: {response.status_code})"

    
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    os.chdir('/home/rnahub/rnahub-web/')

    if list != type(args.file):
        args.file = [args.file]

    for filename in args.file:

            db = RFAM_DB_PATH
            #-o 
            d = os.path.dirname(filename)
            if os.path.exists(d + "/cmscan.out"):
                os.remove(d + "/cmscan.out")                
            cmd = "cmscan --tblout " + d + "/cmscan.tblout -o " + d + "/cmscan.out --noali " + db + " " + filename
            #families=`grep -v '^#' cmscan.tblout | head -n $max_rfam_num | uniq | awk '{print $2}' | sed -z 's/\n/|/g;s/|$/\n/'`
            #echo "Rfam families:" $families
            ic(cmd)
            os.system(cmd)
            with open(d + "/cmscan.out") as f:
                text = f.read()
                            # Define the regular expression
 #            pattern = r"Hit scores:(.*?)------ inclusion threshold ------"
 #            # Extract the match
 #            match = re.search(pattern, text, re.DOTALL)

 #            # Process and print the result
 #            result = match.group(1).strip().replace("""rank     E-value  score  bias  modelname     start    end   mdl trunc   gc  description
 # ----   --------- ------ -----  ------------ ------ ------   --- ----- ----  -----------""", '').strip()  # Extract and clean the matched content
 #            print(result)
 #            print("No match found.")
            hits = []
            for line in text.split("\n"):
                if 'rank' in line:
                    header = line
                # if '-----------' in line:
                #     header += '\n' + line
                if 'cm' in line:
                    if line.split()[1] == '!':
                        hits.append(line)
            if len(hits):
                print(header)
                for hit in hits:
                    print(hit)
                    # ['(1)', '!', '7e-09', '51.8', '0.0', 'PLRV_xrRNA',

            if hits:      
                top_hit_modelname = hits[0].split()[5]
                rfam_accession = get_rfam_accession(top_hit_modelname)
                print()
                print(f"Rfam accession for {top_hit_modelname}: {rfam_accession}")
                job_path = args.job_path
                aligment = process_rfam(rfam_accession, filename)

                cmd = ''.join([f'{EASEL_PATH}/esl-alistat ', aligment ,' > ', job_path, '/rfam_stats.txt'])
                print(cmd)
                exe(cmd)

                print()
                coverge(aligment)
                print()
                rscape(aligment)

            
