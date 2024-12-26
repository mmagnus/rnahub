#!/home/rnahub/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function
import argparse
from icecream import ic
import sys
ic.configureOutput(outputFunction=lambda *a: print(*a, file=sys.stderr))
ic.configureOutput(prefix='> ')
import os
from config import RFAM_DB_PATH
import re

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="", nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    for f in args.file:

            db = RFAM_DB_PATH
            #-o 
            d = os.path.dirname(f)
            if os.path.exists(d + "/cmscan.out"):
                os.remove(d + "/cmscan.out")                
            cmd = "cmscan --tblout " + d + "/cmscan.tblout -o " + d + "/cmscan.out --noali " + db + " " + f
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

            #print(text)
                
