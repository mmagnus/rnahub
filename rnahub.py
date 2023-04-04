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


def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    #parser.add_argument('-', "--", help="", default="")

    parser.add_argument("-v", "--verbose",
                        action="store_true", help="be verbose")
    parser.add_argument("file", help="", default="") # nargs='+')
    return parser


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if list != type(args.file):
        args.file = [args.file]

    import os

    os.system('figlet -f smblock rnahub')
    for f in args.file:
        fh = open(f)
        h = fh.readline().strip()
        seq = fh.readline().strip()
        print(h)
        print(seq)
        cmd = 'blastn -db /home/rnahub/rnahub/db/pdbnt -query ' + f
        print(cmd)
        os.system(cmd)

