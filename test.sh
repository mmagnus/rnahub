#!/bin/bash
eval "$(/home/rnahub/miniconda3/bin/conda shell.bash hook)"
conda activate RF2NA
#source /home/rnahub/miniconda3/bin/activate RF2NA
set -x
#python rnahub.py --blast pdbnt example/example.fa
python3 search_rfam.py --rfam /home/rnahub/rnahub/db/rfam/Rfam.cm --fasta example/example.fa --seed /home/rnahub/rnahub/db/rfam/Rfam.seed.gz

#python rnahub.py --blast refseq example/example.fa
##python rfam.py example/example.fa

