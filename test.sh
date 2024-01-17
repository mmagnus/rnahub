#!/bin/bash
# ./nextflow rnahub.nf

echo 'fasta: example/random.fa' > config.yaml #example.fa
#echo 'fasta: example/example.fa' > config.yaml

function test1 {
   eval "$(/home/rnahub/miniconda3/bin/conda shell.bash hook)"
   conda activate RF2NA
   python ./search_blast.py --db pdbnt example/random.fa
   }

function test2 {       
   eval "$(/home/rnahub/miniconda3/bin/conda shell.bash hook)"
   conda activate RF2NA
   snakemake -p --cores 2 -F ##
   #./search_snakemake.py
}

test2

# eval "$(/home/rnahub/miniconda3/bin/conda shell.bash hook)"
# conda activate RF2NA
# exit

# #!/bin/bash
# eval "$(/home/rnahub/miniconda3/bin/conda shell.bash hook)"
# conda activate RF2NA
# #source /home/rnahub/miniconda3/bin/activate RF2NA
# set -x

# python ./search_blast.py --db pdbnt example/example.fa
# python ./search_blast.py --db refseq example/example.fa
# #python ./search_blast.py --db nt example/example.fa
# python search_rfam.py --rfam /home/rnahub/rnahub/db/rfam/Rfam.cm --fasta example/example.fa --seed /home/rnahub/rnahub/db/rfam/Rfam.seed.gz
