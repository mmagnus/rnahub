# rnahub

```
(base) rnahub@rnahub:~/rnahub$ ./rnahub.py -h
usage: rnahub.py [-h] [--db DB] [--job-name JOB_NAME] [-v] [--slurm] [-f] [--evalue EVALUE] [--evalue-final EVALUE_FINAL] [--iteractions ITERACTIONS] [--cpus CPUS] [--dry] [--dev-skip-nhmmer0] [--dev-skip-nhmmer123]
                 [--dev-skip-cmcalibrate] [--dev-skip-rscape] [--dev-skip-infernal] [--rscape] [--fasta FASTA [FASTA ...]]

j job directory

.. warning ::

   R-scape version, searching for .helixcov files here!

optional arguments:
  -h, --help            show this help message and exit
  --db DB
  --job-name JOB_NAME   by default is input file name (wihout extension)
  -v, --verbose         be verbose
  --slurm               send it to slumrm
  -f, --flanked         run flanked mode (create extra v0 files), syntax in the fasta header '><seq_name> <start>-<end>
  --evalue EVALUE       e-value threshold for all the runs but the final one
  --evalue-final EVALUE_FINAL
                        e-value threshold for the final run
  --iteractions ITERACTIONS
                        number of iterations
  --cpus CPUS           number of cpus for nhmmer
  --dry                 show all cmds, dont run them

  --rscape              rscape only
  --fasta FASTA [FASTA ...]
                        .fa for now, don't use .fasta
						
```
# Workflow

Pre-workflow: 

- run Blast on NR to narrow the searches with nhmmer to e.g., fungi gnomes* OR ask the user.
- add flanks!
  
*we need pre-compiled gnomes

![gkab355fig2](https://github.com/mmagnus/rnahub/assets/118740/f52da725-41f8-4224-aed4-7a393fe6432f)
Gao, W., Jones, T. A. & Rivas, E. Discovery of 17 conserved structural RNAs in fungi. Nucleic Acids Res 49, gkab355- (2021).
https://academic.oup.com/nar/article/49/11/6128/6292099

