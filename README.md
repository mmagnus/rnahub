# rnahub

```
/rnahub.py -h
usage: rnahub.py [-h] [--db DB [DB ...]] [--job-name JOB_NAME] [-v] [--slurm] [--evalue EVALUE] [--lmin LMIN] [--evalue-final EVALUE_FINAL] [--iteractions ITERACTIONS]
                 [--cpus CPUS] [--dry] [--repeatmasker] [--job-folder JOB_FOLDER] [--dev-skip-search] [--dev-skip-nhmmer0] [--dev-skip-nhmmer123] [--dev-skip-cmcalibrate]
                 [--dev-skip-rscape] [--dev-skip-infernal] [--rscape] [--fasta FASTA] [--flanked FLANKED] [--flanks-in-header] [--flanks-start FLANKS_START]
                 [--flanks-end FLANKS_END]

job_path directory (j before)

.. warning ::

   R-scape version, searching for .helixcov files here!

Is hit now its true with 1 covariance.

Re-run it in given folder:

    ~/rnahub/rnahub.py --job-folder . --dev-skip-nhmmer123

optional arguments:
  -h, --help            show this help message and exit
  --db DB [DB ...]
  --job-name JOB_NAME   by default is input file name (wihout extension)
  -v, --verbose         be verbose
  --slurm               send it to slumrm
  --evalue EVALUE       e-value threshold for all the runs but the final one
  --lmin LMIN           esl-alimanip for v0 processing, default 50
  --evalue-final EVALUE_FINAL
                        e-value threshold for the final run
  --iteractions ITERACTIONS
                        number of iterations
  --cpus CPUS           number of cpus for nhmmer
  --dry                 show all cmds, dont run them
  --repeatmasker
  --job-folder JOB_FOLDER
                        create a job folder based on the path to the input fasta sequence, by default 'jobs/'so with example/seq.fa, the job folder is going to be jobs/seq/seq.fa
                        [and other files here]
  --dev-skip-search     skip v0..v3 all nhmmer searches
  --dev-skip-nhmmer0    show all cmds, dont run them
  --dev-skip-nhmmer123  show all cmds, dont run them
  --dev-skip-cmcalibrate
                        show all cmds, dont run them
  --dev-skip-rscape     show all cmds, dont run them
  --dev-skip-infernal   show all cmds, dont run them
  --rscape              rscape only
  --fasta FASTA         .fa for now, don't use .fasta
  --flanked FLANKED     .fa for now, don't use .fasta, flank the sequence including the query sequence
  --flanks-in-header    run flanked mode (create extra v0 files), syntax in the fasta header '><seq_name> <start>-<end>, use this or --flanked fasta file
  --flanks-start FLANKS_START
                        start of flank
  --flanks-end FLANKS_END
                        end of flank
```

# Install
See config_local.py to set up paths to your HMMER, R-Scape, INFERNAL.

     python -m pip install icecream

RC Harvard Cluster:

```
[mmagnus@holylogin08 ~]$ module load python/3.10.12-fasrc01
[mmagnus@holylogin08 ~]$ python
Python 3.10.12 | packaged by conda-forge | (main, Jun 23 2023, 22:40:32) [GCC 12.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>>
[mmagnus@holylogin08 ~]$ python -m pip install icecream
```

# Workflow


Pre-workflow: 

- run Blast on NR to narrow the searches with nhmmer to e.g., fungi gnomes* OR ask the user.
- add flanks!
  
*we need pre-compiled gnomes

![gkab355fig2](https://github.com/mmagnus/rnahub/assets/118740/f52da725-41f8-4224-aed4-7a393fe6432f)
Gao, W., Jones, T. A. & Rivas, E. Discovery of 17 conserved structural RNAs in fungi. Nucleic Acids Res 49, gkab355- (2021).
https://academic.oup.com/nar/article/49/11/6128/6292099

