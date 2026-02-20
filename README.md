# rnahub

rnahub orchestrates iterative nhmmer/Infernal/R-Scape searches for discovering new RNA homologs. 

The CLI wires common structural-RNA tools, organizes job folders, and captures diagnostics so you can re-run specific pipeline stages with minimal bookkeeping.

For the webserver, go to http://rnahub.org

(databases for now are not included)

## Highlights
- Launch nhmmer searches (v0–v3) with automatic job-folder management.
- Run optional steps such as R-Scape, Infernal, RepeatMasker, or flanked-search
  generation without retooling commands.
- Resume work by skipping completed stages, or run in `--dry` mode to inspect
  commands before execution.
- Designed for workstation runs as well as scheduled jobs on SLURM clusters.

## Installation

### Python dependency
Install Python requirements that are not bundled with the standard library. At
the moment only [`icecream`](https://github.com/gruns/icecream) is needed:

```bash
python -m pip install icecream  # installs into the Python you are using
```

### Configure tool paths
Create `config_local.py` in the project root to override paths from
`config.py`. Required entries point to your nhmmer, Infernal, R-Scape, and
auxiliary binaries:

```python
nhmmer = "/home/rnahub/opt/hmmer/hmmer-3.3.2/bin/nhmmer"
RSCAPE_PATH = "/home/rnahub/opt/rscape/rscape_v2.5.2/bin/R-scape"
EASEL_PATH = "/home/rnahub/opt/hmmer/hmmer-3.3.2/bin/"
INFERNAL_PATH = "/usr/bin/"  # used for cmalign
REPEAT_MASKER_PATH = "/home/rnahub/opt/RepeatMasker"
# Optional Rfam resources
RFAM_FILES = "/home/rnahub/rnahub/rfam"
RFAM_DB_PATH = "/home/rnahub/db/rfam/Rfam.cm"
```

Any variable omitted from `config_local.py` falls back to the defaults provided
in `config.py`.

### Harvard RC cluster snippet
Example session for the FAS RC environment:

```bash
module load python/3.10.12-fasrc01
python -m pip install --user icecream
```

## Usage

1. Prepare an input `.fa/.fasta` sequence or `.sto` alignment.
2. (Optional) Pre-filter the search space using BLAST against NR to focus on
   a taxonomic group. Add flanking regions if needed.
3. Run the pipeline, e.g.:

```bash
python rnahub.py \
  --input example/seq.fa \
  --db db/fungi.cm \
  --job-name seq_run \
  --job-folder jobs/seq_run \
  --cpus 8
```

Set `--slurm` to emit SLURM submission scripts, `--dry` to preview commands,
and pass `--dev-skip-*` switches to resume from intermediate stages (for
instance `--dev-skip-nhmmer123` to skip later nhmmer rounds).

When re-running inside a finished job folder:

```bash
python rnahub.py --job-folder . --dev-skip-nhmmer123
```

### Flanked mode
Add flanks via either:
- `--flanked path/to/flanks.fa` (FASTA only, not `.fasta`).
- `--flanks-in-header` where headers follow `><seq_name> <start>-<end>`.

## Command-line reference

```
$ python rnahub.py -h
usage: rnahub.py [-h] [--db DB [DB ...]] [--job-name JOB_NAME] [--verbose] [-v] [--slurm]
                 [--evalue EVALUE] [--lmin LMIN] [--evalue-final EVALUE_FINAL]
                 [--iteractions ITERACTIONS] [--cpus CPUS] [--dry] [--repeatmasker] [--utot]
                 [--job-folder JOB_FOLDER] [--dev-skip-search] [--dev-skip-nhmmer0]
                 [--dev-skip-nhmmer123] [--dev-skip-rscape] [--dev-skip-infernal]
                 [--rscape-path RSCAPE_PATH] [--rscape] [--input INPUT [INPUT ...]]
                 [--flanked FLANKED] [--flanks-in-header] [--flanks-start FLANKS_START]
                 [--flanks-end FLANKS_END]

optional arguments:
  -h, --help            show this help message and exit
  --db DB [DB ...]      databases (.cm/.sto) searched by nhmmer
  --job-name JOB_NAME   defaults to the input file stem
  --verbose             increase logging
  -v, --version         print the version based on the latest Git tag
  --slurm               submit to SLURM instead of running locally
  --evalue EVALUE       e-value threshold for all runs except the final one
  --lmin LMIN           `esl-alimanip` minimum length (default 50)
  --evalue-final EVALUE_FINAL
                        e-value threshold for the last round
  --iteractions ITERACTIONS
                        number of nhmmer iterations
  --cpus CPUS           nhmmer CPU count
  --dry                 show commands without executing them
  --repeatmasker        run RepeatMasker preprocessing
  --utot                treat sequences as DNA (U→T) for RepeatMasker
  --job-folder JOB_FOLDER
                        destination folder (default `jobs/<name>/<input>`)
  --dev-skip-search     skip v0–v3 nhmmer searches
  --dev-skip-nhmmer0    skip only the v0 nhmmer search
  --dev-skip-nhmmer123  skip v1–v3 nhmmer searches
  --dev-skip-rscape     skip the R-Scape stage
  --dev-skip-infernal   skip Infernal stage
  --rscape-path RSCAPE_PATH
                        override R-Scape path
  --rscape              run only the R-Scape step
  --input INPUT [INPUT ...]
                        FASTA/FA or Stockholm input sequences/alignments
  --flanked FLANKED     supply a FASTA with flanking regions
  --flanks-in-header    parse `><seq start-end` spans from headers
  --flanks-start FLANKS_START
                        left flank coordinate
  --flanks-end FLANKS_END
                        right flank coordinate
```

## Workflow tips
- Run BLAST on NR before nhmmer to scope the search space (e.g., fungi).
- Add flanks to recover context for hits.
- Keep a cache of pre-compiled genomes for faster iterative runs.
- Inspect `.helixcov` files produced by newer R-Scape versions.

![gkab355fig2](https://github.com/mmagnus/rnahub/assets/118740/f52da725-41f8-4224-aed4-7a393fe6432f)
Gao, W., Jones, T. A. & Rivas, E. Discovery of 17 conserved structural RNAs in fungi.
Nucleic Acids Res 49, gkab355 (2021). https://academic.oup.com/nar/article/49/11/6128/6292099
