# rnahub

```
./rnahub.py -h
usage: rnahub.py [-h] [--job-name JOB_NAME] [-v] [--slurm] [--evalue EVALUE] [--iteractions ITERACTIONS] [--dry] [--rscape] db file [file ...]

positional arguments:
  db
  file

options:
  -h, --help            show this help message and exit
  --job-name JOB_NAME   by default is input file name (wihout extension)
  -v, --verbose         be verbose
  --slurm               send it to slumrm
  --evalue EVALUE       e-value threshold
  --iteractions ITERACTIONS
                        number of iterations
  --dry                 show all cmds, dont run them
  --rscape              rscape only
```
# db

- `1409_Acomycota_genomes-may19.fa` - Gao, W., Jones, T. A. & Rivas, E. Discovery of 17 conserved structural RNAs in fungi. Nucleic Acids Res 49, gkab355- (2021).
- [for testing] `1409_Acomycota_genomes-may19_HEAD_G1.fa` - head -c 1G of 1409_Acomycota_genomes-may19.fa
- `nt.fa` - blast nt in fasta format
- `nt` - blast nt in blast format
- `bacteria_archaea.fa  `
- `invertebrates_fungi_fugu.fa` 
- `pdbnt`  
- `vertebrates_85.fa`
- `AY302558.1.fasta` 
- `invertebrate_genomes_db.fasta`  
- `nr ff`
- `rfam`  
- `wgao`

Gao, W., Jones, T. A. & Rivas, E. Discovery of 17 conserved structural RNAs in fungi. Nucleic Acids Res 49, gkab355- (2021).  
https://academic.oup.com/nar/article/49/11/6128/6292099
