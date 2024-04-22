# config.py
import os
DEBUG = True
DATABASE_URI = 'mysql://user:password@localhost/mydatabase'
nhmmer = "/n/home06/mmagnus/m/opt/rscape_v2.0.4.b/lib/hmmer/src/nhmmer"
db = "/n/home06/mmagnus/m/rnahub-db/vertebrates_185.fa"
db = "/n/home06/mmagnus/m/rnahub-db/AY302558.1.fasta" #Echovirus_E6.fa"
rscapedir = "/n/home06/mmagnus/m/opt/rscape_v2.0.4.b/"
rscape = os.path.join(rscapedir, "bin/R-scape")

from config_local import *
