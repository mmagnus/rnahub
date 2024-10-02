# config.py
import os
DEBUG = True
nhmmer = "/n/home06/mmagnus/m/opt/rscape_v2.0.4.b/lib/hmmer/src/nhmmer"
db = "/n/home06/mmagnus/m/rnahub-db/vertebrates_185.fa"
db = "/n/home06/mmagnus/m/rnahub-db/AY302558.1.fasta" #Echovirus_E6.fa"
RFAM_DB_PATH = "db/Rfam.cm"
RSCAPE_DIR = "/n/home06/mmagnus/m/opt/rscape_v2.0.4.b/"
RSCAPE_PATH = os.path.join(RSCAPE_DIR, "bin/R-scape")
EASEL_PATH = ''
SCRIPTS_DIR = ''
if os.path.exists('config_local.py'):
    from config_local import *
