# config.py
import os
DEBUG = True
nhmmer = "/home/rnahub/opt/hmmer-3.3.2/bin/nhmmer"
RFAM_DB_PATH = "/home/rnahub/db/Rfam.cm"
RSCAPE_DIR = "/home/rnahub/opt/rscape/rscape_v2.0.5"
RSCAPE_PATH = os.path.join(RSCAPE_DIR, "bin/R-scape")
EASEL_PATH = '/home/rnahub/opt/hmmer-3.3.2/bin/'

if os.path.exists('config_local.py'):
    from config_local import *

