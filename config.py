import os
DEBUG = True
nhmmer = "/home/rnahub/opt/hmmer/hmmer-3.3.2/bin/nhmmer"
RFAM_DB_PATH = "/home/rnahub/db/rfam/Rfam.cm"
RSCAPE_DIR = "/home/rnahub/opt/rscape/rscape_v2.5.2"
RSCAPE_PATH = os.path.join(RSCAPE_DIR, "bin/R-scape")
EASEL_PATH = '/home/rnahub/opt/hmmer/hmmer-3.3.2/bin/'
INFERNAL_PATH = '/usr/bin/'
REPEAT_MASKER_PATH = '/home/rnahub/opt/RepeatMasker'
RFAM_FILES = "/home/rnahub/rnahub/rfam"
try:
    from config_local import *
except ImportError:
    pass

