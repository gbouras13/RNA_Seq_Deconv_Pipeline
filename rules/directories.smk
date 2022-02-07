"""
Database and output locations for Hecatomb
Ensures consistent variable names and file locations for the pipeline, the database download script,
and the addHost script.
"""


DBDIR = 'Databases'

### OUTPUT DIRECTORY
if config['Output'] is None:
    OUTPUT = 'wgs_tcga_out'
else:
    OUTPUT = config['Output']

if config['STAR_DIR'] is None:
    STAR_DIR = "/hpcfs/users/a1667917/STAR_Ref_Genomes/hg38"
else:
    STAR_DIR = config['STAR_DIR']

### DATABASE SUBDIRs
CONPATH = os.path.join(DBDIR, "contaminants")
TAX = os.path.join(DBDIR, "tax", "taxonomy")
TABLES = os.path.join(DBDIR, "tables")
HOSTPATH = os.path.join(DBDIR, "host")

### OUTPUT DIRs
RESULTS = os.path.join(OUTPUT, 'RESULTS')
WORKDIR = os.path.join(OUTPUT, 'PROCESSING')
TMP = os.path.join(WORKDIR, 'TMP')
LOGS = os.path.join(OUTPUT, 'LOGS')

# fastqc
FASTQC = os.path.join(RESULTS, "FASTQC")

# needs to be created before fastqc is run
if not os.path.exists(RESULTS):
  os.makedirs(RESULTS)
  # needs to be created before fastqc is run)
if not os.path.exists(FASTQC):
  os.makedirs(FASTQC)
