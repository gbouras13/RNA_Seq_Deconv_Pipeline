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
STAR_BAMS = os.path.join(RESULTS, 'STAR_BAMS')

# fastqc
FASTQC = os.path.join(RESULTS, "FASTQC")

# needs to be created before fastqc is run
if not os.path.exists(RESULTS):
  os.makedirs(RESULTS)
  # needs to be created before fastqc is run)
if not os.path.exists(FASTQC):
  os.makedirs(FASTQC)

 if not os.path.exists(STAR_BAMS):
   os.makedirs(STAR_BAMS)
