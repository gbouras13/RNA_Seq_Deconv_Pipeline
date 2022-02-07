"""
Database and output locations for the pipeline
"""


DBDIR = 'Databases'

### OUTPUT DIRECTORY
if config['Output'] is None:
  OUTPUT = 'rna_tcga_out'
else:
  OUTPUT = config['Output']

if config['hg38_dir'] is None:
  HG38_dir = '/hpcfs/users/a1667917/STAR_Ref_Genomes'
else:
  HG38_dir = config['hg38_dir']


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
# needs to be created before alignment
if not os.path.exists(STAR_BAMS):
  os.makedirs(STAR_BAMS)
# needs to be created for fastqc 
if not os.path.exists(TMP):
  os.makedirs(TMP)