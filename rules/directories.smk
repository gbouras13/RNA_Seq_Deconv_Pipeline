"""
Database and output locations for the pipeline
"""


DBDIR = 'Databases'
TRUST4DIR = 'Trust_4_Files'

### OUTPUT DIRECTORY
if config['Output'] is None:
  OUTPUT = 'rna_tcga_out'
else:
  OUTPUT = config['Output']

if config['HG38_dir'] is None:
  HG38_dir = '/hpcfs/users/a1667917/STAR_Ref_Genomes'
else:
  HG38_dir = config['HG38_dir']

if config['Salmon_dir'] is None:
  Salmon_dir = '/hpcfs/users/a1667917/Salmon_Ref_Genomes'
else:
  Salmon_dir = config['Salmon_dir']

if config['Kallisto_dir'] is None:
  Kallisto_dir = '/hpcfs/users/a1667917/Kallisto_Ref_Genomes'
else:
  Kallisto_dir = config['Kallisto_dir']


### OUTPUT DIRs
RESULTS = os.path.join(OUTPUT, 'RESULTS')
WORKDIR = os.path.join(OUTPUT, 'PROCESSING')
TMP = os.path.join(WORKDIR, 'TMP')
LOGS = os.path.join(OUTPUT, 'LOGS')
STAR_BAMS = os.path.join(RESULTS, 'STAR_BAMS')
SALMON_OUTPUT = os.path.join(RESULTS, 'SALMON_OUTPUT')
KALLISTO_OUTPUT = os.path.join(RESULTS, 'KALLISTO_OUTPUT')

TRUST4 = os.path.join(TRUST4, 'STAR_BAMS')


# fastqc
FASTQC = os.path.join(RESULTS, "FASTQC")
MULTIQC = os.path.join(FASTQC, "MULTIQC")

# needs to be created before fastqc is run
if not os.path.exists(RESULTS):
  os.makedirs(RESULTS)
# needs to be created before fastqc is run)
if not os.path.exists(FASTQC):
  os.makedirs(FASTQC)
# needs to be created before fastqc is run)
if not os.path.exists(MULTIQC):
  os.makedirs(MULTIQC)
# needs to be created before alignment
if not os.path.exists(STAR_BAMS):
  os.makedirs(STAR_BAMS)
# needs to be created for fastqc 
if not os.path.exists(TMP):
  os.makedirs(TMP)
# needs to be created for kallisto 
if not os.path.exists(KALLISTO_OUTPUT):
  os.makedirs(KALLISTO_OUTPUT)