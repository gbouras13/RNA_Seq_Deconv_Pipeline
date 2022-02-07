"""
The snakefile that runs the pipeline.
# run from login node
snakemake -c 16 -s RNA_Seq_Deconv_Pipeline/rna_seq_runner.smk --use-conda --config Reads=RNA_EGA_Bams Output=RNA_EGA_Out
# from computer node
snakemake -c 16 -s /Users/a1667917/Documents/Pipelines/RNA_Seq_Deconv_Pipeline/rna_seq_runner.smk --use-conda --config Reads=RNA_EGA_Bams Output=RNA_EGA_Out
# HPC

snakemake -s rna_seq_runner.smk --use-conda --config Reads=RNA_EGA_Bams Output=RNA_EGA_Out hg38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' --conda-create-envs-only --conda-frontend conda
# to run
snakemake -s rna_seq_runner.smk --use-conda --config Reads=RNA_EGA_Bams Output=RNA_EGA_Out hg38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' --profile wgs_tcga
"""


### DEFAULT CONFIG FILE
configfile: os.path.join(workflow.basedir,  'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]

if config["hg_38_dir"] is None:
    hg38_dir = "/hpcfs/users/a1667917/STAR_Ref_Genomes"
else:
    hg38_dir = config["hg38_dir"]

### DIRECTORIES
include: "rules/directories.smk"

READS = config['Reads']
OUTPUT = config['Output']

# Parse the samples and read files
include: "rules/samples.smk"
sampleReads = parseSamples(READS)
SAMPLES = sampleReads.keys()

# Import rules and functions
include: "rules/targets.smk"
include: "rules/bam_to_fastq.smk"
include: "rules/qc.smk"
include: "rules/align.smk"

rule all:
    input:
        PreprocessingFiles
        # ## Assembly
        # AssemblyFiles,
        # ## Translated (nt-to-aa) search
        # SecondarySearchFilesAA,
        # ## Untranslated (nt-to-nt) search
        # SecondarySearchFilesNT,
        # ## Contig annotation
        # ContigAnnotFiles,
        # ## Mapping (read-based contig id)
        # MappingFiles,
        # ## Summary
        # SummaryFiles
