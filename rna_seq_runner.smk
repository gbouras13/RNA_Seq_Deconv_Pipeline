"""
The snakefile that runs the pipeline.
# run from login node
snakemake -c 16 -s RNA_Seq_Deconv_Pipeline/rna_seq_runner.smk --use-conda --config Reads=RNA_EGA_Bams Output=RNA_EGA_Out
# from computer node
snakemake -c 16 -s /Users/a1667917/Documents/Pipelines/RNA_Seq_Deconv_Pipeline/rna_seq_runner.smk --use-conda --config Reads=RNA_EGA_Bams Output=RNA_EGA_Out
# HPC
# on login node from pipeline dir
snakemake -s rna_seq_runner.smk -c 1 --use-conda --config Reads=Bams Output=test hg38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' --conda-create-envs-only --conda-frontend conda
# to run
snakemake -s rna_seq_runner.smk --use-conda --config Reads=TCGA_RNA_Total_Bams/ Output=RNA_EGA_Out HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' TCGA=True --profile wgs_tcga
"""


### DEFAULT CONFIG FILE
configfile: os.path.join(workflow.basedir,  'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]
# STAR can only use 16 threads on hpc for some reason
# https://github.com/alexdobin/STAR/issues/1074
MediumJobCpu = config["MediumJobCpu"]


### DIRECTORIES
include: "rules/directories.smk"

READS = config['Reads']
OUTPUT = config['Output']

# Parse the samples and read files
include: "rules/samples.smk"
sampleReads = parseSamples(READS)
SAMPLES = sampleReads.keys()

# TCGA or EGA 

TCGA = config['TCGA']

# Import rules and functions
include: "rules/targets.smk"
include: "rules/bam_to_fastq.smk"
include: "rules/qc.smk"
include: "rules/align.smk"

rule all:
    input:
        PreprocessingFiles,
        AlignFiles,
        FeatureCountFiles
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
