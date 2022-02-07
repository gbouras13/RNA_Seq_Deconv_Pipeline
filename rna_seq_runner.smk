"""
The snakefile that runs the pipeline.
# run from login node
snakemake -c 16 -s RNA_Seq_Deconv_Pipeline/rna_seq_runner.smk --use-conda --config Reads=RNA_EGA_Bams Output=RNA_EGA_Out
# from computer node
snakemake -c 16 -s /Users/a1667917/Documents/Pipelines/RNA_Seq_Deconv_Pipeline/rna_seq_runner.smk --use-conda --config Reads=RNA_EGA_Bams Output=RNA_EGA_Out
snakemake -c 16 -s /Users/a1667917/Documents/Pipelines/RNA_Seq_Deconv_Pipeline/rna_seq_runner.smk --use-conda --config Reads=RNA_EGA_Bams Output=RNA_EGA_Out STAR_DIR='/Users/a1667917/Documents/Pipelines'
"""


### DEFAULT CONFIG FILE
configfile: os.path.join(workflow.basedir,  'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]

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
