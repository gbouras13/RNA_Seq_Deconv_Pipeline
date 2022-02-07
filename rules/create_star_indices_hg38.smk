"""
Snakefile for building STAR index for hg38
snakemake -c 16 -s rules/create_star_indices_hg38.smk --use-conda --config HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes'

snakemake -c 16 -s rules/create_star_indices_hg38.smk --use-conda --config HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' --conda-create-envs-only --conda-frontend conda

"""
import os

# load default config
configfile: os.path.join(workflow.basedir, '..',  'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]

if config['HG38_dir'] is None:
    HG38_dir = "/hpcfs/users/a1667917/STAR_Ref_Genomes"
else:
    HG38_dir = config["HG38_dir"]

# needs to be created before star is run
if not os.path.exists(os.path.join(HG38_dir, 'hg38_50')):
  os.makedirs(os.path.join(HG38_dir, 'hg38_50'))

if not os.path.exists(os.path.join(hg38_dir, 'hg38_75')):
  os.makedirs(os.path.join(HG38_dir, 'hg38_75'))


rule all:
    input:
        os.path.join(HG38_dir, "chr8_50_star_build.log"),
        os.path.join(HG38_dir, "chr8_75_star_build.log")

rule create_50_indecies:
    """create index."""
    input:
        os.path.join(HG38_dir, 'hg38.fa'),
        os.path.join(HG38_dir, 'hg38.ncbiRefSeq.gtf')
    output:
        os.path.join(HG38_dir,"chr8_50_star_build.log")
    params:
        os.path.join(HG38_dir, 'hg38_50')
    threads:
        BigJobCpu
    conda:
        os.path.join('..', 'envs','align.yaml')
    resources:
        mem_mb=BigJobMem
    shell:
        """
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params[0]} \
        --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} \
        --sjdbOverhang 49
        touch {output[0]}
        """

rule create_75_indecies:
    """create index."""
    input:
        os.path.join(HG38_dir, 'hg38.fa'),
        os.path.join(HG38_dir, 'hg38.ncbiRefSeq.gtf')
    output:
        os.path.join(HG38_dir,"chr8_75_star_build.log")
    params:
        os.path.join(HG38_dir, 'hg38_75')
    threads:
        BigJobCpu
    conda:
        os.path.join('..', 'envs','align.yaml')
    resources:
        mem_mb=BigJobMem
    shell:
        """
        STAR --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params[0]} \
        --genomeFastaFiles {input[0]} --sjdbGTFfile {input[1]} \
        --sjdbOverhang 74
        touch {output[0]}
        """
