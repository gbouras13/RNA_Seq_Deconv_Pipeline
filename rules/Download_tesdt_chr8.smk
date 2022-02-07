"""
Snakefile for downloading and building STAR index for hg38
"""
import os

# load default config


configfile: os.path.join(workflow.basedir, '..',  'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]

if config['hg38_dir'] is None:
    hg38_dir = "/Users/a1667917/Documents/Pipelines"
else:
    hg38_dir = config["hg38_dir"]


# needs to be created before star is run
if not os.path.exists(os.path.join(hg38_dir, 'STAR_50')):
  os.makedirs(os.path.join(hg38_dir, 'STAR_50'))

if not os.path.exists(os.path.join(hg38_dir, 'STAR_75')):
  os.makedirs(os.path.join(hg38_dir, 'STAR_75'))

rule all:
    input:
        os.path.join(hg38_dir,"chr8_50_star.log"),
        os.path.join(hg38_dir,"chr8_75_star.log")

rule create_50_indecies:
    """create index."""
    input:
        os.path.join(hg38_dir, 'chr8.fasta'),
        os.path.join(hg38_dir, 'hg38.ncbiRefSeq.gtf')
    output:
        os.path.join(hg38_dir,"chr8_50_star.log")
    params:
        os.path.join(hg38_dir, 'STAR_50')
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
        os.path.join(hg38_dir, 'chr8.fasta'),
        os.path.join(hg38_dir, 'hg38.ncbiRefSeq.gtf')
    output:
        os.path.join(hg38_dir,"chr8_75_star.log")
    params:
        os.path.join(hg38_dir, 'STAR_75')
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
