"""
Snakefile for building STAR index for hg38
nakemake -c 16 -s rules/create_star_indices_hg38.smk --use-conda --config STAR_DIR='/Users/a1667917/Documents/Pipelines'

"""
import os

# load default config
configfile: os.path.join(workflow.basedir, '..',  'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]

if config['hg38_dir'] is None:
    hg38_dir = "/hpcfs/users/a1667917/STAR_Ref_Genomes"
else:
    hg38_dir = config["hg38_dir"]

# needs to be created before star is run
if not os.path.exists(os.path.join(hg38_dir, 'hg38_50')):
  os.makedirs(os.path.join(hg38_dir, 'hg38_50'))

if not os.path.exists(os.path.join(hg38_dir, 'hg38_75')):
  os.makedirs(os.path.join(hg38_dir, 'hg38_75'))

rule create_50_indecies:
    """create index."""
    input:
        os.path.join(hg38_dir, 'chr8.fasta'),
        os.path.join(hg38_dir, 'hg38.ncbiRefSeq.gtf')
    output:
        os.path.join(hg38_dir,"chr8_50_star.log")
    params:
        os.path.join(hg38_dir, 'hg38_50')
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
        os.path.join(hg38_dir, 'hg38_75')
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
