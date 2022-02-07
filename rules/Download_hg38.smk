"""
Snakefile for downloading STAR index for hg38
snakemake -c 16 -s rules/Download_hg38.smk --config hg38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes'
"""
import os

# load default config
configfile: os.path.join(workflow.basedir, '..', 'config', 'config.yaml')

BigJobMem = config["BigJobMem"]
BigJobCpu = config["BigJobCpu"]

if config['hg38_dir'] is None:
    hg38_dir = "/hpcfs/users/a1667917/STAR_Ref_Genomes"
else:
    hg38_dir = config["hg38_dir"]

# needs to be created before (should exist)
if not os.path.exists(os.path.join(hg38_dir)):
  os.makedirs(os.path.join(hg38_dir))


rule all:
    input:
        os.path.join(hg38_dir, 'download_hg38_fasta.dlflag'),
        os.path.join(hg38_dir, 'download_hg38_gtf.dlflag')
        # ,
        # os.path.join(hg38_dir, 'unzip.dlflag')

rule download_hg_38_fasta:
    """Rule to Download hg38 fasta."""
    output:
        os.path.join(hg38_dir,'download_hg38_fasta.dlflag'),
        os.path.join(hg38_dir, 'hg38.fa')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        wget -c "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" -O hg38.fa
        mv "hg38.fa" {output[1]}
        touch {output[0]}
        """

rule download_hg_38:
    """Rule to Download hg38 gtf."""
    output:
        os.path.join(hg38_dir,'download_hg38_gtf.dlflag'),
        os.path.join(hg38_dir, 'hg38.ncbiRefSeq.gtf')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        wget -c "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz" -O hg38.ncbiRefSeq.gtf
        mv "hg38.ncbiRefSeq.gtf" {output[1]}
        touch {output[0]}
        """
#
# rule unzip:
#     """gunzip files."""
#     input:
#         os.path.join(hg38_dir, 'hg38.fa.gz'),
#         os.path.join(hg38_dir, 'hg38.ncbiRefSeq.gtf.gz')
#     output:
#         os.path.join(hg38_dir, 'hg38.fa'),
#         os.path.join(hg38_dir, 'hg38.ncbiRefSeq.gtf'),
#         os.path.join(hg38_dir, 'unzip.dlflag')
#     threads:
#         1
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         gunzip {input[0]} -C {output[0]}
#         gunzip {input[1]} -C {output[1]}
#         touch {output[2]}
#         """
