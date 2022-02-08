"""
Snakefile for downloading files for Salmon indices
snakemake -c 1 -s rules/Download_Salmon_transcripts.smk --config Salmon_dir='/hpcfs/users/a1667917/Salmon_Ref_Genomes' --use-conda --conda-frontend conda
"""
import os

# load default config
configfile: os.path.join(workflow.basedir, '..', 'config', 'config.yaml')

if config['Salmon_dir'] is None:
    Salmon_dir = "/hpcfs/users/a1667917/Salmon_Ref_Genomes"
else:
    Salmon_dir = config["Salmon_dir"]

# needs to be created before (should exist)
if not os.path.exists(os.path.join(Salmon_dir)):
  os.makedirs(os.path.join(Salmon_dir))

rule all:
    input:
        os.path.join(Salmon_dir,'download_salmon_indices.dlflag'),
        os.path.join(Salmon_dir, 'genome_config.yaml') 

rule download_salmon_index:
    """Rule to Download salmon transcripts."""
    output:
        os.path.join(Salmon_dir, 'genome_config.yaml') ,
        os.path.join(Salmon_dir,'download_salmon_indices.dlflag')
    conda:
        os.path.join('..', 'envs','refgenie.yaml')
    threads:
        1
    shell:
        """
        refgenie init -c {output[0]}
        refgenie pull hg38/salmon_sa_index -c {output[0]} --pull-large
        refgenie pull hg38/tgMap -c {output[0]}
        touch {output[1]}
        """


