# RNA_Seq_Pipeline
Snakemake Pipeline to Extract and Align RNA-Seq Data from bam files (TCGA/ENA)

#### This is a Work in Progress.

* input fastqc to be fixed

* Some code (namely, the sample parsing in samples.smk) and inspiration for the general structure has been borrowed from https://github.com/shandley/hecatomb
* Inputs required are the relevant Bam files, which all must be placed in a certain directory - ideally in the Bams/ directory, but can be anywhere (the directory they are in must be specified with Reads={directory}).
* Only software requirement is that snakemake be in the $PATH.
* All snakemake commands are assumed to be run from the pipeline directory (such that for offline use the conda envs are available)

# Usage

1. The hg38 fasta and gtf files need to be downloaded (on login node) to the directory specified with HG38_dir:
```console
snakemake -c 1 -s rules/Download_hg38.smk --config HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes'
```

2. Download the conda envs for offline use:
```console
snakemake -c 1 -s rna_seq_runner.smk --use-conda --config Reads=Bams Output=test hg38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' --conda-create-envs-only --conda-frontend conda
```

3. The STAR indices need to be built - a 95 bp (TCGA), 150 (EGA) read length index will be created (run on the compute node through slurm):

* TCGA consists of 2x48 paired end reads (after trimming - av read length is 95 bp), EGA is 2x75 paired end reads


```console
snakemake -s rules/create_star_indices_hg38.smk --use-conda --config HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' --profile tcga_wgs
```

4. Run the pipeline

```console
snakemake -s rna_seq_runner.smk --use-conda --config Reads={path_to}/Bams Output=RNA_TCGA_Out HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' TCGA=True --profile wgs_tcga
```

* With a Slurm profile (see https://snakemake.readthedocs.io/en/stable/executing/cli.html https://github.com/Snakemake-Profiles/slurm https://fame.flinders.edu.au/blog/2021/08/02/snakemake-profiles-updated)
* You will need to cd to the pipeline directory in your jobscript before running if you want to run this offline (to use the premade conda envs)
* Use the TCGA=True flag for TCGA bams, TCGA=False for ENA

#### On slurm
```console
# creates the indices
sbatch rna_tcga_create_indices.sh
# runs the pipeline
sbatch rna_tcga_align.sh
```
