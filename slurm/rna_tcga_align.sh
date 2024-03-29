#!/bin/bash -l

#SBATCH --job-name=rna_tcga_snkg
#SBATCH --mail-user=george.bouras@adelaide.edu.au
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --err="rna_tcga_align_snk.err"
#SBATCH --output="rna_tcga_align_snk.out"

# Resources allocation request parameters
#SBATCH -p batch
#SBATCH -N 1                                                    # number of tasks (sequential job starts 1 task) (check this if your job unexpectedly uses 2 nodes)
#SBATCH -c 1                                                    # number of cores (sequential job calls a multi-thread program that uses 8 cores)
#SBATCH --time=2-23:00:00                                         # time allocation, which has the format (D-HH:MM), here set to 1 hou                                           # generic resource required (here requires 1 GPUs)
#SBATCH --mem=1GB                                              # specify memory required per node


SNK_DIR="/hpcfs/users/a1667917/Kevin/RNA_Seq_Deconv_Pipeline"
PROF_DIR="/hpcfs/users/a1667917/snakemake_slurm_profile"

cd $SNK_DIR

module load Anaconda3/2020.07
conda activate snakemake_clean_env

snakemake -s rna_seq_runner.smk --use-conda --config TCGA=True  Bams=/hpcfs/users/a1667917/Kevin/HNSC_RNA_BAMS/ \
Output=/hpcfs/users/a1667917/Kevin/RNA_TCGA_Out HG38_dir='/hpcfs/users/a1667917/STAR_Ref_Genomes' Salmon_dir='/hpcfs/users/a1667917/Salmon_Ref_Genomes' --profile $PROF_DIR/wgs_tcga

conda deactivate
