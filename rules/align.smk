rule tcga_align_star:
    """align to hg38 tcga"""
    input:
        os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
    output:
        os.path.join(STAR_BAMS,"{sample}_star_50_"),
        os.path.join(STAR_BAMS,"{sample}_star_50_Aligned.sortedByCoord.out.bam")
    log:
        os.path.join(LOGS,"{sample}_star.log")
    params:
        os.path.join(HG38_dir, 'hg38_50')
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {params[0]} \
            --readFilesIn {input[0]} {input[1]} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {output[0]} \
            --outSAMtype BAM SortedByCoordinate
        """

rule enaalign_star:
    """align to hg38 ena """
    input:
        os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
    output:
        os.path.join(STAR_BAMS,"{sample}_star_75_"),
        os.path.join(STAR_BAMS,"{sample}_star_75_Aligned.sortedByCoord.out.bam")
    log:
        os.path.join(LOGS,"{sample}_star.log")
    params:
        os.path.join(HG38_dir, 'hg38_75')
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --genomeDir {params[0]} \
            --readFilesIn {input[0]} {input[1]} \
            --readFilesCommand gunzip -c \
            --outFileNamePrefix {output[0]} \
            --outSAMtype BAM SortedByCoordinate
        """


#### aggregation rule

rule aggr_align_tcga:
    input:
        expand(os.path.join(STAR_BAMS,"{sample}_star_50_Aligned.sortedByCoord.out.bam"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "star_50_align.txt")
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """

rule aggr_align_ena:
    input:
        expand(os.path.join(STAR_BAMS,"{sample}_star_75_Aligned.sortedByCoord.out.bam"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "star_75_align.txt")
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """




# rule feature_count:
#     """feature_counts """
#     input:
#         os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
#         os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
#     output:
#         os.path.join(TMP,"{sample}_star")
#     log:
#         os.path.join(LOGS,"{sample}_star.log")
#     params:
#         STAR_DIR
#     conda:
#         os.path.join('..', 'envs','align.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         STAR \
#             --runThreadN {threads} \
#             --genomeDir {params[0]} \
#             --readFilesIn {input[0]} {input[1]} \
#             --readFilesCommand gunzip -c \
#             --outFileNamePrefix {output[0]} \
#             --outSAMtype BAM SortedByCoordinate
#         """
# #
# #
#
# featureCounts -Q 10 -s 0 -T ${THREADS} -p -a $GTF -o $QUANTDATA/geneCounts.out ${sampleList}
#
#
#
# cut -f1,7- $QUANTDATA/geneCounts.out | sed 1d > $QUANTDATA/geneCounts.txt
