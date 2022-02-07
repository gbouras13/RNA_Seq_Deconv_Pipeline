rule align_star:
    """align to hg38 """
    input:
        os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
    output:
        os.path.join(RESULTS,"{sample}_star_"),
        os.path.join(RESULTS,"{sample}_star_Aligned.sortedByCoord.out.bam")
    log:
        os.path.join(LOGS,"{sample}_star.log")
    params:
        os.path.join(STAR_DIR, 'hg38_50')
    conda:
        os.path.join('..', 'envs','align.yaml')
    threads:
        2
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

rule aggr_align:
    input:
        expand(os.path.join(RESULTS,"{sample}_star_Aligned.sortedByCoord.out.bam"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "star_align.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """


#
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
