
# from fastqs

# rule map_salmon:
#     """map salmon"""
#     input:
#         os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
#         os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
#     output:
#         os.path.join(SALMON_OUTPUT,"{sample}", "quant.sf")
#     log:
#         os.path.join(LOGS,"{sample}_salmon.log")
#     params:
#         os.path.join(Salmon_dir, 'gencode.v39_decoys_salmon'),
#         os.path.join(SALMON_OUTPUT, "{sample}")
#     conda:
#         os.path.join('..', 'envs','salmon.yaml')
#     threads:
#         BigJobCpu
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         "salmon quant -i {params[0]} -l A -1 {input[0]} -2 {input[1]} --validateMappings -p {threads} -o {params[1]}"

# rule aggr_salmon:
#     input:
#         expand(os.path.join(SALMON_OUTPUT,"{sample}", "quant.sf"), sample = SAMPLES)
#     output:
#         os.path.join(LOGS, "salmon_agr.txt")
#     threads:
#         1
#     resources:
#         mem_mb=BigJobMem
#     shell:
#         """
#         touch {output[0]}
#         """

# from bam

## ena
rule map_salmon_ena:
    """map salmon from star ena"""
    input:
         os.path.join(STAR_BAMS,"{sample}_star_150_Aligned.sortedByCoord.out.bam")
    output:
        os.path.join(SALMON_OUTPUT,"{sample}", "quant.sf")
    log:
        os.path.join(LOGS,"{sample}_salmon.log")
    params:
        os.path.join(Salmon_dir, 'gencode.v39_decoys_salmon'),
        os.path.join(SALMON_OUTPUT, "{sample}")
    conda:
        os.path.join('..', 'envs','salmon.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        "salmon quant -i {params[0]} -l A -a {input[0]} -p {threads} -o {params[1]}"

rule aggr_salmon_ena:
    input:
        expand(os.path.join(SALMON_OUTPUT,"{sample}", "quant.sf"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "salmon_agr_ena.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """


## ena
rule map_salmon_tcga:
    """map salmon from star ena"""
    input:
         os.path.join(STAR_BAMS,"{sample}_star_100_Aligned.sortedByCoord.out.bam")
    output:
        os.path.join(SALMON_OUTPUT,"{sample}", "quant.sf")
    log:
        os.path.join(LOGS,"{sample}_salmon.log")
    params:
        os.path.join(Salmon_dir, 'gencode.v39_decoys_salmon'),
        os.path.join(SALMON_OUTPUT, "{sample}")
    conda:
        os.path.join('..', 'envs','salmon.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        "salmon quant -i {params[0]} -l A -a {input[0]} -p {threads} -o {params[1]}"

rule aggr_salmon_tcga:
    input:
        expand(os.path.join(SALMON_OUTPUT,"{sample}", "quant.sf"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "salmon_agr_tcga.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """