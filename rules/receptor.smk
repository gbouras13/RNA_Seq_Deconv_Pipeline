rule trust4:
    input:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    params:
        os.path.join(TRUST4,"{sample}"),
        TRUST4DIR
    output:
        os.path.join(TRUST4, "{sample}_final.out")
    conda:
        os.path.join('..', 'envs','receptor.yaml')
    threads:
        16
    resources:
        mem_mb=MediumJobMem, 
        time=60
    shell:
        '''
        run-trust4 -f {params[1]}/human_IMGT+C.fa --ref {params[1]}/hg38_bcrtcr.fa -t {threads} \
         -1 {input[0]} -2 {input[1]} -o {params[0]}  
        '''


rule aggr_trust4:
    """aggregated"""
    input:
        expand(os.path.join(TRUST4,"{sample}", "{sample}_final.out"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_trust4.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=3
    shell:
        """
        touch {output[0]}
        """