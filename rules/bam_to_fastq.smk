
#### if you want to only look at the unmapped READS
#### I have decided to look at them all for now

rule bam_sort:
    """Sort Unmapped Reads"""
    input:
        os.path.join(READS, "{sample}.bam")
    output:
        os.path.join(TMP,"{sample}_sorted.bam")
    log:
        os.path.join(LOGS,"{sample}_sort.log")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools sort -@ {threads} {input[0]} > {output[0]} 2> {log}
        """

rule bam_to_fastq:
    """converted unmapped reads to fastq"""
    input:
        os.path.join(READS, "{sample}.bam")
    output:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    log:
        os.path.join(LOGS,"{sample}.bam_to_fastq.log")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        samtools fastq -@ {threads} {input[0]} \
        -1 {output[0]} \
        -2 {output[1]} \
        -0 /dev/null -s /dev/null -n 2> {log}
        """

#### aggregation rule

rule test:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(TMP,"{sample}_R1.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "bam_to_fastq.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """
