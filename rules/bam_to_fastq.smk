
def get_bam(wildcards):
    return sampleBams[wildcards.sample]["bam"]

#### if you want to only look at the mapped READS

rule bam_map_sort_fastq:
    """converted mapped reads to fastq"""
    input:
        get_bam,
    output:
        os.path.join(UNALIGNED_FASTQ,"{sample}_R1.fastq.gz"),
        os.path.join(UNALIGNED_FASTQ,"{sample}_R2.fastq.gz")
    conda:
        os.path.join('..', 'envs','samtools.yaml')
    threads:
        32
    resources:
        mem_mb=32000,
        time=180
    wildcard_constraints:
        sample="[^/]+"
    shell:
        """
        samtools view -u -F 4 -@ {threads} {input[0]} | samtools sort -n -@ {threads} |   
        samtools fastq -@ {threads} \
        -1 {output[0]} \
        -2 {output[1]} \
        -0 /dev/null -s /dev/null -n 
        """



#### aggregation rule

rule aggr_bam_to_fastq:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(UNALIGNED_FASTQ,"{sample}_R1.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "aggr_bam_to_fastq.txt")
    threads:
        1
    resources:
        mem_mb=SmallJobMem,
        time=5
    wildcard_constraints:
        sample="[^/]+"
    shell:
        """
        touch {output[0]}
        """
