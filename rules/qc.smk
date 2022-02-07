rule fastp_Trim:
    """remove adapters etc """
    input:
        os.path.join(TMP,"{sample}_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_R2.fastq.gz")
    output:
        os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_trim_R2.fastq.gz")
    log:
        os.path.join(LOGS,"{sample}_fastp.log")
    conda:
        os.path.join('..', 'envs','qc.yaml')
    threads:
        16
    resources:
        mem_mb=BigJobMem
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output[0]} -O {output[1]} -w {threads}
        """

rule fastqc:
    """fastqc trimmed reads"""
    input:
        expand(os.path.join(TMP,"{sample}_trim_R1.fastq.gz"), sample = SAMPLES),
        expand(os.path.join(TMP,"{sample}_trim_R2.fastq.gz"), sample = SAMPLES),
        TMP
    output:
        directory(FASTQC),
        os.path.join(FASTQC,"multiqc_report.html")
    log:
        os.path.join(LOGS,"fastqc.log")
    conda:
        os.path.join('..', 'envs','qc.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        fastqc -t {threads} -o {output[0]} {input[0]}
        fastqc -t {threads} -o {output[0]} {input[1]}
        multiqc {output[0]} {input[2]} -o {output[0]}
        """

#### aggregation rule

rule test_2:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(TMP,"{sample}_trim_R2.fastq.gz"), sample = SAMPLES)
    output:
        os.path.join(LOGS, "fastp.txt")
    threads:
        1
    resources:
        mem_mb=BigJobMem
    shell:
        """
        touch {output[0]}
        """
