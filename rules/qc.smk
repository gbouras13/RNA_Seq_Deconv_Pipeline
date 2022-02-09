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

rule bbmap:
    """ decontaminate sampe """
    input:
        os.path.join(TMP,"{sample}_trim_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_trim_R2.fastq.gz"),
        os.path.join(HG38_dir, 'GRCh38.primary_assembly.genome.fa')
    output:
        os.path.join(TMP,"{sample}_clean_R1.fastq.gz"),
        os.path.join(TMP,"{sample}_clean_R2.fastq.gz")
    log:
        os.path.join(LOGS,"{sample}_bbmap.log")
    conda:
        os.path.join('..', 'envs','qc.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        bbsplit.sh in1={input[0]} in2={input[1]} ref={input[2]} basename={wildcards.sample} outu1={output[0]}  outu2={output[1]} t={threads}
        """


rule fastqc:
    """fastqc trimmed reads"""
    input:
        fwd = expand(os.path.join(TMP,"{sample}_clean_R1.fastq.gz"), sample = SAMPLES),
        rev = expand(os.path.join(TMP,"{sample}_clean_R2.fastq.gz"), sample = SAMPLES),
        dir = TMP
    output:
        os.path.join(MULTIQC,"multiqc_report.html")
    params:
        fastqc = FASTQC,
        multiqc = MULTIQC
    conda:
        os.path.join('..', 'envs','qc.yaml')
    threads:
        BigJobCpu
    resources:
        mem_mb=BigJobMem
    shell:
        """
        fastqc -t {threads} -o {params.fastqc} {input.fwd}
        fastqc -t {threads} -o {params.fastqc} {input.rev}
        multiqc {params.fastqc} {input.dir} -o {params.multiqc}
        """





#### aggregation rule

rule test_2:
    """Index a .bam file for rapid access with samtools."""
    input:
        expand(os.path.join(TMP,"{sample}_clean_R2.fastq.gz"), sample = SAMPLES)
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
