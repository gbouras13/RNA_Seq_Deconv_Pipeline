"""
All target output files are declared here
"""

# Preprocessing files
PreprocessingFiles = [
    os.path.join(LOGS, "aggr_bam_to_fastq.txt"),
    os.path.join(MULTIQC,"multiqc_report.html")
]


if TCGA:
    AlignFiles = os.path.join(LOGS, "star_100_align.txt")
    FeatureCountFiles = [
    os.path.join(RESULTS,"geneCounts_tcga.out"),
    os.path.join(RESULTS,"geneCounts_tcga.txt"),

    ]
else:
    AlignFiles = os.path.join(LOGS, "star_150_align.txt")
    FeatureCountFiles = [
        os.path.join(RESULTS,"geneCounts_ena.out"),
        os.path.join(RESULTS,"geneCounts_ena.txt")
    ]

SalmonFiles = os.path.join(LOGS, "salmon_agr.txt")
KallistoFiles = os.path.join(LOGS, "kallisto_agr.txt")
Trust4Files = os.path.join(LOGS, "aggr_trust4.txt")

# if TCGA:
#     SalmonFiles = os.path.join(LOGS, "salmon_agr_tcga.txt")
# else:
#     SalmonFiles = os.path.join(LOGS, "salmon_agr_ena.txt")


