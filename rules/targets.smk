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





# else:
#     # Assembly files
#     AssemblyFiles = [
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","assembly.fasta"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","MAPPING","contig_count_table.tsv"),
#         os.path.join(RESULTS,"assembly.properties.tsv")]
#     # Contig annotations
#     ContigAnnotFiles = [
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_phylum_summary.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_class_summary.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_order_summary.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_family_summary.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_genus_summary.tsv"),
#         os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_species_summary.tsv")
#     ]
#     # Mapping files
#     MappingFiles = [
#         os.path.join(MAPPING,"assembly.seqtable.bam"),
#         os.path.join(MAPPING,"assembly.seqtable.bam.bai"),
#         os.path.join(RESULTS,"contigSeqTable.tsv"),
#         os.path.join(SUMDIR,"contigKrona.html")
#     ]
#
#
# # Secondary AA search files
# SecondarySearchFilesAA = [
#     os.path.join(SECONDARY_AA_OUT, "AA_bigtable.tsv"),
# ]
#
# # Secondary NT search files
# SecondarySearchFilesNT = [
#     os.path.join(SECONDARY_NT_OUT, "NT_bigtable.tsv"),
#     os.path.join(RESULTS, "bigtable.tsv"),
# ]
#
# # Summary files
# optionalSummary = [
#     os.path.join(SUMDIR,"Step00_counts.tsv"),
#     os.path.join(SUMDIR,"Step01_counts.tsv"),
#     os.path.join(SUMDIR,"Step02_counts.tsv"),
#     os.path.join(SUMDIR,"Step03_counts.tsv"),
#     os.path.join(SUMDIR,"Step04_counts.tsv"),
#     os.path.join(SUMDIR,"Step05_counts.tsv"),
#     os.path.join(SUMDIR,"Step06_counts.tsv"),
#     os.path.join(SUMDIR,"Step07_counts.tsv"),
#     os.path.join(SUMDIR,"Step08_counts.tsv"),
#     os.path.join(SUMDIR,"Step09_counts.tsv"),
#     os.path.join(SUMDIR,"Step10_counts.tsv"),
#     os.path.join(SUMDIR,"Step11_counts.tsv"),
#     os.path.join(SUMDIR,"Step12_counts.tsv"),
#     os.path.join(SUMDIR,"Step13_counts.tsv"),
#     os.path.join(SUMDIR,"Sankey.svg"),
# ]
# SummaryFiles = [
#     optionalSummary,
#     os.path.join(SUMDIR, 'hecatomb.samples.tsv'),
#     os.path.join(SUMDIR, "taxonLevelCounts.tsv"),
#     os.path.join(SUMDIR, "krona.html")
# ]
