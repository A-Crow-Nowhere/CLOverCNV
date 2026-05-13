# CLOverCNV
A custom tool for calling copy number variations in noisy genomic backgrounds

## CLOverCNV workflow:
1. Clean raw fastq files [fastp_clean.sh]
2. Align the reads using a custom BWA based aligner [bwa_align.sh]
