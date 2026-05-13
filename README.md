# CLOverCNV
A custom tool for calling copy number variations in noisy genomic backgrounds

## CLOverCNV workflow:
1. Clean raw fastq files [fastp_clean.sh]
   > Name your files with the sample name as a prefix, IE: ABC1_sample_clean.fastq, Those names will automaticallly be interpreted in downstream steps. 
3. Align the reads using a custom BWA based aligner [bwa_align.sh]
   > This will make a folder that contains bam folders that contain different types of reads. This is will be called the 'bam folder' by other tools.
4. Use CLOverCNV.
   > This will create many files if you keep the working directories. It can be run one step at a time, or once alltogether (it will overide any single steps)

##
Example Usage for the three steps above:

Fastp: bash fastp_clean.sh --r1 {r1} --r2 {r2} --detect-adapters --cut-front --cut-tail --n-base-limit 20

BWA align: bash bwa_align.sh \
    --r1 {r1} --r2 {r2} --genome Dd2 --out-dir [scratch]/clover/clean/bwa_align_out --sample {sample} --out-key mem2 \
    --threads 8 --sort-threads 4  --sort-mem 256M --rg-id {sample} --rg-sm {sample} --rg-pl ILLUMINA \
    --audit-counts --use-samblaster --emit-discordant --emit-splitters --emit-unmapped --emit-supplementary --emit-primary

CLOverCNV:
bash CLOverCNV.sh \
  --sample {sample} \
  --genome Dd2.fasta \
  --gff ~/MalariAPI/genomes/Dd2.gff \
  --bam [scratch]/clover/cannon/bwa/{sample}.mem2.primary.sorted.bam \
  --bam-dir [scratch]/clover/cannon/bwa/ \
  --out-root [scratch]/clover/cannon/downsampled/CLOverCNV \
  --r-modules 'gcc/11.4.0 openmpi/4.1.4 R/4.3.1' \
  --train-window 75 \
  --train-drop-dup \
  --apply-binsize 75 \
  --lambda 30 --tile-bp 75 --min-tile-bp 75 \
  --min-probes-per-seg 5 \
  --flank 1000 \
  --agg sum  \
  --final-weak-ratio 1.5 \
  --final-strong-ratio 2.0 \
  --final-weak-z 0.3 \
  --final-strong-z 1.0 --final-count-flank 1000 \
  --confidence-method mean \
  --final-fuse --final-fuse-max-gap 3000 \
  --final-count-drop-dup \
  --final-keep-mode both \
  --circular-contigs 'PfDd2_API,PfDd2_MT' \
  --force
