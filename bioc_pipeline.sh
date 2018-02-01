#!/usr/bin/env bash

## Configuration
# Programme executables
CPUS=16
STAR="/array/Patrick_RNA_seq/software/STAR-2.5.3a/bin/Linux_x86_64/STAR"
SAMTOOLS="/array/Patrick_RNA_seq/software/samtools-1.6/samtools"
GENOME="/array/Patrick_RNA_seq/reference/ref_genome/STAR/STAR_index"

# In- and output directories
BASEDIR="/array/Patrick_RNA_seq/data/GEPD_run/bioconductor"
FASTQ="$BASEDIR/FASTQ"
mkdir $BASEDIR/aligned
mkdir $BASEDIR/BAM


## Main
# 1. Align - max-intronlen is set to the 99th centile of Dmel intron lengths
for f in $(cat fastq_files); do
  $STAR --runThreadN $CPUS --alignIntronMax 26745 --genomeDir $GENOME \
  --readFilesIn $FASTQ/$f/$f\_R1_001.fastq \
  $FASTQ/$f/$f\_R2_001.fastq --outFileNamePrefix $BASEDIR/aligned/$f.
done

cd $BASEDIR/aligned
for s in $(echo *.sam); do
  $SAMTOOLS sort -@ $CPUS -o $BASEDIR/BAM/$s.bam $BASEDIR/aligned/$s
done

cd $BASEDIR/BAM
for b in $(echo *.bam); do
  $SAMTOOLS index -@ $CPUS $BASEDIR/BAM/$b
done

exit
