#!/usr/bin/env bash

mkdir /array/Patrick_RNA_seq/data/GEPD_run/bioconductor/aligned

files=`cat fastq_files`

for f in $files; do
  /array/Patrick_RNA_seq/software/STAR-2.5.3a/bin/Linux_x86_64/STAR \
  --runThreadN 4 \
  --alignIntronMax 50000 \
  --genomeDir /array/Patrick_RNA_seq/reference/ref_genome/STAR/STAR_index \
  --readFilesIn /array/Patrick_RNA_seq/data/GEPD_run/bioconductor/FASTQ/$f/$f\_R1_001.fastq \
  /array/Patrick_RNA_seq/data/GEPD_run/bioconductor/FASTQ/$f/$f\_R2_001.fastq \
  --outFileNamePrefix /array/Patrick_RNA_seq/data/GEPD_run/bioconductor/aligned/$f.
done

mkdir /array/Patrick_RNA_seq/data/GEPD_run/bioconductor/BAM

cd /array/Patrick_RNA_seq/data/GEPD_run/bioconductor/aligned

sams=`echo *.sam`

for s in $sams; do
  /array/Patrick_RNA_seq/software/samtools-1.6/samtools \
  sort -@ 4 \
  -o /array/Patrick_RNA_seq/data/GEPD_run/bioconductor/BAM/$s.bam \
  /array/Patrick_RNA_seq/data/GEPD_run/bioconductor/aligned/$s
done

cd /array/Patrick_RNA_seq/data/GEPD_run/bioconductor/BAM

bams=`echo *.bam`

for b in $bams; do
  /array/Patrick_RNA_seq/software/samtools-1.6/samtools \
  index -@ 4 \
  /array/Patrick_RNA_seq/data/GEPD_run/bioconductor/BAM/$b
done

exit
