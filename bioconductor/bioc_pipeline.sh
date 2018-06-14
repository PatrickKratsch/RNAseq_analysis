#!/usr/bin/env bash

## Configuration
# Programme executables
CPUS=24
STAR="/array/Patrick_RNA_seq/software/STAR-2.5.3a/bin/Linux_x86_64/STAR"
SAMTOOLS="/array/Patrick_RNA_seq/software/samtools-1.6/samtools"
GENOME="/array/Patrick_RNA_seq/reference/ref_genome/STAR/STAR_index"

# In- and output directories
BASEDIR="/array/Patrick_RNA_seq/data/GEPD_run/bioconductor"
FASTQ="$BASEDIR/FASTQ"
mkdir $BASEDIR/aligned
mkdir $BASEDIR/BAM

################################################################################
#### Calculate the 99th centile of Drosophila melanogaster before alignment ####
################################################################################

# INTRONFILE="dmel-all-intron-r6.21.fasta"
# INTRONS="ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.21_FB2018_02/fasta/$INTRONFILE.gz"
# wget/curl $INTRONS # depending on Linux vs. Mac, respectively
# gunzip $INTRONFILE.gz
# grep "^>" $INTRONFILE | cut -f7 -d ";" | cut -f2 -d "=" > introns

################################################################################
######## Load 'introns' into R and calculate the 99th centle = 26745 bp ########
################################################################################
MAXINTRON="26745"

## Main
# 1. Align - max-intronlen is set to the 99th centile of Dmel intron lengths
for f in $(cat fastq_files); do
  $STAR --runThreadN $CPUS --alignIntronMax $MAXINTRON --genomeDir $GENOME \
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
