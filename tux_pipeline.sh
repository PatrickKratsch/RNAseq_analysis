#!/usr/bin/env bash

## Configuration
# Programme executables
CPUS=16
HISAT2="/array/Patrick_RNA_seq/software/hisat2-2.1.0/hisat2"
SAMTOOLS="/array/Patrick_RNA_seq/software/samtools-1.6/samtools"
STRINGTIE="/array/Patrick_RNA_seq/software/stringtie-1.3.3b.Linux_x86_64/stringtie"
GENOME="/array/Patrick_RNA_seq/reference/ref_genome/bdgp6_tran/genome_tran"
GTF="/array/Patrick_RNA_seq/reference/gtf/Drosophila_melanogaster.BDGP6.91.gtf"

# In- and output directories
BASEDIR="/array/Patrick_RNA_seq/data/GEPD_run/tuxedo"
FASTQ="$BASEDIR/FASTQ"
mkdir $BASEDIR/aligned
mkdir $BASEDIR/BAM
mkdir $BASEDIR/assembled

## Main
# 1. Align - max-intronlen is set to the 99th centile of Dmel intron lengths
for f in $(cat fastq_files); do
  $HISAT2 -p $CPUS --dta -x $GENOME -1 $FASTQ/$f/$f\_R1_001.fastq \
  -2 $FASTQ/$f/$f\_R2_001.fastq -S $BASEDIR/aligned/$f.sam \
  --rna-strandness RF --max-intronlen 26745
done

# 2. Convert to BAM
cd $BASEDIR/aligned
for s in $(echo *.sam); do
  $SAMTOOLS sort -@ $CPUS -o $BASEDIR/BAM/$s.bam \
  $BASEDIR/aligned/$s
done

# 3. Index BAM
cd $BASEDIR/BAM
# for b in $(echo *.bam); do
  # $SAMTOOLS index -@ $CPUS $BASEDIR/BAM/$b
# done

# 4. Assemble and quantify
for b in $(echo *.bam); do
  $STRINGTIE -p $CPUS -G $GTF -o $BASEDIR/assembled/$b.gtf \
  $BASEDIR/BAM/$b
done








# cd /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/assembled
# touch mergelist.txt
# for gtf_dir in $(ls); do
# mkdir /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/merged_gtf
# /array/Patrick_RNA_seq/software/stringtie-1.3.3b.Linux_x86_64/stringtie \
# -p 4 \
# -G /array/Patrick_RNA_seq/reference/gtf/Drosophila_melanogaster.BDGP6.91.gtf \
# -o /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/merged_gtf/merged.gtf \
# /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/assembled/mergelist.txt

# cd /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/BAM
# mkdir /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/ballgown
# for b in $(echo *.bam); do
  # mkdir /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/ballgown/$b
  # /array/Patrick_RNA_seq/software/stringtie-1.3.3b.Linux_x86_64/stringtie \
  # -e -B -p 8 \
  # -G /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/merged_gtf/merged.gtf \
  # -o /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/ballgown/$b/$b\_merged.gtf \
  # /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/BAM/$b
# done

exit
