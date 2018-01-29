#!/usr/bin/env bash

mkdir /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/aligned
files=`cat fastq_files`
for f in $files; do
  /array/Patrick_RNA_seq/software/hisat2-2.1.0/hisat2 \
  -p 8 --dta \
  -x /array/Patrick_RNA_seq/reference/ref_genome/bdgp6_tran/genome_tran \
  -1 /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/FASTQ/$f/$f\_R1_001.fastq \
  -2 /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/FASTQ/$f/$f\_R2_001.fastq \
  -S /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/aligned/$f.sam \
  --rna-strandness RF --max-intronlen 26745
done

mkdir /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/BAM
cd /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/aligned
sams=`echo *.sam`
for s in $sams; do
  /array/Patrick_RNA_seq/software/samtools-1.6/samtools \
  sort -@ 4 \
  -o /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/BAM/$s.bam \
  /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/aligned/$s
done

cd /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/BAM
bams=`echo *.bam`
for b in $bams; do
  /array/Patrick_RNA_seq/software/samtools-1.6/samtools \
  index -@ 4 \
  /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/BAM/$b
done

mkdir /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/assembled
for b in $bams; do
  /array/Patrick_RNA_seq/software/stringtie-1.3.3b.Linux_x86_64/stringtie \
  -p 8 \
  -G /array/Patrick_RNA_seq/reference/gtf/Drosophila_melanogaster.BDGP6.91.gtf \
  -o /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/assembled/$b.gtf \
  -l $b\_assembled \
  /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/BAM/$b
done

cd /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/assembled
gtf_dirs=ls
touch mergelist.txt
for gtf_dir in $gtf_dirs; do
  cd gtf_dir
  printf '%s\n' *.gtf > ../mergelist.txt
  cd ..
done
mkdir /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/merged_gtf
/array/Patrick_RNA_seq/software/stringtie-1.3.3b.Linux_x86_64/stringtie \
-p 4 \
-G /array/Patrick_RNA_seq/reference/gtf/Drosophila_melanogaster.BDGP6.91.gtf \
-o /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/merged_gtf/merged.gtf \
/array/Patrick_RNA_seq/data/GEPD_run/tuxedo/assembled/mergelist.txt

cd /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/BAM
mkdir /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/ballgown
for b in $bams; do
  mkdir /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/ballgown/$b
  /array/Patrick_RNA_seq/software/stringtie-1.3.3b.Linux_x86_64/stringtie \
  -e -B -p 8 \
  -G /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/merged_gtf/merged.gtf \
  -o /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/ballgown/$b/$b\_merged.gtf \
  /array/Patrick_RNA_seq/data/GEPD_run/tuxedo/BAM/$b
done

exit
