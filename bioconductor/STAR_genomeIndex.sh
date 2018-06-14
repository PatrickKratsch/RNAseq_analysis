#!/usr/bin/env bash

## Configuration
CPUS=24
STAR="/array/Patrick_RNA_seq/software/STAR-2.5.3a/bin/Linux_x86_64/STAR"
BASEDIR="/array/Patrick_RNA_seq/reference/ref_genome/STAR"
GTF="/array/Patrick_RNA_seq/reference/gtf/Drosophila_melanogaster.BDGP6.91.gtf"
mkdir $BASEDIR/STAR_index

## Download Dmel genomic fasta and concatenate
wget ftp://ftp.ensembl.org/pub/release-91/fasta/drosophila_melanogaster/dna/*.dna.chromosome*
wget ftp://ftp.ensembl.org/pub/release-91/fasta/drosophila_melanogaster/dna/*.dna.nonchromosomal*
gunzip *.gz

## Generate genome index
# Set --genomeSAindexNbases to 13, due to BDGEP6 genome length = 142,573,017
# --> min(14, log2(142,573,017)/2 - 1) =~ 12.5 =~ 13
$STAR --runThreadN $CPUS --runMode genomeGenerate \
--genomeDir $BASEDIR/STAR_index \
--genomeFastaFiles $BASEDIR/*.fa \
--sjdbGTFfile $GTF \
--sjdbOverhang 75 \
--genomeSAindexNbases 13

exit
