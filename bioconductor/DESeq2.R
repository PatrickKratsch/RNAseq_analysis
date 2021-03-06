## DESeq2 RNA-seq analysis of GD vs. lP (9 replicates)
## per genotype) - Alignment: STAR 2.5.3a - Pipeline from:
## https://www.bioconductor.org/help/workflows/rnaseqGene/
## paper: https://f1000research.com/articles/4-1070

## Before DESeq2

# Align reads with STAR and convert SAM to BAM with samtools

# Make a sample table in Excel and load
# Include all metadata for the samples you deem necessary
# Use row.names = 1 to indicate that row names are contained
# in first column
sampleTable <- read.csv("sampleTable.csv", row.names = 1)

# Update bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite()

## Start DESeq2 analysis
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")
library("BiocParallel")
library("DESeq2")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("ggbeeswarm")
library("genefilter")
library("AnnotationDbi")
library("org.Dm.eg.db")
library("ReportingTools")
library("data.table")
# Construct the full paths to the BAM files we want to perform counting
# operations on
wd <- getwd()
filenames <- file.path(wd, "BAM", paste0(rownames(sampleTable), ".Aligned.out.sam.bam"))
file.exists(filenames)

# Load BAM files
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])

# Load GTF file and indicate that none of the sequences are circular
gtffile <- file.path("/Volumes/array/Patrick_RNA_seq/reference/gtf", 
                     "Drosophila_melanogaster.BDGP6.91.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())

# Generatea GRangesList of all exons grouped by gene - each element
# of the list is a GRanges object of the exons for a gene
ebg <- exonsBy(txdb, by="gene")

# Read counting - specify to use only one core
register(SerialParam())
# Create a SummarizedExperiment object with counts
# The reason to use invertStrand can be found here:
# https://support.bioconductor.org/p/65844/
# Note that this gives a warning that some sequence levels
# are not present in either x or y --> problem with building
# genome index?
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=FALSE,
                        fragments=TRUE,
                        preprocess.reads=invertStrand)

# Add colData to se
colData(se) <- DataFrame(sampleTable)

# Check whether the first factor of se$genotype is the control
se$genotype

# Make DESeqDataSet with appropriate analysis design
dds <- DESeqDataSet(se, design = ~ batch + genotype)

# Perform filtering to reduce object size of dds
# Not that this means that we disregard any genes for which
# the number of counts across all samples is 1 or lower
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Perform rlog transformation
# rlog tends to work well on small datasets (n < 30)
rld <- rlog(dds, blind = FALSE)

# Perform VST transformation
vsd <- vst(dds, blind = FALSE)

# Visualise the effect of the rlog transformation
# NB: These transformations are used for exploratory
# analysis, but not for differential gene expression
# analysis, for which raw counts are used
dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation) 

# Save this plot with 300 dpi resolution (instead of default 72)
ggsave("transformations.png")

# Assess the similarity between samples by calculating the 
# Euclidian distances between samples - for this, use the 
# rlog-transformed counts, so that all genes contribute equally
# Plot1: Visualise similarity with a heatmap
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$genotype, rld$batch, sep = " - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

# Plot2: Alternatively to the Euclidian distance, we can plot the
# Poisson distance - both these visualisations show the differences
# in counts between all samples
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste(rld$genotype, rld$batch, sep=" - ")
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# Plot3: PCA plot - another way to visualise sample-to-sample distances
plotPCA(rld, intgroup = c("genotype", "batch"))

# Plot4: PCA plot with ggplot2 - extract PCA data from plotPCA and feed
# into ggplot2
pcaData <- plotPCA(rld, intgroup = c("genotype", "batch"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = genotype, shape = batch)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed()

# Save PCA plot with 300 dpi
ggsave("./R_plots/PCA_libraries.png")

# Plot5: MDS plot - multidimensional scaling = similar to PCA 
mds <- as.data.frame(colData(rld))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = genotype, shape = batch)) +
  geom_point(size = 3) + coord_fixed()

# Plot6: MDS plot with PoissonDistance
mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = genotype, shape = batch)) +
  geom_point(size = 3) + coord_fixed()

## Differential expression analysis - we have already
# specified the design above, so can run differential expression
# on raw counts with a single call to DESeq()
dds <- DESeq(dds)

# Build a results table
# Note that by default, the FDR cutoff of results() is 0.1
res <- results(dds)
# Get metadata of res
mcols(res, use.names = TRUE)
# Get more summary of res
summary(res)
# Filter results (res) by setting a lower cutoff for the
# FRD (padj) - set it to < 0.05
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
# Filter res by setting a cutoff for the log2(FC)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

# Note that the padj values are BH-adjusted p values - it is the
# answer to the question: If one were to call all genes that have
# a specific adjusted p value significant, what would be the 
# fraction of false positives in that list? By setting the padj to
# < 0.05, we are accepting that 5 % of the genes that show a significant
# difference in expression based on this padj value will be false positives
# How many genes with a padj < 0.05 are there?
sum(res$padj < 0.05, na.rm=TRUE)

# We now subset the results table to genes with a padj < 0.05, and then sort 
# the results by the log2(FC) estimate to get the significant genes with the
# strongest down-regulation
resSig <- subset(res, padj < 0.05)
head(resSig[order(resSig$padj), ], 9)

# Plot7: Plot counts of specific genes between groups
# Do this for pdf here - use any FB gene id
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = "FBgn0014018", intgroup=c("genotype"))

# Plot8: Repeat Plot7 as a ggplot version
geneCounts <- plotCounts(dds, gene = "FBgn0014018", intgroup = c("genotype", "batch"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = genotype, y = count, color = batch)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

# Plot9: MA-plot - on the y-axis, M stands for minus, on the x-axis,
# A stands for average - note that this plot shows only genes with
# a padj threshold of < 0.1 (default)
res <- lfcShrink(dds, contrast=c("genotype","m","c"), res=res)
plotMA(res, alpha=0.05, ylim = c(-1.2, 1))

# Plot10: If we don't shrink the fold changes (which mainly affects
# genes with very low counts and highly variable counts), the plot 
# looks very different:
res.noshr <- results(dds)
plotMA(res.noshr, ylim = c(-9, 9))

# Plot11: Label specific genes on the MA plot
plotMA(res, alpha=0.05, ylim = c(-1.2, 1))
topGene <- rownames(res)[which.min(res$padj)]
with(res["FBgn0039298", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0039298", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0010041", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0010041", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0028374", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0028374", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0243514", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0243514", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0013812", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0013812", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0013813", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0013813", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0035160", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0035160", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0283437", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0283437", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0033367", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0033367", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0033926", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0033926", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0033928", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0033928", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0035766", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0035766", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0027348", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0027348", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0267792", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0267792", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0014018", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0014018", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})
with(res["FBgn0042201", ], {
  gene_name_dot <- mapIds(org.Dm.eg.db, keys="FBgn0042201", column="SYMBOL", 
                          keytype="FLYBASE", multiVals="first")
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, gene_name_dot, pos=2, col="dodgerblue")
})

# Plot12: Plot a histogram of p values - exclude very lowly expressed
# genes
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

# Plot13: Cluster the 20 genes with highest variance across samples,
# using the rlog-transformed counts
# We then draw a heatmap - this will be done not on absolute expression,
# but rather on the amount by which a specific gene deviates in
# a specific sample from the gene's average across all samples
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)
mat  <- assay(rld)[topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("batch","genotype")])
pheatmap(mat, annotation_col = anno)

# Show how the fractions of pvalues < 0.05 occur over genes
# of differing counts
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
  mean(p < .05, na.rm = TRUE))
barplot(fractionSig, xlab = "mean normalized count",
        ylab = "fraction of small p values")

## Annotating and exporting results
# Annotate FB gene IDs
res$symbol <- mapIds(org.Dm.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="FLYBASE",
                     multiVals="first")

res$FBid <- mapIds(org.Dm.eg.db,
                   keys=row.names(res),
                   column="FLYBASE",
                   keytype="FLYBASE",
                   multiVals="first")

resOrdered <- res[order(res$pvalue), ]

# Exporting results
resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file = "genes.csv")

htmlRep <- HTMLReport(shortName="report", title="report",
                      reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)

# Filter results and only output padj < 0.05
sign <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]
signOrdered <- sign[order(sign$padj), ]
signOrderedDF <- as.data.frame(signOrdered)
write.csv(signOrderedDF, "sig_genes.csv")

## Extract the intersect significant genes of bc and tux pipelines
# Use awk and bash to generate a text file contining the names
# of the genes that are significant in both pipelines
intersect_genes <- read.table("intersect_genes.txt")
signOrderedDT <- setDT(signOrderedDF)
intersectDT <- signOrderedDT[symbol %in% intersect_genes$V1, ]
intersectDTordered <- intersectDT[order(intersectDT$padj), ]
write.csv(intersectDTordered, "intersect_genes.csv")

## New heatmap with the 50 most significant genes (based on FDR)
# Note that we need the loop to get the correct gene symbols
# that correspond to the gene IDs
# This heatmap is based on rld-transformed raw counts - we are
# taking the mean of the rld-transformed counts across all samples,
# and then subtracting that mean from the rld-transformed counts
# for each sample - that way, we get the relative expression of
# each gene per sample relative to the mean of the rld-transformed counts
# for each sample
most_sign_genes <- as.data.frame(cbind(as.character(head(resOrderedDF$FBid, 50)),
                         as.character(head(resOrderedDF$symbol, 50))))
most_sign_genes_ind <- which(row.names(assay(rld)) %in% most_sign_genes$V1)
mat  <- assay(rld)[most_sign_genes_ind, ]
mat  <- mat - rowMeans(mat)

mat_row_names <- row.names(assay(rld))[most_sign_genes_ind]
gene_names_ordered <- list()
i <- 1
for(id in mat_row_names){
  j <- 1
  for(symbol in most_sign_genes$V1){
    if(id == symbol){
      gene_names_ordered[i] <- as.character(most_sign_genes$V2[j])
      i <- i + 1
    }
    j <- j + 1
  }
}
gene_names_ordered <- sapply(gene_names_ordered, paste0, collapse="")

row.names(mat) <- gene_names_ordered
anno <- as.data.frame(colData(rld)[, c("batch","genotype")])
pheatmap(mat, annotation_col = anno)

## New heatmap with selected genes based on DAVID gene ontology
# analysis - focus on 4 different clusters here - use awk to
# extract the Flybase IDs of these genes, save as .txt file, and
# lead into R - Do one heatmap per functional cluster, because otherwise
# the trends are not easily seen
DAVID_cluster <- read.table("gst.txt")
DAVID_ind <- which(resOrderedDF$FBid %in% DAVID_cluster$V1)
DAVID_genes <- as.data.frame(cbind(as.character(resOrderedDF$FBid[DAVID_ind]),
                                       as.character(resOrderedDF$symbol[DAVID_ind])))
DAVID_rld_ind <- which(row.names(assay(rld)) %in% DAVID_genes$V1)
mat  <- assay(rld)[DAVID_rld_ind, ]
mat  <- mat - rowMeans(mat)

mat_row_names <- row.names(assay(rld))[DAVID_rld_ind]
gene_names_ordered <- list()
i <- 1
for(id in mat_row_names){
  j <- 1
  for(symbol in DAVID_genes$V1){
    if(id == symbol){
      gene_names_ordered[i] <- as.character(DAVID_genes$V2[j])
      i <- i + 1
    }
    j <- j + 1
  }
}
gene_names_ordered <- sapply(gene_names_ordered, paste0, collapse="")

row.names(mat) <- gene_names_ordered
anno <- as.data.frame(colData(rld)[, c("batch","genotype")])
pheatmap(mat, annotation_col = anno, show_rownames = T)

