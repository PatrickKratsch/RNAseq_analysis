### RNA-seq analysis
### Analysis based on https://rpubs.com/kapeelc12/Ballgown
library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)
library(data.table)
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(plyr)
library(plotrix)
library(RColorBrewer)

## Reading in ballgown object
# Need to generate pheno_data manually, e.g. in Excel
# ballgown() reads in the Stringtie output (exon-, intron-, transcript data)
# The ballgown object stores the expression data for each feature
pheno_data <- read.csv("pheno_data.csv", header = TRUE)
bg <- ballgown(dataDir = "ballgown", samplePattern = "sam.bam", 
               pData = pheno_data)

## Filter ballgown object
# Only return transcripts with a mean expression across samples above 1 FPKM
bg_filt <- exprfilter(bg, 1, meas = "FPKM")
# Remove all transcripts with a variance across samples of less than 1
bg_filt <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset = TRUE)

## Load gene names for lookup later
# columns 9 and 10 from the transcript data.frame contain gene_id and 
# gene_name, respectively
bg_gene_names <- unique(texpr(bg_filt, 'all')[, 9:10])

## Pull gene expression data.frame from ballgown object and convert matrix
# to data.frame
gene_expression <- as.data.frame(gexpr(bg_filt))

## Change column names
colnames(gene_expression) <- c("17_G1", "17_G2", "17_G5",
                               "17_l1", "17_l3", "17_l6",
                               "18_G1", "18_G2", "18_G3",
                               "18_G4", "18_G5", "18_G6",
                               "18_l1", "18_l2", "18_l3",
                               "18_l4", "18_l5", "18_l6")

## Assign colors to each sample (randomly)
color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
data_colors <- sample(color, 18)

## Load transcript to gene index from ballgown object
transcript_gene_table <- indexes(bg)$t2g

## Plot1: Distribution of transcripts/gene as frequency plot
counts <- table(transcript_gene_table[, "g_id"])
c_one <- length(which(counts == 1))
c_more_than_one <- length(which(counts > 1))
c_max <- max(counts)
hist(counts, breaks = 100, col = "bisque4", xlab = "Transcripts per gene", 
     main = "Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), 
                paste("Genes with more than one transcript =", 
                c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty=NULL)

## Plot2: Distribution of transcript sizes
full_table <- texpr(bg , 'all')
hist(full_table$length, breaks = 100, xlab = "Transcript length (bp)", 
     main="Distribution of transcript lengths", col = "steelblue")

## Min and max FPKM values for each library
i <- 1
maxi <- list()
for(colname in colnames(gene_expression)){
  maxi[i] <- max(gene_expression[, colname])
  i <- i + 1
}

j <- 1
mini <- list()
for(colname in colnames(gene_expression)){
  mini[j] <- min(gene_expression[, colname])
  j <- j + 1
}

## Plot3: Plot maximum gene expression per sample/library (by Patrick Kratschmer)
maxi_df <- as.data.frame(cbind(sample=as.vector(colnames(gene_expression)), 
                               max_expr=as.numeric(maxi), 
                               genotype=c(rep("GD",3),
                                          rep("lP",3),
                                          rep("GD",6),
                                          rep("lP",6))))
ggplot(maxi_df, aes(x=sample, y=max_expr, colour=genotype)) +
  geom_point(size=2) +
  theme_bw()

## Plot4: FPKM distribution across libraries
labels <- colnames(gene_expression)
boxplot(log2(gene_expression+1), col = "steelblue", 
        xaxt = "n", xlab = "", ylab = "log2(FPKM)", 
        main = "FPKM Distribution")
staxlab(1, 1:18, labels, nlines = 3)

## Plot5: Plot pairs of replicates to assess reproducibility of technical 
## and biological replicates
x = gene_expression[,"18_GD1"]
y = gene_expression[,"17_GD1"]
plot(x = log2(x + 1), y=log2(y + 1), 
     pch = 16, col = "blue", cex = 0.25, xlab = "FPKM (18_GD1)", 
     ylab = "FPKM (17_GD1)", 
     main = "Comparison of Expression Values for a Pair of Replicates")
abline(a = 0, b = 1)
rs = cor(x, y)^2
legend("topleft", paste("R squared = ", 
                        round(rs, digits=3), sep=""), lwd=1, col="black")
legend("bottomright", paste("Number of genes = 6654"))

## Plot6: Regenerate Plot4 as density scatter plot
colors = colorRampPalette(c("white", "blue", "#007FFF", 
                            "cyan","#7FFF7F", "yellow", 
                            "#FF7F00", "red", "#7F0000"))
smoothScatter(x=log2(x + 1), y=log2(y + 1), 
              xlab="FPKM (18_GD1)", 
              ylab="FPKM (17_GD1)", 
              main="Comparison of expression values for a pair of replicates", 
              colramp=colors, nbin=200)

## Compare correlation distance of all libaries to one another
# Calculate FPKM sum across libraries, and select overall FPKM > 5
gene_expression[, "sum"] = apply(gene_expression, 1, sum)
i = which(gene_expression[, "sum"] > 5)
r = cor(gene_expression[i,], 
      use = "pairwise.complete.obs", method = "pearson")

## Plot7: Convert correlation to physical distance, and use
## multi-dimentional scaling to visualise this; libraries with
## similar expression patterns should group together in space
data_colours <- c(rep("tomato1", 3), rep("wheat1", 3),
                  rep("tomato1", 6), rep("wheat1", 6))
d = 1 - r
mds = cmdscale(d, k = 2, eig = TRUE)
par(mfrow = c(1, 1))
plot(mds$points, type = "n", xlab = "", ylab = "", 
     main = "MDS Distance Plot for all Libraries", 
     xlim = c(-0.02, 0.02), ylim = c(-0.02, 0.02))
points(mds$points[,1], mds$points[,2], col = "grey", cex =  2, ph = 16)
text(mds$points[, 1], mds$points[, 2], colnames(gene_expression),
     col = data_colours)


## Use standard linear model comparison for transcript expression 
results_transcripts <- stattest(bg_filt, feature = "transcript",   
                                covariate = "genotype", 
                                adjustvars = c("batch"), 
                                getFC = TRUE, meas = "FPKM")

## Do the same for gene expression 
results_genes <- stattest(bg_filt, feature = "gene", 
                          covariate = "genotype", 
                          adjustvars = c("batch"), 
                          getFC = TRUE, meas = "FPKM")

## Plot8: View distribution of differential expression values as histogram
## Display only significant genes according to Ballgown
sig = which(results_genes$qval < 0.05)
results_genes[, "log2(FC)"] = log2(results_genes[, "fc"])
hist(results_genes[sig, "log2(FC)"], breaks = 50, col = "steelblue", 
     xlab = "log2(FC) GD vs lP", 
     main = "Distribution of Differential Expression Values")
abline(v = -1, col = "black", lwd = 2, lty = 2)
abline(v = 1, col = "black", lwd = 2, lty = 2)
legend("topleft", "Fold-change > 2", lwd = 2, lty = 2)
legend("topright", "Number of DEG = 772")

## Plot9: Display all expression values and mark the ones that are significant
# Note that these data are based on bg_filt (because gene_expression
# was made by fetching gexpr(bg_filt), hence all genes present here are
# pre-filtered!)
gene_expression[, "GD"] = apply(gene_expression[, c(1:3, 7:12)], 1, mean)
gene_expression[, "lP"] = apply(gene_expression[, c(4:6, 13:18)], 1, mean)
x = log2(gene_expression[, "GD"] + 1)
y = log2(gene_expression[, "lP"] + 1)
plot(x = x, y = y, pch = 16, cex = 0.25, xlab = "GD log2(FPKM)", 
     ylab="lP log2(FPKM)", main="GD vs lP")
abline(a = 0, b = 1)
xsig = x[sig]
ysig = y[sig]
points(x = xsig, y = ysig, col = "red", pch = 16, cex = 0.5)
legend("topleft", "Significant", col = "red", pch = 16)
legend("bottomright", "Number of DEG = 722\nNumber of genes = 6654")

## Plot10: Vulcano plot (from Mingyang and Alan)
# Add a threshold column reporting the labels for qvals
vulcano_genes <- results_genes %>% mutate(threshold=ifelse(qval < 0.01, "A", 
                                                           ifelse(qval < 0.05 & qval >= 0.01, "B", "C")))
## Construct the plot object
ggplot(vulcano_genes, aes(x=log2(fc), y=-log10(qval), colour=threshold)) +
  geom_point(shape=21, size=1.25, aes(fill=factor(threshold))) +
  xlim(c(-3.5, 3.5)) +
  xlab("log2(FC)") + ylab("-log10(qval)") +
  theme_bw() +
  theme(legend.position="none") +
  scale_fill_manual(values = c(A= "red", B="blue", C="grey")) +
  scale_colour_manual(values = c(A= "red", B="blue", C="grey"))

## Add gene names and IDs to results_transcripts data.frame 
results_transcripts <- data.frame(geneNames=ballgown::geneNames(bg_filt),
                                  geneIDs=ballgown::geneIDs(bg_filt), 
                                  results_transcripts)

## Write results into a csv ##
write.csv(results_transcripts, "all_transcripts.csv")

## Get results_genes with gene names ##
indices <- match(results_genes$id, texpr(bg, 'all')$gene_id)
gene_names_for_result <- texpr(bg, 'all')$gene_name[indices]
results_genes_wnames <- data.frame(geneNames = gene_names_for_result, 
                                   results_genes)

## Filter for significant (qval<0.5) genes
# and for genes with names only
# Then arrange by ascending qval
sig_genes <- results_genes_wnames[sig, ]
name <- which(sig_genes$geneNames != ".")
sig_genes <- sig_genes[name, ]
sig_genes <- arrange(sig_genes, qval)
write.csv(results_genes_wnames, "all_genes.csv")
write.csv(sig_genes, "sig_genes.csv")

## Plot11: Heatmap (from Mingyang and Alan) based on log2(FPKM)
# Plot only top 50 genes by lowest qvals AND
# top 50 genes by highest fc
top_sig_genes <- sig_genes[1:50, "id"]
top_sig_gene_names <- sig_genes[1:50, "geneNames"]

fc_sig_genes <- arrange(sig_genes, fc)[1:50, "id"]
fc_sig_gene_names <- arrange(sig_genes, fc)[1:50, "geneNames"]

mydist <- function(c) {dist(c,method="euclidian")}
myclust <- function(c) {hclust(c,method="average")}

main_title <- "sig DE Transcripts"
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
par(cex.main=0.8)
data <- log2(as.matrix(gene_expression[fc_sig_genes,1:18])+1)
heatmap.2(data, hclustfun=myclust, distfun=mydist, na.rm = TRUE, scale="none", 
          dendrogram="both", margins=c(6,7), Rowv=TRUE, Colv=TRUE, symbreaks=FALSE, 
          key=TRUE, symkey=FALSE, density.info="none", trace="none", main=main_title, 
          cexRow=1, cexCol=1, labRow=fc_sig_gene_names, col=rev(morecols(100)))

## Look at individual FPKMs ##
fpkm <- texpr(bg, meas="FPKM")
which(ballgown::geneNames(bg) == "DptA")
fpkm[20127, ]

