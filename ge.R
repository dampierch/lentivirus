######
######
###### Expression Analysis: Genes
######
######


# setup

## libraries
library(DESeq2)
library(org.Hs.eg.db)
library(edgeR)
library(ggplot2)

## options
cond <- "Par"
cond_name <- "par"

## colData
colDat <- read.table("/nv/vol326/cphg_caseylab/lentivirus/lenti_colca_counts.csv", sep=",", header=TRUE)
colDat <- data.frame(t(colDat))
colDat <- colDat[2:nrow(colDat),]
colDat$line <- as.factor(substr(rownames(colDat),2,4))
colDat$condition <- as.factor(substr(rownames(colDat),7,9))
colDat$side <- as.factor(substr(rownames(colDat),5,5))
colDat$vector <- as.factor(substr(rownames(colDat),11,13))
keep_rows <- colDat[,"condition"]==cond | colDat[,"condition"]=="SCR"
keep_cols <- c("line","condition","side","vector")
colDat <- colDat[keep_rows,keep_cols]

## countData
cts <- read.table("/nv/vol326/cphg_caseylab/lentivirus/Lenti_counts.csv", sep=",", header=TRUE, row.names=1)
keep_cols <- rownames(colDat)
cts <- cts[,keep_cols]

## design
design <- formula(~ line + condition)

## data object
dds <- DESeqDataSetFromMatrix(countData=cts, colData=colDat, design=design)

## pre-filter lowly expressed genes to speed up downstream analysis
  # want 5-10 counts in a library to prove gene is expressed in that library
  # want minimum expression in at least the number of libraries in the smallest group of interest or 1/3 of samples
  # strict filtering performed internally later; this is simply for size reduction and speed
countLimit <- 10
s_num <- ncol(dds)
s_lim <- (1/2) * s_num
keep <- rowSums( counts(dds) > countLimit ) >= s_lim
dds <- dds[keep, ]


# analyze
dds <- DESeq(dds)
file <- paste0("/scratch/chd5n/lentivirus/",cond_name,"Fit.Rdata")
save(dds, file=file)


# test for DE
res <- results(dds, alpha=0.05)
res$symbol <- mapIds(org.Hs.eg.db, keys=row.names(res), column=c("SYMBOL"), keytype=c("ENSEMBL"), multiVals=c("first"))
res <- res[ , c(7,1:6)]


# extract gene lists

## summarize DE results
summary(res)
print(sum(res$padj < 0.05, na.rm=TRUE))

## rank and report top DE genes
resRank <- res[order(res$padj), ]
resRank$rank <- seq(1,nrow(resRank))
file <- paste0("/scratch/chd5n/lentivirus/",cond_name,"Res.tsv")
write.table(resRank, file=file, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)


# viz

## vst
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("line", "condition"), returnData=TRUE)
percentVar <- round( 100 * attr(pcaData, "percentVar") )

## function
pca_double <- function(pcaData, aes_name, ig_name) {
  color <- sym(aes_name[1])
  shape <- sym(aes_name[2])
  sz=5
  percentVar <- percentVar
  condition <- cond
  main <- paste("PCA plot for",condition,"vs SCR")
  ggplot( pcaData, aes(PC1, PC2, color=!!color, shape=!!shape) ) +
    geom_point(size=sz) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() + theme_classic() +
    theme(plot.title=element_text(hjust=0.5)) +
    ggtitle(main)
}

## open device to save plots
file <- paste0("/scratch/chd5n/lentivirus/",cond_name,"Res.pdf")
pdf(file=file)

## pca
print(pca_double(pcaData, c("line", "condition"), "Tx and Line"))

# ## heatmap of count matrix
# select <- order( rowMeans( counts(ddsFND,normalized=TRUE) ), decreasing=TRUE)[1:500]
# df <- as.data.frame(colData(ddsFND)[,c("center","sampTypeCln")])
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df, main="Heatmap, no cluster", col=colors)
# pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, main="Heatmap, cluster by sampType, center", col=colors)
# pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE, annotation_col=df, main="Heatmap, cluster by sampType, center, gene", col=colors)

## biological coefficient of variation (concept from edgeR)
  # the variance relative to the mean = BCV = sqrt(dispersion)
  # genewise BCV vs gene abundance (in log2 cpm)
D <- dispersions(dds)
C <- counts(dds, normalized=TRUE)
A <- aveLogCPM(C, dispersion=D)
plot(A, sqrt(D), pch=16, cex=0.2, main="Genewise biological coefficient of variation vs abundance", xlab="Abundance (avg log CPM)", ylab="Biological coefficient of variation")

## dispersions
plotDispEsts(dds, main="Genewise dispersion vs abundance")

## MA plot
  # M=log ratio (difference between log counts)
  # A=mean average (average of individual count averages)
DESeq2::plotMA(res, ylim=c(-5,5), main="MA plot for Wald test on DE >0, FDR 5%")

## counts for top gene
colors <- c("red","blue","purple")
plotCounts(dds, gene=which.min(res$padj), intgroup="condition", main=paste("Counts for", res[which.min(res$padj),"symbol"]), col=colors[colData(dds)[,"line"]], cex=2)
plotCounts(dds, gene=c("ENSG00000196167"), intgroup="condition", main=paste("Counts for", res["ENSG00000196167","symbol"]), col=colors[colData(dds)[,"line"]], cex=2)
plotCounts(dds, gene=c("ENSG00000214290"), intgroup="condition", main=paste("Counts for", res["ENSG00000214290","symbol"]), col=colors[colData(dds)[,"line"]], cex=2)

## close device to save plots
dev.off()


# revision of plots for top genes/genes of interest

cntData <- plotCounts(dds, gene=c("ENSG00000214290"), intgroup=c("condition","line","side","vector"), returnData=TRUE)

ggplot(cntData, aes(condition, count, color=line, shape=vector)) +
  geom_point(size=4) +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5)) +
  ggtitle(paste("Counts for", res["ENSG00000214290","symbol"]))

main=paste("Counts for", res["ENSG00000196167","symbol"]), col=colors[colData(dds)[,"line"]], cex=2)
