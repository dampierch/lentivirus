######
######
###### Expression Analysis: Transcripts
######
######


# setup

## libraries
library(tximport)
library(DESeq2)
library(org.Hs.eg.db)
library(edgeR)
library(ggplot2)

## colData
colDat <- read.table("/nv/vol326/cphg_caseylab/lentivirus/dict.csv", sep=",", header=TRUE)
colDat$line <- as.factor(substr(colDat[,"casey_id"],1,3))
colDat$side <- as.factor(substr(colDat[,"casey_id"],4,4))
colDat$condition <- as.factor(substr(colDat[,"casey_id"],6,8))
colDat$passage <- as.factor(substr(colDat[,"casey_id"],10,12))
  # these factors are not production grade yet
  #colDat$race <- 210 = white, 218 = black, 129 = white
  #colDat$sex <- 210 = female, 218 = female, 129 = female
  #colDat$age <- 210 = 47, 218 = 50, 129 = 43

## countData

### define path to quant files
quant_path <- "/nv/vol326/cphg_caseylab/lentivirus/raw_data"
file_pre <- colDat$nwgc_id
file_suf <- "transcriptome_hits.merged.isoforms.results"
file_paths <- file.path(quant_path, paste(file_pre, file_suf, sep="."))
names(file_paths) <- colDat$nwgc_id
all(file.exists(file_paths))
### convert to counts
txi <- tximport(file_paths, type="rsem", txIn=TRUE, txOut=TRUE, countsFromAbundance=c("no"), txIdCol=c("transcript_id"), abundanceCol=c("TPM"), countsCol=c("expected_count"), lengthCol=c("effective_length"))

## design
design <- formula(~ line + condition)

## data object
dds <- DESeqDataSetFromTximport(txi, colData=colDat, design=design)

## pre-filter lowly expressed tx to speed up downstream analysis
  # want 5-10 counts in a library to prove tx is expressed in that library
  # want minimum expression in at least the number of libraries in the smallest group of interest or 1/3 of samples
  # strict filtering performed internally later; this is simply for size reduction and speed
countLimit <- 5
s_num <- ncol(dds)
s_lim <- (1/3) * s_num
keep <- rowSums( counts(dds) > countLimit ) >= s_lim
dds <- dds[keep, ]


# analyze
dds <- DESeq(dds)
file <- paste0("/scratch/chd5n/lentivirus/TxFit.Rdata")
save(dds, file=file)


# test for DE
test_groups <- c("Par", "CC1", "CC2")
res <- list()
for (i in 1:length(test_groups)) {
  res[[i]] <- results(dds, alpha=0.1, contrast=c("condition", test_groups[i], "SCR"))
}
names(res) <- test_groups


# extract results

## summarize DE results
lapply(res, summary)

## rank and report top DE transcripts
resRank <- list()
for (i in 1:length(res)) {
  resRank[[i]] <- res[[i]][order(res[[i]][,"padj"]), ]
  resRank[[i]][,"rank"] <- seq(1,nrow(resRank[[i]]))
}
names(resRank) <- test_groups
file <- paste0("/scratch/chd5n/lentivirus/TxRes.Rdata")
save(resRank, file=file)

for (i in 1:length(resRank)) {
  file <- paste0("/scratch/chd5n/lentivirus/TxRes",names(resRank)[i],".tsv")
  write.table(resRank[[i]], file=file, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
}

## check for rank of target transcripts
colca1_tx <- c("ENST00000620864.1", "ENST00000540738.3", "ENST00000355430.4", "ENST00000532918.4", "ENST00000526150.1")
colca2_tx <- c("ENST00000398035.6", "ENST00000526216.1", "ENST00000614153.4", "ENST00000610738.5", "ENST00000638573.1", "ENST00000528846.5", "ENST00000639470.1")
colca1_rnk <- list()
colca2_rnk <- list()

for (i in 1:length(resRank)) {
  colca1_rnk[[i]] <- vector()
  for (j in 1:length(colca1_tx)) {
    if (length(which(rownames(resRank[[i]])==colca1_tx[j])>0)) {
      #colca1_rnk[[i]] <- append(colca1_rnk[[i]], which(rownames(resRank[[i]])==colca1_tx[j]))
      colca1_rnk[[i]] <- append(colca1_rnk[[i]], which(rownames(resRank[[i]])==colca1_tx[j]))
    } else {
      colca1_rnk[[i]] <- append(colca1_rnk[[i]], NA)
    }
  }
  colca2_rnk[[i]] <- vector()
  for (j in 1:length(colca2_tx)) {
    if (length(which(rownames(resRank[[i]])==colca2_tx[j])>0)) {
      colca2_rnk[[i]] <- append(colca2_rnk[[i]], which(rownames(resRank[[i]])==colca2_tx[j]))
    } else {
      colca2_rnk[[i]] <- append(colca2_rnk[[i]], NA)
    }
  }
}
names(colca1_rnk) <- names(resRank)
names(colca2_rnk) <- names(resRank)


#################
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
