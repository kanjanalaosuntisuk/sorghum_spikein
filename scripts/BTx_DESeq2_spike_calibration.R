##### scripts for DESeq2 analysis with spike-in normalization (Athanasiadou et al., 2019)

# clear workspace
rm(list = ls())
# load packages
library(tidyverse)
library(EDASeq)
library(edgeR)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
setwd("~/sorghum_spikein")

####################################################
##### spike counts ---------------------------------
# load spike-in read count table
BTxspike <- readRDS("RData/BTxspike_new.RData")
# filter out non-detected spikes;  keep reads > 0 in at least 1 samples
spike_filter <- apply(BTxspike, 1, function(x) length(x[x>0])>=1) 
spike_filtered <-  BTxspike[spike_filter,] 
spikes <- rownames(spike_filtered)

##### gene read count data ------------------------
# load gene read count table
BTxgene <- readRDS("RData/BTxgene.RData")
# make a design matrix and DGEList object
conditions <- as.factor(rep(c("controlam", "controlpm", "chillingam", "chillingpm"), each=4))
conditions <- factor(conditions, levels = c("controlam", "controlpm", "chillingam", "chillingpm"))
design <- model.matrix(~ 0 + conditions)
y <- DGEList(BTxgene, group = conditions)
# filter low-expressed genes
keep.exprs <- filterByExpr(y, design)
y_a <- y[keep.exprs,, keep.lib.sizes=FALSE]
# get filterd dataframe
genefiltered <- y_a$counts
genefilteredList <- rownames(genefiltered)

####################################################
# read spikein_filtered read counts ------------------
matSI <- as.matrix(spike_filtered)
# read gene_filtered read counts
matRNA <- genefiltered

# read spike-in concentration txt files
spike_conc <- read.delim("spike_conc_BTx.txt", header = TRUE)
str(spike_conc)
spike_conc %>%
  column_to_rownames(var = "spikeID") -> spike_conc

# select spike_filtered in spike_conc
spikeconc <- spike_conc[rownames(matSI),]

# inspect library size
LibSize_spike <- colSums(matSI) # spike-in library sizes
LibSize_RNA <- colSums(matRNA) # RNA library sizes
Ratio <- LibSize_spike/LibSize_RNA
data.frame(LibSize_spike, LibSize_RNA, Ratio)

# extract attomol of spikes
Nvec <- spikeconc[,"attomol"]
names(Nvec) <- rownames(spikeconc)

##### compute calibration factor using spike-in proportions
# fraction of total spike-in counts accounted for by each spike-in
f.vec <- rowSums(matSI)/sum(matSI) 
# index of spike-in with largest proportion, to be used as the reference spike-in
INDmax <- which(f.vec == max(f.vec)) 
# fraction (empirical proportion) for reference
f.ref <- print(f.vec[INDmax]) 
# get corresponding amol of f.reff
n.ref <- print(Nvec[names(f.vec[INDmax])])
# for counts to nominal amol conversion
nu <- (f.ref/n.ref)*LibSize_spike 
# normalization giving amol
Zmat <- matRNA%*%diag(1/nu) 
colnames(Zmat) <- colnames(matRNA)

# add calibration constant to matRNA
matRNAplus1 <- matRNA + 1 
matPlus1 <- matRNAplus1%*%diag(1/nu) 
colnames(matPlus1) <- colnames(Zmat)

# RLE plots before correcting for unwanted variation
colors <- brewer.pal(4, "Set2") 
colLib <- rep(colors, each=4)
plotRLE(matPlus1, col = colLib, outline = FALSE, ylim = c(-4, 4))
plotPCA(matPlus1, col=colLib, cex=1.2)

##### Compute δj correction factors to account for library preparation errors and apply to νj normalized abundances
meanLogZplus1 <- t(apply(log(matPlus1), 1, function(x) tapply(x, conditions, mean)))
MeanMatLogZplus1 <- t(apply(meanLogZplus1, 1, function(x) rep(x, each=4)))
aux <- colMeans(log(matPlus1) - MeanMatLogZplus1)
# The correction factors
delta <- exp(aux) 
Zmat <- Zmat%*%diag(1/delta)
colnames(Zmat) <- colnames(matRNA)

# RLE plots for adjusted abundance values, counts normalized by (νjδj)
matPlus1Corrected <- matPlus1%*%diag(1/delta)
colnames(matPlus1Corrected) <- colnames(Zmat)
plotRLE(matPlus1Corrected, col=colLib, outline=FALSE, ylim=c(-4, 4)) 
plotPCA(matPlus1Corrected, col=colLib, cex=1.2)

####################################################################
##### DESeq2 DE analysis -------------------------------------------
# calculate normalization factor
library(EnvStats)
nu
delta
nudelta <- nu*delta
geomean <- geoMean(nudelta, na.rm = FALSE)
nudelta_cal <- nudelta/geomean
# make DESeqDataSet object
geneset <- newSeqExpressionSet(as.matrix(genefiltered), 
                                phenoData = data.frame(conditions, row.names=colnames(genefiltered)))
genedds <- DESeqDataSetFromMatrix(countData = counts(geneset), 
                                   colData = pData(geneset), design = ~ conditions)

# add the size factors of spike-ins to gene_dds
sizeFactors(genedds)
sizeFactors(genedds) <- nudelta_cal
colData(genedds)

# run differential expression test
genedds2 <- estimateDispersions(genedds)
genedds3 <- nbinomWaldTest(genedds2)
resultsNames(genedds3)

# save RData
saveRDS(genedds3, file = "RData/DESeq2_nudelta_genedds.RData")

# make RLE plot from normalized counts
nudeltanormalized_counts <- DESeq2::counts(genedds3, normalized=TRUE)
nudeltanormalized_counts_rc <- nudeltanormalized_counts/geomean
x <- as.factor(rep(c("controlam", "controlpm", "chillingam", "chillingpm"), each = 4)) 
x <- factor(x, levels = c("controlam", "controlpm", "chillingam", "chillingpm"))
colors <- brewer.pal(4, "Set2") 
samplenames <- colnames(nudeltanormalized_counts)
# save RData
saveRDS(nudeltanormalized_counts_rc, file = "RData/DESeq2_nudelta_normcounts.RData")

# save plot
png(filename="DESeq2_nudelta_RLE.png", width=7, height=4, units="in", res=300)
plotRLE(nudeltanormalized_counts_rc, outline=FALSE, ylim=c(-4, 4), col=colors[x],
        xaxt = "n", ylab = "Relative log expression") 
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
dev.off()
png(filename="DESeq2_nudelta_PCA.png", width=7, height=4, units="in", res=300)
plotPCA(nudeltanormalized_counts_rc, col=colors[x], cex=0.75)
dev.off()

png(filename="DESeq2_nudelta_RLE_new.png", width=7, height=4, units="in", res=300)
plotRLE(nudeltanormalized_counts, outline=FALSE, ylim=c(-4, 4), col=colors[x], 
        xaxt = "n", ylab = "Relative log expression") 
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
dev.off()
png(filename="DESeq2_nudelta_PCA_new.png", width=7, height=4, units="in", res=300)
plotPCA(nudeltanormalized_counts, col=colors[x], cex=0.75)
dev.off()

# make side-by-side RLE and PCA plot
png(filename="DESeq2_nudelta_RLE_PCAplot.png", width=7, height=3, units="in", res=300)
par(mfrow = c(1, 2))
par(mar = c(5, 4, 1, 1))
plotRLE(nudeltanormalized_counts, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression", main = "Spike-in normalization", cex.main=1, 
        cex.lab = 1, cex.axis=0.8)
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.8)
plotPCA(nudeltanormalized_counts, col=colors[x], cex=0.7, cex.lab = 1, cex.axis = 0.8)
dev.off()

##############################################################
##### controlam VS controlpm -------------------------------------------
# get controlam vs controlpm results with FDR < 0.05 and logFC 0.5
ctampm_nudeltares <- results(genedds3, contrast= c("conditions","controlpm","controlam"))
ctampm_nudeltares005 <- subset(ctampm_nudeltares, ctampm_nudeltares$padj < 0.05)
ctampm_nudeltares005lfcup <- subset(ctampm_nudeltares005, ctampm_nudeltares005$log2FoldChange >= 0.5)
ctampm_nudeltares005lfcdown <- subset(ctampm_nudeltares005, ctampm_nudeltares005$log2FoldChange <= -0.5)
ctampm_nudeltares005lfc <- rbind(ctampm_nudeltares005lfcup, ctampm_nudeltares005lfcdown)

# Log fold change shrinkage for visualization and ranking
resLFCnudelta <- lfcShrink(genedds3, coef = "conditions_controlpm_vs_controlam", type="apeglm")
resLFCnudelta005 <- subset(resLFCnudelta, resLFCnudelta$padj < 0.05)
resLFCnudelta005lfcup <- subset(resLFCnudelta005 , resLFCnudelta005$log2FoldChange >= 0.5)
resLFCnudelta005lfcdown <- subset(resLFCnudelta005, resLFCnudelta005$log2FoldChange <= -0.5)
resLFCnudelta005lfc <- rbind(resLFCnudelta005lfcup, resLFCnudelta005lfcdown)

# get gene list
resLFCnudelta005lfcuplist <- rownames(resLFCnudelta005lfcup)
resLFCnudelta005lfcdownlist <- rownames(resLFCnudelta005lfcdown)
resLFCnudelta005lfclist <- rownames(resLFCnudelta005lfc)

# export gene list
write.table(resLFCnudelta005lfcuplist, file = "output/DEGs_DESeq2_nudelta_ctrlampm_PM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(resLFCnudelta005lfcdownlist, file = "output/DEGs_DESeq2_nudelta_ctrlampm_AM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# save RData
saveRDS(resLFCnudelta005lfc, file = "RData/DEGs_nudelta_full_ctrlampm.RData")
saveRDS(resLFCnudelta, file = "RData/DEGs_nudelta_full_ctrlampm_full.RData")
# export logFC data
#resLFCnudelta005lfc <- readRDS(file = "RData/DEGs_nudelta_full_ctrlampm.RData")
resLFCnudelta005lfc <- as.data.frame(resLFCnudelta005lfc)
resLFCnudelta005lfc %>%
  rownames_to_column(var = "locusID") -> resLFCnudelta005lfc
write.table(resLFCnudelta005lfc, file = "output/DEGs_DESeq2_nudelta_ctrlampm_full.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

##############################################################
##### controlam VS chillingam -------------------------------------------
# get controlam vs chillingam results with FDR < 0.05 and logFC 0.5
trtam_nudeltares <- results(genedds3, contrast= c("conditions","chillingam","controlam"))
trtam_nudeltares005 <- subset(trtam_nudeltares, trtam_nudeltares$padj < 0.05)
trtam_nudeltares005lfcup <- subset(trtam_nudeltares005, trtam_nudeltares005$log2FoldChange >= 0.5)
trtam_nudeltares005lfcdown <- subset(trtam_nudeltares005, trtam_nudeltares005$log2FoldChange <= -0.5)
trtam_nudeltares005lfc <- rbind(trtam_nudeltares005lfcup, trtam_nudeltares005lfcdown)

# Log fold change shrinkage for visualization and ranking
trtam_nudelta <- lfcShrink(genedds3, coef = "conditions_chillingam_vs_controlam", type="apeglm")
trtam_nudelta005 <- subset(trtam_nudelta, trtam_nudelta$padj < 0.05)
trtam_nudelta005lfcup <- subset(trtam_nudelta005 , trtam_nudelta005$log2FoldChange >= 0.5)
trtam_nudelta005lfcdown <- subset(trtam_nudelta005, trtam_nudelta005$log2FoldChange <= -0.5)
trtam_nudelta005lfc <- rbind(trtam_nudelta005lfcup, trtam_nudelta005lfcdown)

# get gene list
trtam_nudelta005lfcuplist <- rownames(trtam_nudelta005lfcup)
trtam_nudelta005lfcdownlist <- rownames(trtam_nudelta005lfcdown)
trtam_nudelta005lfclist <- rownames(trtam_nudelta005lfc)

# export gene list
write.table(trtam_nudelta005lfcuplist, file = "output/DEGs_DESeq2_nudelta_ctrlcoldam_cold.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(trtam_nudelta005lfcdownlist, file = "output/DEGs_DESeq2_nudelta_ctrlcoldam_ctrl.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# save RData
saveRDS(trtam_nudelta005lfc, file = "RData/DEGs_nudelta_full_ctrlcoldam.RData")

###########################################################
##### coldam VS coldpm -------------------------------------------
# relevel condition
genedds3$conditions
genedds3$conditions <- relevel(genedds3$conditions, ref = "chillingam")
# run differential expression test
genedds3 <- nbinomWaldTest(genedds3)
resultsNames(genedds3)

# get coldam vs coldpm results with FDR < 0.05 and logFC 0.5
cdampm_nudeltares <- results(genedds3, contrast= c("conditions","chillingpm","chillingam"))
cdampm_nudeltares005 <- subset(cdampm_nudeltares, cdampm_nudeltares$padj <= 0.05)
cdampm_nudeltares005lfcup <- subset(cdampm_nudeltares005, cdampm_nudeltares005$log2FoldChange >= 0.5)
cdampm_nudeltares005lfcdown <- subset(cdampm_nudeltares005, cdampm_nudeltares005$log2FoldChange <= -0.5)
cdampm_nudeltares005lfc <- rbind(cdampm_nudeltares005lfcup, cdampm_nudeltares005lfcdown)

# Log fold change shrinkage for visualization and ranking
cdampm_resLFCnudelta <- lfcShrink(genedds3, coef = "conditions_chillingpm_vs_chillingam", type="apeglm")
cdampm_resLFCnudelta005 <- subset(cdampm_resLFCnudelta, cdampm_resLFCnudelta$padj < 0.05)
cdampm_resLFCnudelta005lfcup <- subset(cdampm_resLFCnudelta005 , cdampm_resLFCnudelta005$log2FoldChange >= 0.5)
cdampm_resLFCnudelta005lfcdown <- subset(cdampm_resLFCnudelta005, cdampm_resLFCnudelta005$log2FoldChange <= -0.5)
cdampm_resLFCnudelta005lfc <- rbind(cdampm_resLFCnudelta005lfcup, cdampm_resLFCnudelta005lfcdown)

# get gene list
cdampm_resLFCnudelta005lfcuplist <- rownames(cdampm_resLFCnudelta005lfcup)
cdampm_resLFCnudelta005lfcdownlist <- rownames(cdampm_resLFCnudelta005lfcdown)
cdampm_resLFCnudelta005lfclist <- rownames(cdampm_resLFCnudelta005lfc)

# export gene list
write.table(cdampm_resLFCnudelta005lfcuplist, file = "output/DEGs_DESeq2_nudelta_coldampm_PM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cdampm_resLFCnudelta005lfcdownlist, file = "output/DEGs_DESeq2_nudelta_coldampm_AM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# save RData
saveRDS(cdampm_resLFCnudelta005lfc, file = "RData/DEGs_nudelta_full_coldampm.RData")
saveRDS(cdampm_resLFCnudelta, file = "RData/DEGs_nudelta_full_coldampm_full.RData")
# export logFC data
#cdampm_resLFCnudelta005lfc <- readRDS(file = "RData/DEGs_nudelta_full_coldampm.RData")
cdampm_resLFCnudelta005lfc <- as.data.frame(cdampm_resLFCnudelta005lfc)
cdampm_resLFCnudelta005lfc %>%
  rownames_to_column(var = "locusID") -> cdampm_resLFCnudelta005lfc
write.table(cdampm_resLFCnudelta005lfc, file = "output/DEGs_DESeq2_nudelta_coldampm_full.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

###########################################################
##### controlpm VS chillingpm -------------------------------------------
# relevel condition
genedds3$conditions
genedds3$conditions <- relevel(genedds3$conditions, ref = "controlpm")
# run differential expression test
genedds3 <- nbinomWaldTest(genedds3)
resultsNames(genedds3)

# get coldam vs coldpm results with FDR < 0.05 and logFC 0.5
trtpm_nudeltares <- results(genedds3, contrast= c("conditions","chillingpm","controlpm"))
trtpm_nudeltares005 <- subset(trtpm_nudeltares, trtpm_nudeltares$padj <= 0.05)
trtpm_nudeltares005lfcup <- subset(trtpm_nudeltares005, trtpm_nudeltares005$log2FoldChange >= 0.5)
trtpm_nudeltares005lfcdown <- subset(trtpm_nudeltares005, trtpm_nudeltares005$log2FoldChange <= -0.5)
trtpm_nudeltares005lfc <- rbind(trtpm_nudeltares005lfcup, trtpm_nudeltares005lfcdown)

# Log fold change shrinkage for visualization and ranking
trtpm_nudelta <- lfcShrink(genedds3, coef = "conditions_chillingpm_vs_controlpm", type="apeglm")
trtpm_nudelta005 <- subset(trtpm_nudelta, trtpm_nudelta$padj < 0.05)
trtpm_nudelta005lfcup <- subset(trtpm_nudelta005 , trtpm_nudelta005$log2FoldChange >= 0.5)
trtpm_nudelta005lfcdown <- subset(trtpm_nudelta005, trtpm_nudelta005$log2FoldChange <= -0.5)
trtpm_nudelta005lfc <- rbind(trtpm_nudelta005lfcup, trtpm_nudelta005lfcdown)

# get gene list
trtpm_nudelta005lfcuplist <- rownames(trtpm_nudelta005lfcup)
trtpm_nudelta005lfcdownlist <- rownames(trtpm_nudelta005lfcdown)
trtpm_nudelta005lfclist <- rownames(trtpm_nudelta005lfc)

# export gene list
write.table(trtpm_nudelta005lfcuplist, file = "output/DEGs_DESeq2_nudelta_ctrlcoldpm_cold.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(trtpm_nudelta005lfcdownlist, file = "output/DEGs_DESeq2_nudelta_ctrlcoldpm_ctrl.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# save RData
saveRDS(trtpm_nudelta005lfc, file = "RData/DEGs_nudelta_full_ctrlcoldpm.RData")

##########################################################
#### compare DEGs
resLFCnudelta005lfclist # controlam vs controlpm
trtam_nudelta005lfclist # controlam vs chillingam
cdampm_resLFCnudelta005lfclist # chillingam vs chillingpm
trtpm_nudelta005lfclist # controlpm vs chillingpm

# make a venn diagram
gplots::venn(list(ctrlAM_PM = resLFCnudelta005lfclist,
                  ctrl_chillAM = trtam_nudelta005lfclist,
                  chillAM_PM = cdampm_resLFCnudelta005lfclist,
                  ctrl_chillPM = trtpm_nudelta005lfclist))

# make a venn diagram to compare DE genes
library(VennDiagram)
colors <- brewer.pal(4, "Pastel1")
venn.diagram(list(set1 = resLFCnudelta005lfclist,
                  set2 = trtam_nudelta005lfclist,
                  set3 = cdampm_resLFCnudelta005lfclist,
                  set4 = trtpm_nudelta005lfclist), 
             filename = "./venn_DESeq2_nudelta_DEGs.png", 
             fill = colors, lwd = 2, lty = 'blank',
             category.names = c("AM-PM in ctrl", "ctrl-chill in AM", 
                                "AM-PM in chill", "ctrl-chill in PM"),
             cat.pos = c(350, 10, 0, 0),
             height = 9, width = 9, units = 'in')
