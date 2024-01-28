##### scripts for DESeq2 analysis with RUV normalization

# clear workspace
rm(list = ls())
# load packages
#BiocManager::install("apeglm")
library(tidyverse)
library(RUVSeq)
library(DESeq2)
library(apeglm)
setwd("~/sorghum_spikein")

####################################################################
##### RUV normalization
# load spike-in read count table
BTxspike <- readRDS("RData/BTxspike_new.RData")
# filter out non-detected spikes;  keep reads > 0 in at least 1 samples
spike_filter <- apply(BTxspike, 1, function(x) length(x[x>0])>=1) 
spike_filtered <-  BTxspike[spike_filter,] 
spikes <- rownames(spike_filtered)

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

# merge gene and spike read counts
BTxgenespikefiltered <- rbind(genefiltered, spike_filtered)

# make a list of columns and set object
RUVg_set <- newSeqExpressionSet(as.matrix(BTxgenespikefiltered), 
                                phenoData = data.frame(conditions, row.names=colnames(BTxgenespikefiltered)))

# normalize with spike-in
RUVg_set <- betweenLaneNormalization(RUVg_set, which="median") 
RUVg_set2 <- RUVg(RUVg_set, spikes, k=1)

# filter out spike-in data
RUVg_set3 <- RUVg_set2[rownames(RUVg_set2) %in% genefilteredList, ]
# make RLE plot
x <- as.factor(rep(c("controlam", "controlpm", "coldam", "coldpm"), each = 4)) 
x <- factor(x, levels = c("controlam", "controlpm", "coldam", "coldpm"))
colors <- brewer.pal(4, "Set2")
plotRLE(RUVg_set3, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(RUVg_set3, col=colors[x], cex=1.2)

# make a dds object
dds_ruv <- DESeqDataSetFromMatrix(countData = counts(RUVg_set3), 
                              colData = pData(RUVg_set3), design = ~ W_1 + conditions)
# run DESeq
dds_ruv2 <- DESeq(dds_ruv) 
resultsNames(dds_ruv2)

# save RData
saveRDS(dds_ruv2, file = "RData/DESeq2_RUV_genedds.RData")

# make RLE plot
normalized_counts <- counts(dds_ruv2, normalized=TRUE)
x <- as.factor(rep(c("controlam", "controlpm", "chillingam", "chillingpm"), each = 4)) 
x <- factor(x, levels = c("controlam", "controlpm", "chillingam", "chillingpm"))
colors <- brewer.pal(4, "Set2")  
samplenames <- colnames(normalized_counts)
plotRLE(normalized_counts, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression")
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
plotPCA(normalized_counts, col=colors[x], cex=0.75)
# save data
write.csv(normalized_counts, file = "normcount_DESeq2_RUV.txt")

# save plot
png(filename="DESeq2_RUV_RLE.png", width=7, height=4, units="in", res=300)
plotRLE(normalized_counts, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression")
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
dev.off()

png(filename="DESeq2_RUV_PCA.png", width=7, height=4, units="in", res=300)
plotPCA(normalized_counts, col=colors[x], cex=0.75)
dev.off()

# make side-by-side RLE and PCA plot
png(filename="DESeq2_RUV_RLE_PCAplot.png", width=7.5, height=3, units="in", res=300)
par(mfrow = c(1, 2))
par(mar = c(5, 4, 1, 1))
plotRLE(normalized_counts, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression", cex.lab = 0.7, cex.axis=0.7)
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
plotPCA(normalized_counts, col=colors[x], cex=0.7, cex.axis = 0.7, cex.lab = 0.7)
dev.off()

##############################################################
##### controlam VS controlpm ---------------------------------------------
# extract controlam vs controlpm results with FDR < 0.05
ctampm_RUVres <- results(dds_ruv2, contrast= c("conditions","controlpm","controlam"))
ctampm_RUVres005 <- subset(ctampm_RUVres, ctampm_RUVres$padj < 0.05)
ctampm_RUVres005lfcup <- subset(ctampm_RUVres005, ctampm_RUVres005$log2FoldChange >= 0.5)
ctampm_RUVres005lfcdown <- subset(ctampm_RUVres005, ctampm_RUVres005$log2FoldChange <= -0.5)
ctampm_RUVres005lfc <- rbind(ctampm_RUVres005lfcup, ctampm_RUVres005lfcdown)

# Log fold change shrinkage for visualization and ranking
resLFC_RUV <- lfcShrink(dds_ruv2, coef = "conditions_controlpm_vs_controlam", type="apeglm")
resLFC_RUV005 <- subset(resLFC_RUV, resLFC_RUV$padj < 0.05)
resLFC_RUV005lfcup <- subset(resLFC_RUV005, resLFC_RUV005$log2FoldChange >= 0.5)
resLFC_RUV005lfcdown <- subset(resLFC_RUV005, resLFC_RUV005$log2FoldChange <= -0.5)
resLFC_RUV005lfc <- rbind(resLFC_RUV005lfcup, resLFC_RUV005lfcdown)

# get gene list
resLFC_RUV005lfcuplist <- rownames(resLFC_RUV005lfcup)
resLFC_RUV005lfcdownlist <- rownames(resLFC_RUV005lfcdown)
resLFC_RUV005lfclist <- rownames(resLFC_RUV005lfc)

# save RData
saveRDS(resLFC_RUV005lfc, file = "RData/DEGs_deseq_RUV_full_ctrlampm.RData")
# export gene list
write.table(resLFC_RUV005lfcuplist, file = "output/DEGs_DESeq2_RUV_ctrlampm_PM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(resLFC_RUV005lfcdownlist, file = "output/DEGs_DESeq2_RUV_ctrlampm_AM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

####################################################################
##### controlam VS chillingam -------------------------------------------------
# extract controlam vs controlpm results with FDR < 0.05 and logFC 0.5
trtam_RUVres <- results(dds_ruv2, contrast= c("conditions","chillingam","controlam"))
trtam_RUVres05 <- subset(trtam_RUVres, trtam_RUVres$padj < 0.05)
trtam_RUVreslfcup <- subset(trtam_RUVres05, trtam_RUVres05$log2FoldChange >= 0.5)
trtam_RUVreslfcdown <- subset(trtam_RUVres05, trtam_RUVres05$log2FoldChange <= -0.5)
trtam_RUVreslfc <- rbind(trtam_RUVreslfcup, trtam_RUVreslfcdown)

# Log fold change shrinkage for visualization and ranking
trtam_RUV <- lfcShrink(dds_ruv2, coef = "conditions_chillingam_vs_controlam", type="apeglm")
trtam_RUV05 <- subset(trtam_RUV, trtam_RUV$padj < 0.05)
trtam_RUV05lfcup <- subset(trtam_RUV05 , trtam_RUV05$log2FoldChange >= 0.5)
trtam_RUV05lfcdown <- subset(trtam_RUV05, trtam_RUV05$log2FoldChange <= -0.5)
trtam_RUV05lfc <- rbind(trtam_RUV05lfcup, trtam_RUV05lfcdown)

# get gene list
trtam_RUV05lfcuplist <- rownames(trtam_RUV05lfcup)
trtam_RUV05lfcdownlist <- rownames(trtam_RUV05lfcdown)
trtam_RUV05lfclist <- rownames(trtam_RUV05lfc)

# save RData
saveRDS(trtam_RUV05lfc, file = "RData/DEGs_deseq_RUV_full_ctrlcoldam.RData")
# export gene list
write.table(trtam_RUV05lfcuplist, file = "output/DEGs_DESeq2_RUV_ctrlcoldam_cold.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(trtam_RUV05lfcdownlist, file = "output/DEGs_DESeq2_RUV_ctrlcoldam_ctrl.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

##############################################################
##### coldam VS coldpm
# re-level condition
dds_ruv2$conditions
dds_ruv2$conditions <- relevel(dds_ruv2$conditions, ref = "chillingam")
# run DESeq
dds_ruv2 <- nbinomWaldTest(dds_ruv2)
resultsNames(dds_ruv2)

# extract coldam vs coldpm results with FDR < 0.05
cdampm_RUVres <- results(dds_ruv2, contrast= c("conditions","chillingpm","chillingam"))
cdampm_RUVres005 <- subset(cdampm_RUVres, cdampm_RUVres$padj < 0.05)
cdampm_RUVres005lfcup <- subset(cdampm_RUVres005, cdampm_RUVres005$log2FoldChange >= 0.5)
cdampm_RUVres005lfcdown <- subset(cdampm_RUVres005, cdampm_RUVres005$log2FoldChange <= -0.5)
cdampm_RUVres005lfc <- rbind(cdampm_RUVres005lfcup, cdampm_RUVres005lfcdown)

# Log fold change shrinkage for visualization and ranking
cdampm_resLFC_RUV <- lfcShrink(dds_ruv2, coef = "conditions_chillingpm_vs_chillingam", type="apeglm")
cdampm_resLFC_RUV005 <- subset(cdampm_resLFC_RUV, cdampm_resLFC_RUV$padj < 0.05)
cdampm_resLFC_RUV005lfcup <- subset(cdampm_resLFC_RUV005, cdampm_resLFC_RUV005$log2FoldChange >= 0.5)
cdampm_resLFC_RUV005lfcdown <- subset(cdampm_resLFC_RUV005, cdampm_resLFC_RUV005$log2FoldChange <= -0.5)
cdampm_resLFC_RUV005lfc <- rbind(cdampm_resLFC_RUV005lfcup, cdampm_resLFC_RUV005lfcdown)

# get gene list
cdampm_resLFC_RUV005lfcuplist <- rownames(cdampm_resLFC_RUV005lfcup)
cdampm_resLFC_RUV005lfcdownlist <- rownames(cdampm_resLFC_RUV005lfcdown)
cdampm_resLFC_RUV005lfclist <- rownames(cdampm_resLFC_RUV005lfc)

# save RData
saveRDS(cdampm_resLFC_RUV005lfc, file = "RData/DEGs_deseq_RUV_full_coldampm.RData")
# export gene list
write.table(cdampm_resLFC_RUV005lfcuplist, file = "output/DEGs_DESeq2_RUV_coldampm_PM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cdampm_resLFC_RUV005lfcdownlist, file = "output/DEGs_DESeq2_RUV_coldampm_AM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

##############################################################
##### controlpm VS chillingpm -------------------------------------------------
# re-level condition
dds_ruv2$conditions
dds_ruv2$conditions <- relevel(dds_ruv2$conditions, ref = "controlpm")
# run DESeq
dds_ruv2 <- nbinomWaldTest(dds_ruv2)
resultsNames(dds_ruv2)

# extract controlam vs controlpm results with FDR < 0.05 and logFC 0.5
trtpm_RUVres <- results(dds_ruv2, contrast= c("conditions","chillingpm","controlpm"))
trtpm_RUVres005 <- subset(trtpm_RUVres, trtpm_RUVres$padj < 0.05)
trtpm_RUVreslfcup <- subset(trtpm_RUVres005, trtpm_RUVres005$log2FoldChange >= 0.5)
trtpm_RUVreslfcdown <- subset(trtpm_RUVres005, trtpm_RUVres005$log2FoldChange <= -0.5)
trtpm_RUVreslfc <- rbind(trtpm_RUVreslfcup, trtpm_RUVreslfcdown)

# Log fold change shrinkage for visualization and ranking
trtpm_RUV <- lfcShrink(dds_ruv2, coef = "conditions_chillingpm_vs_controlpm", type="apeglm")
trtpm_RUV05 <- subset(trtpm_RUV, trtpm_RUV$padj < 0.05)
trtpm_RUV05lfcup <- subset(trtpm_RUV05 , trtpm_RUV05$log2FoldChange >= 0.5)
trtpm_RUV05lfcdown <- subset(trtpm_RUV05, trtpm_RUV05$log2FoldChange <= -0.5)
trtpm_RUV05lfc <- rbind(trtpm_RUV05lfcup, trtpm_RUV05lfcdown)

# get gene list
trtpm_RUV05lfcuplist <- rownames(trtpm_RUV05lfcup)
trtpm_RUV05lfcdownlist <- rownames(trtpm_RUV05lfcdown)
trtpm_RUV05lfclist <- rownames(trtpm_RUV05lfc)

# save RData
saveRDS(trtpm_RUV05lfc, file = "RData/DEGs_deseq_RUV_full_ctrlcoldpm.RData")
# export gene list
write.table(trtpm_RUV05lfcuplist, file = "output/DEGs_DESeq2_RUV_ctrlcoldpm_cold.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(trtpm_RUV05lfcdownlist, file = "output/DEGs_DESeq2_RUV_ctrlcoldpm_ctrl.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


