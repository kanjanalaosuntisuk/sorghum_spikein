##### scripts for DESeq2 analysis with spike-in normalization 
# reference: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# clear workspace
rm(list = ls())
# load packages
# BiocManager::install("apeglm")
library(tidyverse)
library(RUVSeq)
library(DESeq2)
library(apeglm)
setwd("~/sorghum_spikein")

####################################################################
##### spike-in normalization by size factor ----
# load spike-in read count table
BTxspike <- readRDS("RData/BTxspike_new.RData")
# filter out non-detected spikes;  keep reads > 0 in at least 1 samples
spike_filter <- apply(BTxspike, 1, function(x) length(x[x>0])>=1) 
spike_filtered <-  BTxspike[spike_filter,] 
spikes <- rownames(spike_filtered)

# make DESeqDataSet
conditions <- as.factor(rep(c("controlam", "controlpm", "chillingam", "chillingpm"), each=4))
conditions <- factor(conditions, levels = c("controlam", "controlpm", "chillingam", "chillingpm"))
spike_set <- newSeqExpressionSet(as.matrix(spike_filtered), 
                                 phenoData = data.frame(conditions, row.names=colnames(spike_filtered)))
spike_dds <- DESeqDataSetFromMatrix(countData = counts(spike_set), 
                                    colData = pData(spike_set), design = ~ conditions)
# estimate size factor from spike in data
spike_cds <- estimateSizeFactors(spike_dds)
spike_sf <- sizeFactors(spike_cds)

##### gene read count data ----
# load gene read count table
BTxgene <- readRDS("RData/BTxgene.RData")
# make a design matrix and DGEList object
design <- model.matrix(~ 0 + conditions)
y <- DGEList(BTxgene, group = conditions)
# filter low-expressed genes
keep.exprs <- filterByExpr(y, design)
y_a <- y[keep.exprs,, keep.lib.sizes=FALSE]
# get filterd dataframe
genefiltered <- y_a$counts
genefilteredList <- rownames(genefiltered)

# make DESeqDataSet object
gene_set <- newSeqExpressionSet(as.matrix(genefiltered), 
                                phenoData = data.frame(conditions, row.names=colnames(genefiltered)))
gene_dds <- DESeqDataSetFromMatrix(countData = counts(gene_set), 
                                   colData = pData(gene_set), design = ~ conditions)
# add the size factors of spike-ins to gene_dds
sizeFactors(gene_dds)
sizeFactors(gene_dds) <- spike_sf
colData(gene_dds)

# run differential expression test
gene_dds2 <- estimateDispersions(gene_dds)
gene_dds3 <- nbinomWaldTest(gene_dds2)
resultsNames(gene_dds3)

# save RData
saveRDS(gene_dds3, file = "RData/DESeq2_spikesf_genedds.RData")

# make RLE plot
normalized_counts <- counts(gene_dds3, normalized=TRUE)
x <- as.factor(rep(c("controlam", "controlpm", "chillingam", "chillingpm"), each = 4)) 
x <- factor(x, levels = c("controlam", "controlpm", "chillingam", "chillingpm"))
colors <- brewer.pal(4, "Set2")
samplenames <- colnames(normalized_counts)
plotRLE(normalized_counts, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression")
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
plotPCA(normalized_counts, col=colors[x], cex=0.75)
write.csv(normalized_counts, file = "normcount_DESeq2_spikesf.txt")

# save plot
png(filename="DESeq2_spikesf_RLE.png", width=7, height=4, units="in", res=300)
plotRLE(normalized_counts, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression")
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
dev.off()

png(filename="DESeq2_spikesf_PCA.png", width=7, height=4, units="in", res=300)
plotPCA(normalized_counts, col=colors[x], cex=0.75)
dev.off()

# make side-by-side RLE and PCA plot
png(filename="DESeq2_spikesf_RLE_PCAplot.png", width=7.5, height=3, units="in", res=300)
par(mfrow = c(1, 2))
par(mar = c(5, 4, 1, 1))
plotRLE(normalized_counts, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression", cex.lab = 0.7, cex.axis=0.7)
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
plotPCA(normalized_counts, col=colors[x], cex=0.7, cex.axis = 0.7, cex.lab = 0.7)
dev.off()

##############################################################
##### controlam VS controlpm
# get controlam vs controlpm results with FDR < 0.05 and logFC 0.5
ctampm_spikesfres <- results(gene_dds3, contrast= c("conditions","controlpm","controlam"))
ctampm_spikesfres005 <- subset(ctampm_spikesfres, ctampm_spikesfres$padj < 0.05)
ctampm_spikesfres005lfcup <- subset(ctampm_spikesfres005, ctampm_spikesfres005$log2FoldChange >= 0.5)
ctampm_spikesfres005lfcdown <- subset(ctampm_spikesfres005, ctampm_spikesfres005$log2FoldChange <= -0.5)
ctampm_spikesfres005lfc <- rbind(ctampm_spikesfres005lfcup, ctampm_spikesfres005lfcdown)

# Log fold change shrinkage for visualization and ranking
spikesfresLFC <- lfcShrink(gene_dds3, coef = "conditions_controlpm_vs_controlam", type="apeglm")
spikesfresLFC005 <- subset(spikesfresLFC, spikesfresLFC$padj < 0.05)
spikesfresLFC005lfcup <- subset(spikesfresLFC005, spikesfresLFC005$log2FoldChange >= 0.5)
spikesfresLFC005lfcdown <- subset(spikesfresLFC005, spikesfresLFC005$log2FoldChange <= -0.5)
spikesfresLFC005lfc <- rbind(spikesfresLFC005lfcup, spikesfresLFC005lfcdown)

# get gene list
spikesfresLFC005lfcuplist <- rownames(spikesfresLFC005lfcup)
spikesfresLFC005lfcdownlist <- rownames(spikesfresLFC005lfcdown)
spikesfresLFC005lfclist <- rownames(spikesfresLFC005lfc)

# save RData
saveRDS(spikesfresLFC005lfc, file = "RData/DEGs_deseq_spikesf_full_ctrlampm.RData")
# export gene list
write.table(spikesfresLFC005lfcuplist, file = "output/DEGs_DESeq2_spikesf_ctrlampm_PM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(spikesfresLFC005lfcdownlist, file = "output/DEGs_DESeq2_spikesf_ctrlampm_AM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

####################################################################
##### controlam VS chillingam -------------------------------------------------
# extract controlam vs controlpm results with FDR < 0.05 and logFC 0.5
trtam_spikesfres <- results(gene_dds3, contrast= c("conditions","chillingam","controlam"))
trtam_spikesfres005 <- subset(trtam_spikesfres, trtam_spikesfres$padj < 0.05)
trtam_spikesfreslfcup <- subset(trtam_spikesfres005, trtam_spikesfres005$log2FoldChange >= 0.5)
trtam_spikesfreslfcdown <- subset(trtam_spikesfres005, trtam_spikesfres005$log2FoldChange <= -0.5)
trtam_spikesfreslfc <- rbind(trtam_spikesfreslfcup, trtam_spikesfreslfcdown)

# Log fold change shrinkage for visualization and ranking
trtam_spikesf <- lfcShrink(gene_dds3, coef = "conditions_chillingam_vs_controlam", type="apeglm")
trtam_spikesf05 <- subset(trtam_spikesf, trtam_spikesf$padj < 0.05)
trtam_spikesf05lfcup <- subset(trtam_spikesf05 , trtam_spikesf05$log2FoldChange >= 0.5)
trtam_spikesf05lfcdown <- subset(trtam_spikesf05, trtam_spikesf05$log2FoldChange <= -0.5)
trtam_spikesf05lfc <- rbind(trtam_spikesf05lfcup, trtam_spikesf05lfcdown)

# get gene list
trtam_spikesf05lfcuplist <- rownames(trtam_spikesf05lfcup)
trtam_spikesf05lfcdownlist <- rownames(trtam_spikesf05lfcdown)
trtam_spikesf05lfclist <- rownames(trtam_spikesf05lfc)

# save RData
saveRDS(trtam_spikesf05lfc, file = "RData/DEGs_deseq_spikesf_full_ctrlcoldam.RData")
# export gene list
write.table(trtam_spikesf05lfcuplist, file = "output/DEGs_DESeq2_spikesf_ctrlcoldam_cold.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(trtam_spikesf05lfcdownlist, file = "output/DEGs_DESeq2_spikesf_ctrlcoldam_ctrl.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

##############################################################
##### coldam VS coldpm
# re-level condition
gene_dds3$conditions
gene_dds3$conditions <- relevel(gene_dds3$conditions, ref = "chillingam")

# run differential expression test
gene_dds3 <- nbinomWaldTest(gene_dds3)
resultsNames(gene_dds3)

# get coldam vs coldpm results with FDR < 0.05 and logFC 0.5
cdampm_spikesfres <- results(gene_dds3, contrast= c("conditions","chillingpm","chillingam"))
cdampm_spikesfres005 <- subset(cdampm_spikesfres, cdampm_spikesfres$padj < 0.05)
cdampm_spikesfres005lfcup <- subset(cdampm_spikesfres005, cdampm_spikesfres005$log2FoldChange >= 0.5)
cdampm_spikesfres005lfcdown <- subset(cdampm_spikesfres005, cdampm_spikesfres005$log2FoldChange <= -0.5)
cdampm_spikesfres005lfc <- rbind(cdampm_spikesfres005lfcup, cdampm_spikesfres005lfcdown)

# Log fold change shrinkage for visualization and ranking
cdampm_spikesfresLFC <- lfcShrink(gene_dds3, coef = "conditions_chillingpm_vs_chillingam", type="apeglm")
cdampm_spikesfresLFC005 <- subset(cdampm_spikesfresLFC, cdampm_spikesfresLFC$padj < 0.05)
cdampm_spikesfresLFC005lfcup <- subset(cdampm_spikesfresLFC005, cdampm_spikesfresLFC005$log2FoldChange >= 0.5)
cdampm_spikesfresLFC005lfcdown <- subset(cdampm_spikesfresLFC005, cdampm_spikesfresLFC005$log2FoldChange <= -0.5)
cdampm_spikesfresLFC005lfc <- rbind(cdampm_spikesfresLFC005lfcup, cdampm_spikesfresLFC005lfcdown)

# get gene list
cdampm_spikesfresLFC005lfcuplist <- rownames(cdampm_spikesfresLFC005lfcup)
cdampm_spikesfresLFC005lfcdownlist <- rownames(cdampm_spikesfresLFC005lfcdown)
cdampm_spikesfresLFC005lfclist <- rownames(cdampm_spikesfresLFC005lfc)

# save RData
saveRDS(cdampm_spikesfresLFC005lfc, file = "RData/DEGs_deseq_spikesf_full_coldampm.RData")
# export gene list
write.table(cdampm_spikesfresLFC005lfcuplist, file = "output/DEGs_DESeq2_spikesf_coldampm_PM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cdampm_spikesfresLFC005lfcdownlist, file = "output/DEGs_DESeq2_spikesf_coldampm_AM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

##############################################################
##### controlpm VS chillingpm -------------------------------------------------
# re-level condition
gene_dds3$conditions
gene_dds3$conditions <- relevel(gene_dds3$conditions, ref = "controlpm")
# run DESeq
gene_dds3 <- nbinomWaldTest(gene_dds3)
resultsNames(gene_dds3)

# extract controlam vs controlpm results with FDR < 0.05 and logFC 0.5
trtpm_spikesfres <- results(gene_dds3, contrast= c("conditions","chillingpm","controlpm"))
trtpm_spikesfres005 <- subset(trtpm_spikesfres, trtpm_spikesfres$padj < 0.05)
trtpm_spikesfreslfcup <- subset(trtpm_spikesfres005, trtpm_spikesfres005$log2FoldChange >= 0.5)
trtpm_spikesfreslfcdown <- subset(trtpm_spikesfres005, trtpm_spikesfres005$log2FoldChange <= -0.5)
trtpm_spikesfreslfc <- rbind(trtpm_spikesfreslfcup, trtpm_spikesfreslfcdown)

# Log fold change shrinkage for visualization and ranking
trtpm_spikesf <- lfcShrink(gene_dds3, coef = "conditions_chillingpm_vs_controlpm", type="apeglm")
trtpm_spikesf05 <- subset(trtpm_spikesf, trtpm_spikesf$padj < 0.05)
trtpm_spikesf05lfcup <- subset(trtpm_spikesf05 , trtpm_spikesf05$log2FoldChange >= 0.5)
trtpm_spikesf05lfcdown <- subset(trtpm_spikesf05, trtpm_spikesf05$log2FoldChange <= -0.5)
trtpm_spikesf05lfc <- rbind(trtpm_spikesf05lfcup, trtpm_spikesf05lfcdown)

# get gene list
trtpm_spikesf05lfcuplist <- rownames(trtpm_spikesf05lfcup)
trtpm_spikesf05lfcdownlist <- rownames(trtpm_spikesf05lfcdown)
trtpm_spikesf05lfclist <- rownames(trtpm_spikesf05lfc)

# save RData
saveRDS(trtpm_spikesf05lfc, file = "RData/DEGs_deseq_spikesf_full_ctrlcoldpm.RData")
# export gene list
write.table(trtpm_spikesf05lfcuplist, file = "output/DEGs_DESeq2_spikesf_ctrlcoldpm_cold.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(trtpm_spikesf05lfcdownlist, file = "output/DEGs_DESeq2_spikesf_ctrlcoldpm_ctrl.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
