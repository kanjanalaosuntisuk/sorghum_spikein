##### scripts for DESeq2 analysis (Median of Ratio - traditional normalization method)
# reference: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

# clear workspace
rm(list = ls())
# load packages
# BiocManager::install("apeglm")
library(tidyverse)
library(RColorBrewer)
library(RUVSeq)
library(DESeq2)
library(apeglm)
setwd("~/sorghum_spikein")

####################################################################
##### Median of Ratio normalization (no spike-in)
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
# make DESeqDataSet object
gene_set_nonorm <- newSeqExpressionSet(as.matrix(genefiltered), 
                                       phenoData = data.frame(conditions, row.names=colnames(genefiltered)))
gene_dds_nonorm <- DESeqDataSetFromMatrix(countData = counts(gene_set_nonorm), 
                                          colData = pData(gene_set_nonorm), design = ~ conditions)

# run DESeq
sizeFactors(gene_dds_nonorm)
gene_dds_nonorm2 <- DESeq(gene_dds_nonorm) 
resultsNames(gene_dds_nonorm2)

# save RData (containing normalized read counts)
saveRDS(gene_dds_nonorm2, file = "RData/DESeq2_nospike_genedds.RData")

# make RLE plot
normalized_counts_deseq <- counts(gene_dds_nonorm2, normalized = TRUE)
x <- as.factor(rep(c("controlam", "controlpm", "chillingam", "chillingpm"), each = 4)) 
x <- factor(x, levels = c("controlam", "controlpm", "chillingam", "chillingpm"))
colors <- brewer.pal(4, "Set2")
samplenames <- colnames(normalized_counts_deseq)
plotRLE(normalized_counts_deseq, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression")
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
plotPCA(normalized_counts_deseq, col=colors[x], cex=0.75)
write.csv(normalized_counts_deseq, file = "normcount_DESeq2_nospike.txt")

# save plot
png(filename="DESeq2_nospike_RLE.png", width=7, height=4, units="in", res=300)
plotRLE(normalized_counts_deseq, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression")
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
dev.off()

png(filename="DESeq2_nospike_PCA.png", width=7, height=4, units="in", res=300)
plotPCA(normalized_counts_deseq, col=colors[x], cex=0.75)
dev.off()

# make side-by-side RLE and PCA plot
png(filename="DESeq2_nospike_RLE_PCAplot.png", width=7, height=3, units="in", res=300)
par(mfrow = c(1, 2))
par(mar = c(5, 4, 1, 1))
plotRLE(normalized_counts_deseq, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression", main = "Median-of-Ratio normalization", cex.main=1, 
        cex.lab = 1, cex.axis=0.8)
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.8)
plotPCA(normalized_counts_deseq, col=colors[x], cex=0.7, cex.lab = 1, cex.axis = 0.8)
dev.off()

####################################################################
##### controlam VS controlpm -------------------------------------------------
# extract controlam vs controlpm results with FDR (false discovery rate) < 0.05 and logFC 0.5
ctampm_nonormres <- results(gene_dds_nonorm2, contrast= c("conditions","controlpm","controlam"))
ctampm_nonormres005 <- subset(ctampm_nonormres, ctampm_nonormres$padj < 0.05)
ctampm_nonormreslfcup <- subset(ctampm_nonormres005, ctampm_nonormres005$log2FoldChange >= 0.5)
ctampm_nonormreslfcdown <- subset(ctampm_nonormres005, ctampm_nonormres005$log2FoldChange <= -0.5)
ctampm_nonormreslfc <- rbind(ctampm_nonormreslfcup, ctampm_nonormreslfcdown)

# Log fold change shrinkage for visualization and ranking
resLFC_nonorm <- lfcShrink(gene_dds_nonorm2, coef = "conditions_controlpm_vs_controlam", type="apeglm")
resLFC_nonorm005 <- subset(resLFC_nonorm, resLFC_nonorm$padj < 0.05)
resLFC_nonorm005lfcup <- subset(resLFC_nonorm005, resLFC_nonorm005$log2FoldChange >= 0.5)
resLFC_nonorm005lfcdown <- subset(resLFC_nonorm005, resLFC_nonorm005$log2FoldChange <= -0.5)
resLFC_nonorm005lfc <- rbind(resLFC_nonorm005lfcup, resLFC_nonorm005lfcdown)

# get gene list
resLFC_nonorm005lfcuplist <- rownames(resLFC_nonorm005lfcup)
resLFC_nonorm005lfcdownlist <- rownames(resLFC_nonorm005lfcdown)
resLFC_nonorm005lfclist <- rownames(resLFC_nonorm005lfc)

# export gene list
write.table(resLFC_nonorm005lfcuplist, file = "output/DEGs_DESeq2_nospike_ctrlampm_PM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(resLFC_nonorm005lfcdownlist, file = "output/DEGs_DESeq2_nospike_ctrlampm_AM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# save RData
saveRDS(resLFC_nonorm005lfc, file = "RData/DEGs_deseq_full_ctrlampm.RData")
saveRDS(resLFC_nonorm, file = "RData/DEGs_deseq_full_ctrlampm_full.RData")
# export full table
resLFC_nonorm005lfc <- readRDS(file = "RData/DEGs_deseq_full_ctrlampm.RData")
resLFC_nonorm005lfc <- as.data.frame(resLFC_nonorm005lfc)
resLFC_nonorm005lfc %>%
  rownames_to_column(var = "locusID") -> resLFC_nonorm005lfc
write.table(resLFC_nonorm005lfc, file = "output/DEGs_DESeq2_nospike_ctrlampm_full.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

####################################################################
##### controlam VS chillingam -------------------------------------------------
# extract controlam vs controlpm results with FDR < 0.05 and logFC 0.5
trtam_nonormres <- results(gene_dds_nonorm2, contrast= c("conditions","chillingam","controlam"))
trtam_nonormres005 <- subset(trtam_nonormres, trtam_nonormres$padj < 0.05)
trtam_nonormreslfcup <- subset(trtam_nonormres005, trtam_nonormres005$log2FoldChange >= 0.5)
trtam_nonormreslfcdown <- subset(trtam_nonormres005, trtam_nonormres005$log2FoldChange <= -0.5)
trtam_nonormreslfc <- rbind(trtam_nonormreslfcup, trtam_nonormreslfcdown)

# Log fold change shrinkage for visualization and ranking
trtam_nonorm <- lfcShrink(gene_dds_nonorm2, coef = "conditions_chillingam_vs_controlam", type="apeglm")
trtam_nonorm005 <- subset(trtam_nonorm, trtam_nonorm$padj < 0.05)
trtam_nonorm005lfcup <- subset(trtam_nonorm005 , trtam_nonorm005$log2FoldChange >= 0.5)
trtam_nonorm005lfcdown <- subset(trtam_nonorm005, trtam_nonorm005$log2FoldChange <= -0.5)
trtam_nonorm005lfc <- rbind(trtam_nonorm005lfcup, trtam_nonorm005lfcdown)

# get gene list
trtam_nonorm005lfcuplist <- rownames(trtam_nonorm005lfcup)
trtam_nonorm005lfcdownlist <- rownames(trtam_nonorm005lfcdown)
trtam_nonorm005lfclist <- rownames(trtam_nonorm005lfc)

# export gene list
write.table(trtam_nonorm005lfcuplist, file = "output/DEGs_DESeq2_nospike_ctrlcoldam_cold.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(trtam_nonorm005lfcdownlist, file = "output/DEGs_DESeq2_nospike_ctrlcoldam_ctrl.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# save RData
saveRDS(trtam_nonorm005lfc, file = "RData/DEGs_deseq_full_ctrlcoldam.RData")

####################################################################
##### coldam VS coldpm ------------------------------------------------------
# re-level condition
gene_dds_nonorm2
gene_dds_nonorm2$conditions <- relevel(gene_dds_nonorm2$conditions, ref = "chillingam")
# run DESeq
gene_dds_nonorm2 <- nbinomWaldTest(gene_dds_nonorm2)

# extract coldam vs coldpm results with FDR < 0.05 and logFC 0.5
cdampm_nonormres <- results(gene_dds_nonorm2, contrast= c("conditions","chillingpm","chillingam"))
cdampm_nonormres005 <- subset(cdampm_nonormres, cdampm_nonormres$padj < 0.05)
cdampm_nonormreslfcup <- subset(cdampm_nonormres005, cdampm_nonormres005$log2FoldChange >= 0.5)
cdampm_nonormreslfcdown <- subset(cdampm_nonormres005, cdampm_nonormres005$log2FoldChange <= -0.5)
cdampm_nonormreslfc <- rbind(cdampm_nonormreslfcup, cdampm_nonormreslfcdown)

# Log fold change shrinkage for visualization and ranking
cdampm_resLFC_nonorm <- lfcShrink(gene_dds_nonorm2, coef = "conditions_chillingpm_vs_chillingam", type="apeglm")
cdampm_resLFC_nonorm005 <- subset(cdampm_resLFC_nonorm, cdampm_resLFC_nonorm$padj < 0.05)
cdampm_resLFC_nonorm005lfcup <- subset(cdampm_resLFC_nonorm005, cdampm_resLFC_nonorm005$log2FoldChange >= 0.5)
cdampm_resLFC_nonorm005lfcdown <- subset(cdampm_resLFC_nonorm005, cdampm_resLFC_nonorm005$log2FoldChange <= -0.5)
cdampm_resLFC_nonorm005lfc <- rbind(cdampm_resLFC_nonorm005lfcup, cdampm_resLFC_nonorm005lfcdown)

# get gene list
cdampm_resLFC_nonorm005lfcuplist <- rownames(cdampm_resLFC_nonorm005lfcup)
cdampm_resLFC_nonorm005lfcdownlist <- rownames(cdampm_resLFC_nonorm005lfcdown)
cdampm_resLFC_nonorm005lfclist <- rownames(cdampm_resLFC_nonorm005lfc)

# export gene list
write.table(cdampm_resLFC_nonorm005lfcuplist, file = "output/DEGs_DESeq2_nospike_coldampm_PM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(cdampm_resLFC_nonorm005lfcdownlist, file = "output/DEGs_DESeq2_nospike_coldampm_AM.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# save RData
saveRDS(cdampm_resLFC_nonorm005lfc, file = "RData/DEGs_deseq_full_coldampm.RData")
saveRDS(cdampm_resLFC_nonorm, file = "RData/DEGs_deseq_full_coldampm_full.RData")
# export logFC data
cdampm_resLFC_nonorm005lfc <- readRDS(file = "RData/DEGs_deseq_full_coldampm.RData")
cdampm_resLFC_nonorm005lfc <- as.data.frame(cdampm_resLFC_nonorm005lfc)
cdampm_resLFC_nonorm005lfc %>%
  rownames_to_column(var = "locusID") -> cdampm_resLFC_nonorm005lfc
write.table(cdampm_resLFC_nonorm005lfc, file = "output/DEGs_DESeq2_nospike_coldampm_full.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

####################################################################
##### controlpm VS chillingpm -------------------------------------------------
# re-level condition
gene_dds_nonorm2$conditions
gene_dds_nonorm2$conditions <- relevel(gene_dds_nonorm2$conditions, ref = "controlpm")
# run DESeq
gene_dds_nonorm2 <- nbinomWaldTest(gene_dds_nonorm2)
resultsNames(gene_dds_nonorm2)

# extract controlam vs controlpm results with FDR < 0.05 and logFC 0.5
trtpm_nonormres <- results(gene_dds_nonorm2, contrast= c("conditions","chillingpm","controlpm"))
trtpm_nonormres005 <- subset(trtpm_nonormres, trtpm_nonormres$padj < 0.05)
trtpm_nonormreslfcup <- subset(trtpm_nonormres005, trtpm_nonormres005$log2FoldChange >= 0.5)
trtpm_nonormreslfcdown <- subset(trtpm_nonormres005, trtpm_nonormres005$log2FoldChange <= -0.5)
trtpm_nonormreslfc <- rbind(trtpm_nonormreslfcup, trtpm_nonormreslfcdown)

# Log fold change shrinkage for visualization and ranking
trtpm_nonorm <- lfcShrink(gene_dds_nonorm2, coef = "conditions_chillingpm_vs_controlpm", type="apeglm")
trtpm_nonorm005 <- subset(trtpm_nonorm, trtpm_nonorm$padj < 0.05)
trtpm_nonorm005lfcup <- subset(trtpm_nonorm005 , trtpm_nonorm005$log2FoldChange >= 0.5)
trtpm_nonorm005lfcdown <- subset(trtpm_nonorm005, trtpm_nonorm005$log2FoldChange <= -0.5)
trtpm_nonorm005lfc <- rbind(trtpm_nonorm005lfcup, trtpm_nonorm005lfcdown)

# get gene list
trtpm_nonorm005lfcuplist <- rownames(trtpm_nonorm005lfcup)
trtpm_nonorm005lfcdownlist <- rownames(trtpm_nonorm005lfcdown)
trtpm_nonorm005lfclist <- rownames(trtpm_nonorm005lfc)

# export gene list
write.table(trtpm_nonorm005lfcuplist, file = "output/DEGs_DESeq2_nospike_ctrlcoldpm_cold.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(trtpm_nonorm005lfcdownlist, file = "output/DEGs_DESeq2_nospike_ctrlcoldpm_ctrl.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
# save RData
saveRDS(trtpm_nonorm005lfc, file = "RData/DEGs_deseq_full_ctrlcoldpm.RData")

##########################################################
#### compare DEGs
resLFC_nonorm005lfclist # comtrolam vs comtrolpm
trtam_nonorm005lfclist # controlam vs chillingam
cdampm_resLFC_nonorm005lfclist # chillingam vs chillingpm
trtpm_nonorm005lfclist # controlpm vs chillingpm

# make a venn diagram
gplots::venn(list(ctrlAM_PM = resLFC_nonorm005lfclist,
                  ctrl_chillAM = trtam_nonorm005lfclist,
                  chillAM_PM = cdampm_resLFC_nonorm005lfclist,
                  ctrl_chillPM = trtpm_nonorm005lfclist))

# make a venn diagram to compare DE genes
library(VennDiagram)
colors <- brewer.pal(4, "Pastel1")
venn.diagram(list(set1 = resLFC_nonorm005lfclist,
                  set2 = trtam_nonorm005lfclist,
                  set3 = cdampm_resLFC_nonorm005lfclist,
                  set4 = trtpm_nonorm005lfclist), 
             filename = "./venn_DESeq2_nospike_DEGs.png", 
             fill = colors, lwd = 2, lty = 'blank',
             category.names = c("AM-PM in ctrl", "ctrl-chill in AM", 
                                "AM-PM in chill", "ctrl-chill in PM"),
             cat.pos = c(350, 10, 0, 0),
             height = 9, width = 9, units = 'in')

