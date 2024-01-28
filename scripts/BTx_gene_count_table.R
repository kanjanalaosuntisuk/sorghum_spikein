##### scripts for making gene read count table -----

# clear workspace
rm(list = ls())
# load packages
library(tidyverse)
library(RColorBrewer)
library(edgeR)
library(EDASeq)
#setwd("~/sorghum_spikein")

#################################################################
### load gene count table
BTxgene <- readRDS(file = "RData/BTxgene.RData")

#################################################################
##### filter non-expressed genes by filterByExpr function
##### reference: https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
# make a design matrix and DGEList object
conditions <- as.factor(rep(c("controlam", "controlpm", "chillingam", "chillingpm"), each=4))
conditions <- factor(conditions, levels = c("controlam", "controlpm", "chillingam", "chillingpm"))
design <- model.matrix(~ 0 + conditions)
y <- DGEList(BTxgene, group = conditions)

# data transformation of non-filtered gene read counts
cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)
L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
c(L, M)
summary(lcpm)

# filter low-expressed genes
table(rowSums(y$counts==0)==16)
keep.exprs <- filterByExpr(y, design)
y_a <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y_a)
# get filterd dataframe
genefiltered <- y_a$counts
genefilteredList <- rownames(genefiltered)
#write.csv(genefiltered, file = "unnormcounts_BTx.txt")
# save RData
saveRDS(genefiltered, file = "RData/genefiltered.RData")
#saveRDS(genefilteredList, file = "RData/genefilteredList_forTopGO.RData")
#write.table(genefilteredList, file = "Sb_genelist.txt", sep=",", 
#            row.names = FALSE, quote=FALSE, col.names=FALSE)

# make side-by-side RLE and PCA plot
#png(filename="BTxgene_RLE_PCAplot.png", width=7, height=3, units="in", res=300)
par(mfrow = c(1, 2))
par(mar = c(5, 4, 1, 1))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], xaxt = "n", 
        ylab = "Relative log expression", main = "Before normalization", cex.main=1, 
        cex.lab = 1, cex.axis=0.8)
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.8)
plotPCA(set, col=colors[x], cex=0.7, cex.lab = 1, cex.axis = 0.8)
#dev.off()

#####################################################
### sum total filtered reads in each sample
readsum_genefiltered <- colSums(genefiltered)
readsum_genefiltered_df <- as.data.frame(readsum_genefiltered)
readsum_genefiltered_df %>%
  rownames_to_column(var = "rowname") -> readsum_genefiltered_df
readsum_genefiltered_df %>%
  mutate(conditions = as.factor(rep(c("controlam", "controlpm", "chillingam", "chillingpm"), each=4))) -> readsum_genefiltered_df
readsum_genefiltered_df$conditions <- factor(readsum_genefiltered_df$conditions, 
                                             levels = c("controlam", "controlpm", "chillingam", "chillingpm"))
readsum_genefiltered_df$rowname <- factor(readsum_genefiltered_df$rowname, 
                                          levels = readsum_genefiltered_df$rowname[order(readsum_genefiltered_df$conditions)])
