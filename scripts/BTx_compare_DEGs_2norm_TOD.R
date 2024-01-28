##### scripts for comparing DE genes from DESeq2 analysis
##### Author: Kanjana Laosuntisuk
##### Date created: May 16, 2020
##### Last modified: June 16, 2020

# clear workspace
rm(list = ls())
library(tidyverse)
library(gplots)
library(VennDiagram)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)
setwd("~/sorghum_spikein")

###########################################################
##### DESeq2: controlam VS controlpm ---------------------------------------------
###########################################################
# load dds from nudelta normalization
nudelta_genedds <- readRDS(file = "RData/DESeq2_nudelta_genedds.RData")
# load dds from nospike normalization
nospike_genedds <- readRDS(file = "RData/DESeq2_nospike_genedds.RData")

### get the normalized counts from nospike genedds ----
nospike_normalizedcounts <- DESeq2::counts(nospike_genedds, normalized=TRUE)
nospike_normalizedcounts <- log2(nospike_normalizedcounts + 1)
nospike_normalizedcounts <- as.data.frame(nospike_normalizedcounts)
# select only control data
nospike_normalizedcounts %>%
  dplyr::select(control_am1, control_am2, control_am3, control_am4, 
                control_pm1, control_pm2, control_pm3, control_pm4) -> nospike_normalizedcounts_ctrl

### get the normalized counts from nudelta genedds
nudelta_normalizedcounts <- DESeq2::counts(nudelta_genedds, normalized=TRUE)
#geomean = 37.89637
#nudelta_normalizedcounts <- nudelta_normalizedcounts/geomean
nudelta_normalizedcounts <- log2(nudelta_normalizedcounts + 1)
nudelta_normalizedcounts <- as.data.frame(nudelta_normalizedcounts)
# select only control data
nudelta_normalizedcounts %>%
  dplyr::select(control_am1, control_am2, control_am3, control_am4, 
                control_pm1, control_pm2, control_pm3, control_pm4) -> nudelta_normalizedcounts_ctrl

### combine control dataset
x <- nospike_normalizedcounts_ctrl
names(x)[names(x) == "control_am1"] <- "MedianofRatio_control_am1"
names(x)[names(x) == "control_am2"] <- "MedianofRatio_control_am2"
names(x)[names(x) == "control_am3"] <- "MedianofRatio_control_am3"
names(x)[names(x) == "control_am4"] <- "MedianofRatio_control_am4"
names(x)[names(x) == "control_pm1"] <- "MedianofRatio_control_pm1"
names(x)[names(x) == "control_pm2"] <- "MedianofRatio_control_pm2"
names(x)[names(x) == "control_pm3"] <- "MedianofRatio_control_pm3"
names(x)[names(x) == "control_pm4"] <- "MedianofRatio_control_pm4"
x %>% rownames_to_column(var = "geneID") -> x
y <- nudelta_normalizedcounts_ctrl
names(y)[names(y) == "control_am1"] <- "NuDelta_control_am1"
names(y)[names(y) == "control_am2"] <- "NuDelta_control_am2"
names(y)[names(y) == "control_am3"] <- "NuDelta_control_am3"
names(y)[names(y) == "control_am4"] <- "NuDelta_control_am4"
names(y)[names(y) == "control_pm1"] <- "NuDelta_control_pm1"
names(y)[names(y) == "control_pm2"] <- "NuDelta_control_pm2"
names(y)[names(y) == "control_pm3"] <- "NuDelta_control_pm3"
names(y)[names(y) == "control_pm4"] <- "NuDelta_control_pm4"
y %>% rownames_to_column(var = "geneID") -> y
xy <- inner_join(x, y)

###########################################################
### load AM upregulated gene list----
# load nospike DEGs
nospike_AMup <- read.delim("output/DEGs_DESeq2_nospike_ctrlampm_AM.txt", header = FALSE)
nospike_AMuplist <- nospike_AMup$V1
nudelta_AMup <- read.delim("output/DEGs_DESeq2_nudelta_ctrlampm_AM.txt", header = FALSE)
nudelta_AMuplist <- nudelta_AMup$V1

# make a venn diagram
gplots::venn(list(Median_of_Ratio = nospike_AMuplist,
                  NuDelta = nudelta_AMuplist))

# genes that DE in both method
AMup_both <- intersect(nospike_AMuplist, nudelta_AMuplist) #genes = 2537
# get DEGs only in nudelta only
AMup_nudelta_only <- dplyr::setdiff(nudelta_AMuplist, nospike_AMuplist) #genes = 0
# get DEGs only in median of ratio only
AMup_nospike_only <- dplyr::setdiff(nospike_AMuplist, nudelta_AMuplist) #genes = 1397

### make a heatmap of DEGs in both method ----
AMup_both_xy <- xy[xy$geneID %in% AMup_both, ]
AMup_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> AMup_both_xy2
my_sample_col <- data.frame(condition = rep(c("control_AM", "control_PM", 
                                              "control_AM", "control_PM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(AMup_both_xy2)
my_colour = list(
  condition = c(control_AM = "#F8766D", control_PM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))
pheatmap(as.matrix(AMup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour,
         fontsize = 20)
# save plot
png(filename="pheatmap_DEG_ctrlampm_AMup_both.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(AMup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour,
         fontsize = 20)
dev.off()

### make a heatmap of DEGs in nospike only ----
AMup_nospike_only_xy <- xy[xy$geneID %in% AMup_nospike_only, ]
AMup_nospike_only_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> AMup_nospike_only_xy2
my_sample_col <- data.frame(condition = rep(c("control_AM", "control_PM", 
                                              "control_AM", "control_PM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(AMup_nospike_only_xy2)
my_colour = list(
  condition = c(control_AM = "#F8766D", control_PM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))
pheatmap(as.matrix(AMup_nospike_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour,
         fontsize = 20)
# save plot
png(filename="pheatmap_DEG_ctrlampm_AMup_nospike.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(AMup_nospike_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour, fontsize = 20)
dev.off()

### check gene annotation ----
Sb_annotation <- read.delim("Sbicolor_454_v3.1.1.annotation_info.txt", header = TRUE)
Sb_annotation %>%
  dplyr::select(locusName, transcriptName, peptideName, Pfam,
                Panther, KOG, ec, KO, GO) -> Sb_annotation_Sb
# convert factor to character & replace blank with NA
Sb_annotation_Sb %>% 
  mutate_if(is.factor, as.character) %>%
  mutate_all(na_if, "") -> Sb_annotation_Sb2

# check NAs
map(Sb_annotation_Sb2, ~sum(is.na(.)))
final[complete.cases(final[ , 5:6]),]
Sb_annotation_Sb3 <- Sb_annotation_Sb2[complete.cases(Sb_annotation_Sb2[, c("Panther", "GO")]), ]
# dim # 21922 x 9
map(Sb_annotation_Sb3, ~sum(is.na(.)))
# get locusName that have both Panther and GO IDs
Sb_annotation_Sb3 %>% 
  dplyr::select(locusName) %>%
  distinct() -> Sb_annotation_Sb_name #dim 15225 x 1
Sb_name_AMup_both <- Sb_annotation_Sb_name[Sb_annotation_Sb_name$locusName %in% AMup_both,] # genes = 1364
Sb_name_nospike_only <- Sb_annotation_Sb_name[Sb_annotation_Sb_name$locusName %in% AMup_nospike_only,] # genes = 730


###########################################################
### load PM upregulated gene list----
# load nospike DEGs
nospike_PMup <- read.delim("output/DEGs_DESeq2_nospike_ctrlampm_PM.txt", header = FALSE)
nospike_PMuplist <- nospike_PMup$V1
nudelta_PMup <- read.delim("output/DEGs_DESeq2_nudelta_ctrlampm_PM.txt", header = FALSE)
nudelta_PMuplist <- nudelta_PMup$V1

# make a venn diagram
gplots::venn(list(Median_of_Ratio = nospike_PMuplist,
                  NuDelta = nudelta_PMuplist))

# genes that DE in both method
PMup_both <- intersect(nospike_PMuplist, nudelta_PMuplist) #genes = 4182
# get DEGs only in nudelta only
PMup_nudelta_only <- dplyr::setdiff(nudelta_PMuplist, nospike_PMuplist) #genes = 2380
# get DEGs only in median of ratio only
PMup_nospike_only <- dplyr::setdiff(nospike_PMuplist, nudelta_PMuplist) #genes = 0

### make a heatmap of DEGs in both method ----
PMup_both_xy <- xy[xy$geneID %in% PMup_both, ]
PMup_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> PMup_both_xy2
my_sample_col <- data.frame(condition = rep(c("control_AM", "control_PM", 
                                              "control_AM", "control_PM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(PMup_both_xy2)
my_colour = list(
  condition = c(control_AM = "#F8766D", control_PM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))
pheatmap(as.matrix(PMup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
# save plot
png(filename="pheatmap_DEG_ctrlampm_PMup_both.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(PMup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
dev.off()

### make a heatmap of DEGs in nudelta only ----
PMup_nudelta_only_xy <- xy[xy$geneID %in% PMup_nudelta_only, ]
PMup_nudelta_only_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> PMup_nudelta_only_xy2
my_sample_col <- data.frame(condition = rep(c("control_AM", "control_PM", 
                                              "control_AM", "control_PM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(PMup_nudelta_only_xy2)
my_colour = list(
  condition = c(control_AM = "#F8766D", control_PM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))
pheatmap(as.matrix(PMup_nudelta_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
# save plot
png(filename="pheatmap_DEG_ctrlampm_PMup_nudelta.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(PMup_nudelta_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
dev.off()

### check gene annotation ----
Sb_name_PMup_both <- Sb_annotation_Sb_name[Sb_annotation_Sb_name$locusName %in% PMup_both,] # genes = 2257
Sb_name_nudelta_only <- Sb_annotation_Sb_name[Sb_annotation_Sb_name$locusName %in% PMup_nudelta_only,] # genes = 1252

### total DEGs ----
# make a venn diagram to compare DE genes
colors <- brewer.pal(4, "Pastel1")
colors <- c("#FBB4AE", "#B3CDE3", "#FBB4AE", "#B3CDE3")
venn.diagram(list(Medianof_Ratio_AM_up = nospike_AMuplist,
                  Median_of_Ratio_PM_up = nospike_PMuplist,
                  SpikeIn_AM_up = nudelta_AMuplist,
                  SpikeIn_PM_up = nudelta_PMuplist), 
             filename = "./venn_DESeq2_2norm_ctrlampm.png", 
             fill = colors, lwd = 2, lty = 'blank', cex = 2, cat.cex = 0.5, margin = 0.3,
             main = "Control_AM vs Control_PM", disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(340, 20, 0, 0), main.cex = 2)
# compare DEGs only in AM 
color2 <- c("#FBB4AE", "#B3CDE3")
venn.diagram(list(`Median of Ratio` = nospike_AMuplist,
                  `Spike-in` = nudelta_AMuplist), 
             filename = "./venn_DESeq2_2norm_ctrlampm_AMup.png", 
             fill = color2, lwd = 2, lty = 'blank', cex = 1.8, cat.cex = 1.5, margin = 0.3,
             main = "Control_AM vs Control_PM", sub = "AM up-regulated genes",
             disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(0, 300), main.cex = 2, sub.cex = 1.5)
# compare DEGs only in PM 
color2 <- c("#FBB4AE", "#B3CDE3")
venn.diagram(list(`Median of Ratio` = nospike_PMuplist,
                  `Spike-in` = nudelta_PMuplist), 
             filename = "./venn_DESeq2_2norm_ctrlampm_PMup.png", 
             fill = color2, lwd = 2, lty = 'blank', cex = 1.8, cat.cex = 1.5, margin = 0.3,
             main = "Control_AM vs Control_PM", sub = "PM up-regulated genes",
             disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(300, 0), main.cex = 2, sub.cex = 1.5)

### make a heatmap of DEGs (AM-up + PM-up) in both method ----
DEG_both_xy <- rbind(AMup_both_xy, PMup_both_xy)
DEG_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> DEG_both_xy2
my_sample_col <- data.frame(condition = rep(c("control_AM", "control_PM", 
                                              "control_AM", "control_PM"), c(4,4,4,4)),
                            method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)))
rownames(my_sample_col) <- colnames(DEG_both_xy2)
pheatmap(as.matrix(DEG_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE)
# save plot
png(filename="pheatmap_DEG_ctrlampm_both.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(DEG_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE)
dev.off()

##############################################################
##############################################################
##### Volcano plots ------------------------------------------
### load logFC data
nospike_logFC <- readRDS(file = "RData/DEGs_deseq_full_ctrlampm_full.RData")
nudelta_logFC <- readRDS(file = "RData/DEGs_nudelta_full_ctrlampm_full.RData")

# load nospike DEGs
nospike_DEG <- readRDS(file = "RData/DEGs_deseq_full_ctrlampm.RData")
nospike_DEG_list <- rownames(nospike_DEG)
# load nudelta DEGs
nudelta_DEG <- readRDS(file = "RData/DEGs_nudelta_full_ctrlampm.RData")
nudelta_DEG_list <- rownames(nudelta_DEG)

# DEGs in both methods
DEG_ctrlampm_both <- intersect(nospike_DEG_list, nudelta_DEG_list) #genes = 6719
# DEGs in nospike only
DEG_ctrlampm_nospike <- dplyr::setdiff(nospike_DEG_list, nudelta_DEG_list) #genes = 1397
# DEGs in nudelta only
DEG_ctrlampm_nudelta <- dplyr::setdiff(nudelta_DEG_list, nospike_DEG_list) #genes = 2380

### DEGs in both methods ----
# make a volcano plot of median of ratio
DEG_ctrlampm_both_nospike_logFC <- nospike_logFC[rownames(nospike_logFC) %in% DEG_ctrlampm_both, ]
# make a volcano plot
volcano1 <- EnhancedVolcano(DEG_ctrlampm_both_nospike_logFC,
                            lab = rownames(DEG_ctrlampm_both_nospike_logFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs in the Median of Ratio method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_ctrlampm_DEGboth_nospike.png", width=8.5, height=11, units="in", res=300)
volcano1 
dev.off()

# make a volcano plot of nudelta
DEG_ctrlampm_both_nudelta_logFC <- nudelta_logFC[rownames(nudelta_logFC) %in% DEG_ctrlampm_both, ]
# make a volcano plot
volcano2 <- EnhancedVolcano(DEG_ctrlampm_both_nudelta_logFC,
                            lab = rownames(DEG_ctrlampm_both_nudelta_logFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs in the NuDelta method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_ctrlampm_DEGboth_nudelta.png", width=8.5, height=11, units="in", res=300)
volcano2 
dev.off()

### DEGs in nospike method only ----
# make a volcano plot of median of ratio
DEG_ctrlampm_onlynospike_nospikelogFC <- nospike_logFC[rownames(nospike_logFC) %in% DEG_ctrlampm_nospike, ]
# make a volcano plot
volcano3 <- EnhancedVolcano(DEG_ctrlampm_onlynospike_nospikelogFC,
                            lab = rownames(DEG_ctrlampm_onlynospike_nospikelogFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs only in the Median of Ratio method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_ctrlampm_DEGnospike_nospike.png", width=8.5, height=11, units="in", res=300)
volcano3
dev.off()
# make a volcano plot of nudelta
DEG_ctrlampm_onlynospike_nudeltalogFC <- nudelta_logFC[rownames(nudelta_logFC) %in% DEG_ctrlampm_nospike, ]
# make a volcano plot
volcano4 <- EnhancedVolcano(DEG_ctrlampm_onlynospike_nudeltalogFC,
                            lab = rownames(DEG_ctrlampm_onlynospike_nudeltalogFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs only in the Median of Ratio method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_ctrlampm_DEGnospike_nudelta.png", width=7, height=7, units="in", res=300)
volcano4
dev.off()

# get no. of genes with p-value < 0.05, logFC >= 0.5
onlynospike_nudeltalogFC <- as.data.frame(DEG_ctrlampm_onlynospike_nudeltalogFC)
onlynospike_nudeltalogFC %>%
  rownames_to_column(var = "geneID") %>%
  filter(padj >= 0.05) -> nudelta_lowpadj #dim 1249 x 6
onlynospike_nudeltalogFC %>%
  rownames_to_column(var = "geneID") %>%
  filter(log2FoldChange >= -0.5) -> nudelta_lowFC #dim 1052 x 6
# make a venn diagram
gplots::venn(list(low_padj = nudelta_lowpadj$geneID,
                  low_logFC = nudelta_lowFC$geneID))

####################################
### DEGs in nudelta method only ----
# make a volcano plot of median of ratio
DEG_ctrlampm_onlynudelta_nospikelogFC <- nospike_logFC[rownames(nospike_logFC) %in% DEG_ctrlampm_nudelta, ]
# make a volcano plot
volcano5 <- EnhancedVolcano(DEG_ctrlampm_onlynudelta_nospikelogFC,
                            lab = rownames(DEG_ctrlampm_onlynudelta_nospikelogFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs only in the NuDelta method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 5)
# save plot
png(filename="volcano_ctrlampm_DEGnudelta_nospike.png", width=7, height=7, units="in", res=300)
volcano5
dev.off()

# get no. of genes with p-value < 0.05, logFC >= 0.5
onlynudelta_nospikelogFC <- as.data.frame(DEG_ctrlampm_onlynudelta_nospikelogFC)
onlynudelta_nospikelogFC %>%
  rownames_to_column(var = "geneID") %>%
  filter(padj >= 0.05) -> nospike_lowpadj #dim 2191 x 6
onlynudelta_nospikelogFC %>%
  rownames_to_column(var = "geneID") %>%
  filter(log2FoldChange <= 0.5) -> nospike_lowFC #dim 2380 x 6
# make a venn diagram
gplots::venn(list(low_padj = nospike_lowpadj$geneID,
                  low_logFC = nospike_lowFC$geneID))

# make a volcano plot of nudelta
DEG_ctrlampm_onlynudelta_nudeltalogFC <- nudelta_logFC[rownames(nudelta_logFC) %in% DEG_ctrlampm_nudelta, ]
# make a volcano plot
volcano6 <- EnhancedVolcano(DEG_ctrlampm_onlynudelta_nudeltalogFC,
                            lab = rownames(DEG_ctrlampm_onlynudelta_nudeltalogFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs only in the NuDelta method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_ctrlampm_DEGnudelta_nudelta.png", width=8.5, height=11, units="in", res=300)
volcano6
dev.off()


###########################################################
###########################################################
##### DESeq2: chillingam VS chillingpm ---------------------------------------------
###########################################################
# load dds from nudelta normalization
nudelta_genedds <- readRDS(file = "RData/DESeq2_nudelta_genedds.RData")
# load dds from nospike normalization
nospike_genedds <- readRDS(file = "RData/DESeq2_nospike_genedds.RData")

### get the normalized counts from nospike genedds ----
nospike_normalizedcounts <- DESeq2::counts(nospike_genedds, normalized=TRUE)
nospike_normalizedcounts <- log2(nospike_normalizedcounts + 1)
nospike_normalizedcounts <- as.data.frame(nospike_normalizedcounts)
# select only chilling data
nospike_normalizedcounts %>%
  dplyr::select(chilling_am1, chilling_am2, chilling_am3, chilling_am4, 
                chilling_pm1, chilling_pm2, chilling_pm3, chilling_pm4) -> nospike_normalizedcounts_cold

### get the normalized counts from nudelta genedds
nudelta_normalizedcounts <- DESeq2::counts(nudelta_genedds, normalized=TRUE)
#geomean = 37.89637
#nudelta_normalizedcounts <- nudelta_normalizedcounts/geomean
nudelta_normalizedcounts <- log2(nudelta_normalizedcounts + 1)
nudelta_normalizedcounts <- as.data.frame(nudelta_normalizedcounts)
# select only control data
nudelta_normalizedcounts %>%
  dplyr::select(chilling_am1, chilling_am2, chilling_am3, chilling_am4, 
                chilling_pm1, chilling_pm2, chilling_pm3, chilling_pm4) -> nudelta_normalizedcounts_cold

### combine control dataset
x <- nospike_normalizedcounts_cold
names(x)[names(x) == "chilling_am1"] <- "MedianofRatio_chilling_am1"
names(x)[names(x) == "chilling_am2"] <- "MedianofRatio_chilling_am2"
names(x)[names(x) == "chilling_am3"] <- "MedianofRatio_chilling_am3"
names(x)[names(x) == "chilling_am4"] <- "MedianofRatio_chilling_am4"
names(x)[names(x) == "chilling_pm1"] <- "MedianofRatio_chilling_pm1"
names(x)[names(x) == "chilling_pm2"] <- "MedianofRatio_chilling_pm2"
names(x)[names(x) == "chilling_pm3"] <- "MedianofRatio_chilling_pm3"
names(x)[names(x) == "chilling_pm4"] <- "MedianofRatio_chilling_pm4"
x %>% rownames_to_column(var = "geneID") -> x
y <- nudelta_normalizedcounts_cold
names(y)[names(y) == "chilling_am1"] <- "NuDelta_chilling_am1"
names(y)[names(y) == "chilling_am2"] <- "NuDelta_chilling_am2"
names(y)[names(y) == "chilling_am3"] <- "NuDelta_chilling_am3"
names(y)[names(y) == "chilling_am4"] <- "NuDelta_chilling_am4"
names(y)[names(y) == "chilling_pm1"] <- "NuDelta_chilling_pm1"
names(y)[names(y) == "chilling_pm2"] <- "NuDelta_chilling_pm2"
names(y)[names(y) == "chilling_pm3"] <- "NuDelta_chilling_pm3"
names(y)[names(y) == "chilling_pm4"] <- "NuDelta_chilling_pm4"
y %>% rownames_to_column(var = "geneID") -> y
xy <- inner_join(x, y)

###########################################################
### load AM upregulated gene list----
# load nospike DEGs
nospike_AMup <- read.delim("output/DEGs_DESeq2_nospike_coldampm_AM.txt", header = FALSE)
nospike_AMuplist <- nospike_AMup$V1
nudelta_AMup <- read.delim("output/DEGs_DESeq2_nudelta_coldampm_AM.txt", header = FALSE)
nudelta_AMuplist <- nudelta_AMup$V1

# make a venn diagram
gplots::venn(list(Median_of_Ratio = nospike_AMuplist,
                  NuDelta = nudelta_AMuplist))

# genes that DE in both method
AMup_both <- intersect(nospike_AMuplist, nudelta_AMuplist) #genes = 1136
# get DEGs only in nudelta only
AMup_nudelta_only <- dplyr::setdiff(nudelta_AMuplist, nospike_AMuplist) #genes = 1
# get DEGs only in median of ratio only
AMup_nospike_only <- dplyr::setdiff(nospike_AMuplist, nudelta_AMuplist) #genes = 1397

### make a heatmap of DEGs in both method ----
AMup_both_xy <- xy[xy$geneID %in% AMup_both, ]
AMup_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> AMup_both_xy2
my_sample_col <- data.frame(condition = rep(c("chilling_AM", "chilling_PM", 
                                              "chilling_AM", "chilling_PM"), c(4,4,4,4)),
                            method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)))
rownames(my_sample_col) <- colnames(AMup_both_xy2)
my_colour = list(
  condition = c(chilling_AM = "#F8766D", chilling_PM = "#00BFC4"),
  method = c(Median_of_Ratio = "#00BE67", NuDelta = "#C77CFF"))
pheatmap(as.matrix(AMup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
# save plot
png(filename="pheatmap_DEG_coldampm_AMup_both.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(AMup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
dev.off()

### make a heatmap of DEGs in nospike only ----
AMup_nospike_only_xy <- xy[xy$geneID %in% AMup_nospike_only, ]
AMup_nospike_only_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> AMup_nospike_only_xy2
my_sample_col <- data.frame(condition = rep(c("chilling_AM", "chilling_PM", 
                                              "chilling_AM", "chilling_PM"), c(4,4,4,4)),
                            method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)))
rownames(my_sample_col) <- colnames(AMup_nospike_only_xy2)
my_colour = list(
  condition = c(chilling_AM = "#F8766D", chilling_PM = "#00BFC4"),
  method = c(Median_of_Ratio = "#00BE67", NuDelta = "#C77CFF"))
pheatmap(as.matrix(AMup_nospike_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
# save plot
png(filename="pheatmap_DEG_coldampm_AMup_nospike.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(AMup_nospike_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
dev.off()

### check gene annotation ----
Sb_annotation <- read.delim("Sbicolor_454_v3.1.1.annotation_info.txt", header = TRUE)
Sb_annotation %>%
  dplyr::select(locusName, transcriptName, peptideName, Pfam,
                Panther, KOG, ec, KO, GO) -> Sb_annotation_Sb
# convert factor to character & replace blank with NA
Sb_annotation_Sb %>% 
  mutate_if(is.factor, as.character) %>%
  mutate_all(na_if, "") -> Sb_annotation_Sb2

# check NAs
map(Sb_annotation_Sb2, ~sum(is.na(.)))
final[complete.cases(final[ , 5:6]),]
Sb_annotation_Sb3 <- Sb_annotation_Sb2[complete.cases(Sb_annotation_Sb2[, c("Panther", "GO")]), ]
# dim # 21922 x 9
map(Sb_annotation_Sb3, ~sum(is.na(.)))
# get locusName that have both Panther and GO IDs
Sb_annotation_Sb3 %>% 
  dplyr::select(locusName) %>%
  distinct() -> Sb_annotation_Sb_name #dim 15225 x 1
Sb_name_AMup_both <- Sb_annotation_Sb_name[Sb_annotation_Sb_name$locusName %in% AMup_both,] # genes = 582
Sb_name_nospike_only <- Sb_annotation_Sb_name[Sb_annotation_Sb_name$locusName %in% AMup_nospike_only,] # genes = 1122
Sb_name_nudelta_only <- Sb_annotation_Sb_name[Sb_annotation_Sb_name$locusName %in% AMup_nudelta_only,] # genes = 0


###########################################################
### load PM upregulated gene list----
# load nospike DEGs
nospike_PMup <- read.delim("output/DEGs_DESeq2_nospike_coldampm_PM.txt", header = FALSE)
nospike_PMuplist <- nospike_PMup$V1
nudelta_PMup <- read.delim("output/DEGs_DESeq2_nudelta_coldampm_PM.txt", header = FALSE)
nudelta_PMuplist <- nudelta_PMup$V1

# make a venn diagram
gplots::venn(list(Median_of_Ratio = nospike_PMuplist,
                  NuDelta = nudelta_PMuplist))

# genes that DE in both method
PMup_both <- intersect(nospike_PMuplist, nudelta_PMuplist) #genes = 3329
# get DEGs only in nudelta only
PMup_nudelta_only <- dplyr::setdiff(nudelta_PMuplist, nospike_PMuplist) #genes = 6141
# get DEGs only in median of ratio only
PMup_nospike_only <- dplyr::setdiff(nospike_PMuplist, nudelta_PMuplist) #genes = 0

### make a heatmap of DEGs in both method ----
PMup_both_xy <- xy[xy$geneID %in% PMup_both, ]
PMup_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> PMup_both_xy2
my_sample_col <- data.frame(condition = rep(c("chilling_AM", "chilling_PM", 
                                              "chilling_AM", "chilling_PM"), c(4,4,4,4)),
                            method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)))
rownames(my_sample_col) <- colnames(PMup_both_xy2)
my_colour = list(
  condition = c(chilling_AM = "#F8766D", chilling_PM = "#00BFC4"),
  method = c(Median_of_Ratio = "#00BE67", NuDelta = "#C77CFF"))
pheatmap(as.matrix(PMup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
# save plot
png(filename="pheatmap_DEG_coldampm_PMup_both.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(PMup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
dev.off()

### make a heatmap of DEGs in nudelta only ----
PMup_nudelta_only_xy <- xy[xy$geneID %in% PMup_nudelta_only, ]
PMup_nudelta_only_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> PMup_nudelta_only_xy2
my_sample_col <- data.frame(condition = rep(c("chilling_AM", "chilling_PM", 
                                              "chilling_AM", "chilling_PM"), c(4,4,4,4)),
                            method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)))
rownames(my_sample_col) <- colnames(PMup_nudelta_only_xy2)
my_colour = list(
  condition = c(chilling_AM = "#F8766D", chilling_PM = "#00BFC4"),
  method = c(Median_of_Ratio = "#00BE67", NuDelta = "#C77CFF"))
pheatmap(as.matrix(PMup_nudelta_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
# save plot
png(filename="pheatmap_DEG_coldampm_PMup_nudelta.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(PMup_nudelta_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour)
dev.off()

### check gene annotation ----
Sb_name_PMup_both <- Sb_annotation_Sb_name[Sb_annotation_Sb_name$locusName %in% PMup_both,] # genes = 1685
Sb_name_PMnudelta_only <- Sb_annotation_Sb_name[Sb_annotation_Sb_name$locusName %in% PMup_nudelta_only,] # genes = 3278

### total DEGs ----
# make a venn diagram to compare DE genes
colors <- brewer.pal(4, "Pastel1")
colors <- c("#FBB4AE", "#B3CDE3", "#FBB4AE", "#B3CDE3")
venn.diagram(list(Median_of_Ratio_AM_up = nospike_AMuplist,
                  Median_of_Ratio_PM_up = nospike_PMuplist,
                  SpikeIn_AM_up = nudelta_AMuplist,
                  SpikeIn_PM_up = nudelta_PMuplist), 
             filename = "./venn_DESeq2_2norm_coldampm.png", 
             fill = colors, lwd = 2, lty = 'blank', cex = 2, cat.cex = 0.5, margin = 0.3,
             main = "Chilling_AM vs Chilling_PM", disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(340, 20, 0, 0), main.cex = 2)
# compare DEGs only in AM 
color2 <- c("#FBB4AE", "#B3CDE3")
venn.diagram(list(`Median of Ratio` = nospike_AMuplist,
                  `Spike-in` = nudelta_AMuplist), 
             filename = "./venn_DESeq2_2norm_coldampm_AMup.png", 
             fill = color2, lwd = 2, lty = 'blank', cex = 1.8, cat.cex = 1.5, margin = 0.3,
             main = "Chilling_AM vs Chilling_PM", sub = "AM up-regulated genes",
             disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(300,0), main.cex = 2, sub.cex = 1.5)
# compare DEGs only in PM 
color2 <- c("#FBB4AE", "#B3CDE3")
venn.diagram(list(`Median of Ratio` = nospike_PMuplist,
                  `Spike-in` = nudelta_PMuplist), 
             filename = "./venn_DESeq2_2norm_coldampm_PMup.png", 
             fill = color2, lwd = 2, lty = 'blank', cex = 1.8, cat.cex = 1.5, margin = 0.3,
             main = "Chilling_AM vs Chilling_PM", sub = "PM up-regulated genes",
             disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(300,0), main.cex = 2, sub.cex = 1.5)

### make a heatmap of DEGs (AM-up + PM-up) in both method ----
cold_both_xy <- rbind(AMup_both_xy, PMup_both_xy)
cold_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> cold_both_xy2
my_sample_col <- data.frame(condition = rep(c("chilling_AM", "chilling_PM", 
                                              "chilling_AM", "chilling_PM"), c(4,4,4,4)),
                            method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)))
rownames(my_sample_col) <- colnames(cold_both_xy2)
pheatmap(as.matrix(cold_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE)
# save plot
png(filename="pheatmap_DEG_coldampm_both.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(cold_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE)
dev.off()

### make a heatmap of DEGs (AM-up + PM-up) in nudelta method ----
AMup_nudelta_only_xy <- xy[xy$geneID %in% AMup_nudelta_only, ]
cold_nudelta_xy <- rbind(AMup_nudelta_only_xy, PMup_nudelta_only_xy)
cold_nudelta_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> cold_nudelta_xy2
my_sample_col <- data.frame(condition = rep(c("chilling_AM", "chilling_PM", 
                                              "chilling_AM", "chilling_PM"), c(4,4,4,4)),
                            method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)))
rownames(my_sample_col) <- colnames(cold_nudelta_xy2)
pheatmap(as.matrix(cold_nudelta_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE)
# save plot
png(filename="pheatmap_DEG_coldampm_nudelta.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(cold_nudelta_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE)
dev.off()

##############################################################
##############################################################
##### Volcano plots ------------------------------------------
### load logFC data
nospike_logFC <- readRDS(file = "RData/DEGs_deseq_full_coldampm_full.RData")
nudelta_logFC <- readRDS(file = "RData/DEGs_nudelta_full_coldampm_full.RData")

# load nospike DEGs
nospike_DEG <- readRDS(file = "RData/DEGs_deseq_full_coldampm.RData")
nospike_DEG_list <- rownames(nospike_DEG)
# load nudelta DEGs
nudelta_DEG <- readRDS(file = "RData/DEGs_nudelta_full_coldampm.RData")
nudelta_DEG_list <- rownames(nudelta_DEG)

# DEGs in both methods
DEG_coldampm_both <- intersect(nospike_DEG_list, nudelta_DEG_list) #genes = 4465
# DEGs in nospike only
DEG_coldampm_nospike <- dplyr::setdiff(nospike_DEG_list, nudelta_DEG_list) #genes = 2012
# DEGs in nudelta only
DEG_coldampm_nudelta <- dplyr::setdiff(nudelta_DEG_list, nospike_DEG_list) #genes = 6142
DEG_coldampm_nudelta <- dplyr::setdiff(nudelta_DEG_list, nospike_DEG_list) #genes = 6142

### DEGs in both methods ----
# make a volcano plot of median of ratio
DEG_coldampm_both_nospike_logFC <- nospike_logFC[rownames(nospike_logFC) %in% DEG_coldampm_both, ]
# make a volcano plot
volcano1 <- EnhancedVolcano(DEG_coldampm_both_nospike_logFC,
                            lab = rownames(DEG_coldampm_both_nospike_logFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs in the Median of Ratio method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_coldampm_DEGboth_nospike.png", width=8.5, height=11, units="in", res=300)
volcano1 
dev.off()

# make a volcano plot of nudelta
DEG_coldampm_both_nudelta_logFC <- nudelta_logFC[rownames(nudelta_logFC) %in% DEG_coldampm_both, ]
# make a volcano plot
volcano2 <- EnhancedVolcano(DEG_coldampm_both_nudelta_logFC,
                            lab = rownames(DEG_coldampm_both_nudelta_logFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs in the NuDelta method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_coldampm_DEGboth_nudelta.png", width=8.5, height=11, units="in", res=300)
volcano2 
dev.off()

### DEGs in nospike method only ----
# make a volcano plot of median of ratio
DEG_coldampm_onlynospike_nospikelogFC <- nospike_logFC[rownames(nospike_logFC) %in% DEG_coldampm_nospike, ]
# make a volcano plot
volcano3 <- EnhancedVolcano(DEG_coldampm_onlynospike_nospikelogFC,
                            lab = rownames(DEG_coldampm_onlynospike_nospikelogFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs only in the Median of Ratio method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_coldampm_DEGnospike_nospike.png", width=8.5, height=11, units="in", res=300)
volcano3
dev.off()
# make a volcano plot of nudelta
DEG_coldampm_onlynospike_nudeltalogFC <- nudelta_logFC[rownames(nudelta_logFC) %in% DEG_coldampm_nospike, ]
# make a volcano plot
volcano4 <- EnhancedVolcano(DEG_coldampm_onlynospike_nudeltalogFC,
                            lab = rownames(DEG_coldampm_onlynospike_nudeltalogFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs only in the Median of Ratio method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_coldampm_DEGnospike_nudelta.png", width=7, height=7, units="in", res=300)
volcano4
dev.off()

# get no. of genes with p-value > 0.05, |logFC| <= 0.5
onlynospike_nudeltalogFC <- as.data.frame(DEG_coldampm_onlynospike_nudeltalogFC)
onlynospike_nudeltalogFC %>%
  rownames_to_column(var = "geneID") %>%
  filter(padj >= 0.05) -> nudelta_lowpadj #dim 1928 x 6
onlynospike_nudeltalogFC %>%
  rownames_to_column(var = "geneID") %>%
  filter(log2FoldChange >= -0.5) -> nudelta_lowFC #dim 1738 x 6
# make a venn diagram
gplots::venn(list(low_padj = nudelta_lowpadj$geneID,
                  low_logFC = nudelta_lowFC$geneID))

####################################
### DEGs in nudelta method only ----
# make a volcano plot of median of ratio
DEG_coldampm_onlynudelta_nospikelogFC <- nospike_logFC[rownames(nospike_logFC) %in% DEG_coldampm_nudelta, ]
DEG_coldampm_onlynudelta_nospikelogFC2 <- DEG_coldampm_onlynudelta_nospikelogFC[rownames(DEG_coldampm_onlynudelta_nospikelogFC) != "Sobic.002G092400", ]
# make a volcano plot
volcano5 <- EnhancedVolcano(DEG_coldampm_onlynudelta_nospikelogFC2,
                            lab = rownames(DEG_coldampm_onlynudelta_nospikelogFC2),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs only in the NuDelta method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_coldampm_DEGnudelta_nospike.png", width=7, height=7, units="in", res=300)
volcano5
dev.off()

# get no. of genes with p-value > 0.05, |logFC| <= 0.5
onlynudelta_nospikelogFC <- as.data.frame(DEG_coldampm_onlynudelta_nospikelogFC2)
onlynudelta_nospikelogFC %>%
  rownames_to_column(var = "geneID") %>%
  filter(padj >= 0.05) -> nospike_lowpadj #dim 5919 x 6
onlynudelta_nospikelogFC %>%
  rownames_to_column(var = "geneID") %>%
  filter(log2FoldChange <= 0.5) -> nospike_lowFC #dim 5265 x 6
# make a venn diagram
gplots::venn(list(low_padj = nospike_lowpadj$geneID,
                  low_logFC = nospike_lowFC$geneID))

# make a volcano plot of nudelta
DEG_coldampm_onlynudelta_nudeltalogFC <- nudelta_logFC[rownames(nudelta_logFC) %in% DEG_coldampm_nudelta, ]
# make a volcano plot
volcano6 <- EnhancedVolcano(DEG_coldampm_onlynudelta_nudeltalogFC,
                            lab = rownames(DEG_coldampm_onlynudelta_nudeltalogFC),
                            x = 'log2FoldChange',
                            y = 'padj', xlim = c(-8, 8), 
                            title = 'DEGs only in the NuDelta method',
                            pCutoff = 0.05, FCcutoff = 0.5, pointSize = 1.5, labSize = 2.5)
# save plot
png(filename="volcano_coldampm_DEGnudelta_nudelta.png", width=8.5, height=11, units="in", res=300)
volcano6
dev.off()




