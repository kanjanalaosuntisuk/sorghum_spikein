##### scripts for comparing DEGs from DESeq2 analysis
##### Median of Ratio vs Spike-in normalization

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
##### DESeq2: controlam VS chillingam ---------------------------------------------
###########################################################
# load dds from nudelta normalization
nudelta_genedds <- readRDS(file = "RData/DESeq2_nudelta_genedds.RData")
# load dds from nospike normalization
nospike_genedds <- readRDS(file = "RData/DESeq2_nospike_genedds.RData")

### get the normalized counts from nospike genedds ----
nospike_normalizedcounts <- DESeq2::counts(nospike_genedds, normalized=TRUE)
nospike_normalizedcounts <- log2(nospike_normalizedcounts + 1)
nospike_normalizedcounts <- as.data.frame(nospike_normalizedcounts)
# select only AM data
nospike_normalizedcounts %>%
  dplyr::select(control_am1, control_am2, control_am3, control_am4, 
                chilling_am1, chilling_am2, chilling_am3, chilling_am4) -> nospike_normalizedcounts_AM

### get the normalized counts from nudelta genedds
nudelta_normalizedcounts <- DESeq2::counts(nudelta_genedds, normalized=TRUE)
#geomean = 37.89637
#nudelta_normalizedcounts <- nudelta_normalizedcounts/geomean
nudelta_normalizedcounts <- log2(nudelta_normalizedcounts + 1)
nudelta_normalizedcounts <- as.data.frame(nudelta_normalizedcounts)
# select only control data
nudelta_normalizedcounts %>%
  dplyr::select(control_am1, control_am2, control_am3, control_am4, 
                chilling_am1, chilling_am2, chilling_am3, chilling_am4) -> nudelta_normalizedcounts_AM

### combine control dataset
x <- nospike_normalizedcounts_AM
names(x)[names(x) == "control_am1"] <- "MedianofRatio_control_am1"
names(x)[names(x) == "control_am2"] <- "MedianofRatio_control_am2"
names(x)[names(x) == "control_am3"] <- "MedianofRatio_control_am3"
names(x)[names(x) == "control_am4"] <- "MedianofRatio_control_am4"
names(x)[names(x) == "chilling_am1"] <- "MedianofRatio_chilling_am1"
names(x)[names(x) == "chilling_am2"] <- "MedianofRatio_chilling_am2"
names(x)[names(x) == "chilling_am3"] <- "MedianofRatio_chilling_am3"
names(x)[names(x) == "chilling_am4"] <- "MedianofRatio_chilling_am4"
x %>% rownames_to_column(var = "geneID") -> x
y <- nudelta_normalizedcounts_AM
names(y)[names(y) == "control_am1"] <- "NuDelta_control_am1"
names(y)[names(y) == "control_am2"] <- "NuDelta_control_am2"
names(y)[names(y) == "control_am3"] <- "NuDelta_control_am3"
names(y)[names(y) == "control_am4"] <- "NuDelta_control_am4"
names(y)[names(y) == "chilling_am1"] <- "NuDelta_chilling_am1"
names(y)[names(y) == "chilling_am2"] <- "NuDelta_chilling_am2"
names(y)[names(y) == "chilling_am3"] <- "NuDelta_chilling_am3"
names(y)[names(y) == "chilling_am4"] <- "NuDelta_chilling_am4"
y %>% rownames_to_column(var = "geneID") -> y
xy <- inner_join(x, y)

###########################################################
### load cold upregulated gene list----
# load nospike DEGs
nospike_coldup <- read.delim("output/DEGs_DESeq2_nospike_ctrlcoldam_cold.txt", header = FALSE)
nospike_colduplist <- nospike_coldup$V1
nudelta_coldup <- read.delim("output/DEGs_DESeq2_nudelta_ctrlcoldam_cold.txt", header = FALSE)
nudelta_colduplist <- nudelta_coldup$V1

# make a venn diagram
gplots::venn(list(Median_of_Ratio = nospike_colduplist,
                  NuDelta = nudelta_colduplist))

# genes that DE in both method
coldup_both <- intersect(nospike_colduplist, nudelta_colduplist) #genes = 2876
# get DEGs only in nudelta only
coldup_nudelta_only <- dplyr::setdiff(nudelta_colduplist, nospike_colduplist) #genes = 1
# get DEGs only in median of ratio only
coldup_nospike_only <- dplyr::setdiff(nospike_colduplist, nudelta_colduplist) #genes = 863

### make a heatmap of DEGs in both method ----
coldup_both_xy <- xy[xy$geneID %in% coldup_both, ]
coldup_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> coldup_both_xy2
my_sample_col <- data.frame(condition = rep(c("control_AM", "chilling_AM", 
                                              "control_AM", "chilling_AM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(coldup_both_xy2)
my_colour = list(
  condition = c(control_AM = "#F8766D", chilling_AM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))

# save plot
#png(filename="pheatmap_DEG_ctrlcoldam_coldup_both.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(coldup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour, fontsize = 20)
#dev.off()

### make a heatmap of DEGs in nospike only ----
coldup_nospike_only_xy <- xy[xy$geneID %in% coldup_nospike_only, ]
coldup_nospike_only_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> coldup_nospike_only_xy2
my_sample_col <- data.frame(condition = rep(c("control_AM", "chilling_AM", 
                                              "control_AM", "chilling_AM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(coldup_nospike_only_xy2)
my_colour = list(
  condition = c(control_AM = "#F8766D", chilling_AM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))

# save plot
#png(filename="pheatmap_DEG_ctrlcoldam_coldup_nospike.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(coldup_nospike_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour, fontsize = 20)
#dev.off()


###########################################################
### load ctrl upregulated gene list----
# load nospike DEGs
nospike_ctrlup <- read.delim("output/DEGs_DESeq2_nospike_ctrlcoldam_ctrl.txt", header = FALSE)
nospike_ctrluplist <- nospike_ctrlup$V1
nudelta_ctrlup <- read.delim("output/DEGs_DESeq2_nudelta_ctrlcoldam_ctrl.txt", header = FALSE)
nudelta_ctrluplist <- nudelta_ctrlup$V1

# make a venn diagram
gplots::venn(list(Median_of_Ratio = nospike_ctrluplist,
                  NuDelta = nudelta_ctrluplist))

# genes that DE in both method
ctrlup_both <- intersect(nospike_ctrluplist, nudelta_ctrluplist) #genes = 3411
# get DEGs only in nudelta only
ctrlup_nudelta_only <- dplyr::setdiff(nudelta_ctrluplist, nospike_ctrluplist) #genes = 1160
# get DEGs only in median of ratio only
ctrlup_nospike_only <- dplyr::setdiff(nospike_ctrluplist, nudelta_ctrluplist) #genes = 0

### make a heatmap of DEGs in both method ----
ctrlup_both_xy <- xy[xy$geneID %in% ctrlup_both, ]
ctrlup_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> ctrlup_both_xy2
my_sample_col <- data.frame(condition = rep(c("control_AM", "chilling_AM", 
                                              "control_AM", "chilling_AM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(ctrlup_both_xy2)
my_colour = list(
  condition = c(control_AM = "#F8766D", chilling_AM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))

# save plot
#png(filename="pheatmap_DEG_ctrlcoldam_ctrlup_both.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(ctrlup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour, fontsize = 20)
#dev.off()

### make a heatmap of DEGs in nudelta only ----
ctrlup_nudelta_only_xy <- xy[xy$geneID %in% ctrlup_nudelta_only, ]
ctrlup_nudelta_only_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> ctrlup_nudelta_only_xy2
my_sample_col <- data.frame(condition = rep(c("control_AM", "chilling_AM", 
                                              "control_AM", "chilling_AM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(ctrlup_nudelta_only_xy2)
my_colour = list(
  condition = c(control_AM = "#F8766D", chilling_AM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))

# save plot
#png(filename="pheatmap_DEG_ctrlcoldam_ctrlup_nudelta.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(ctrlup_nudelta_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour, fontsize = 20)
#dev.off()

### total DEGs ----
# make a venn diagram to compare DE genes
colors <- brewer.pal(4, "Pastel1")
colors <- c( "#B3CDE3", "#FBB4AE", "#B3CDE3", "#FBB4AE")
venn.diagram(list(Median_of_Ratio_chilling_up = nospike_colduplist,
                  Median_of_Ratio_chilling_down = nospike_ctrluplist,
                  SpikeIn_chilling_up = nudelta_colduplist,
                  SpikeIn_chilling_down = nudelta_ctrluplist), 
             filename = "./venn_DESeq2_2norm_ctrlcoldam.png", 
             fill = colors, lwd = 2, lty = 'blank', cex = 2, cat.cex = 0.3, margin = 0.5,
             main = "Control_AM vs Chilling_AM", disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(320, 40, 320, 40), main.cex = 2)
# compare DEGs only in cold
color2 <- c("#FBB4AE", "#B3CDE3")
venn.diagram(list(`Median of Ratio` = nospike_colduplist,
                  `Spike-in` = nudelta_colduplist), 
             filename = "./venn_DESeq2_2norm_ctrlcoldam_coldup.png", 
             fill = color2, lwd = 2, lty = 'blank', cex = 1.8, cat.cex = 1.5, margin = 0.5,
             main = "Control_AM vs Chilling_AM", sub = "chilling up-regulated genes",
             disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(290,0), main.cex = 2, sub.cex = 1.5)
# compare DEGs only in ctrl
color2 <- c("#FBB4AE", "#B3CDE3")
venn.diagram(list(`Median of Ratio` = nospike_ctrluplist,
                  `Spike-in` = nudelta_ctrluplist), 
             filename = "./venn_DESeq2_2norm_ctrlcoldam_ctrlup.png", 
             fill = color2, lwd = 2, lty = 'blank', cex = 1.8, cat.cex = 1.5, margin = 0.5,
             main = "Control_AM vs Chilling_AM", sub = "chilling down-regulated genes",
             disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(300,0), main.cex = 2, sub.cex = 1.5)


###########################################################
###########################################################
##### DESeq2: controlpm VS chillingpm ---------------------------------------------
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
  dplyr::select(control_pm1, control_pm2, control_pm3, control_pm4, 
                chilling_pm1, chilling_pm2, chilling_pm3, chilling_pm4) -> nospike_normalizedcounts_pm

### get the normalized counts from nudelta genedds
nudelta_normalizedcounts <- DESeq2::counts(nudelta_genedds, normalized=TRUE)
#geomean = 37.89637
#nudelta_normalizedcounts <- nudelta_normalizedcounts/geomean
nudelta_normalizedcounts <- log2(nudelta_normalizedcounts + 1)
nudelta_normalizedcounts <- as.data.frame(nudelta_normalizedcounts)
# select only control data
nudelta_normalizedcounts %>%
  dplyr::select(control_pm1, control_pm2, control_pm3, control_pm4, 
                chilling_pm1, chilling_pm2, chilling_pm3, chilling_pm4) -> nudelta_normalizedcounts_pm

### combine control dataset
x <- nospike_normalizedcounts_pm
names(x)[names(x) == "control_pm1"] <- "MedianofRatio_control_pm1"
names(x)[names(x) == "control_pm2"] <- "MedianofRatio_control_pm2"
names(x)[names(x) == "control_pm3"] <- "MedianofRatio_control_pm3"
names(x)[names(x) == "control_pm4"] <- "MedianofRatio_control_pm4"
names(x)[names(x) == "chilling_pm1"] <- "MedianofRatio_chilling_pm1"
names(x)[names(x) == "chilling_pm2"] <- "MedianofRatio_chilling_pm2"
names(x)[names(x) == "chilling_pm3"] <- "MedianofRatio_chilling_pm3"
names(x)[names(x) == "chilling_pm4"] <- "MedianofRatio_chilling_pm4"
x %>% rownames_to_column(var = "geneID") -> x
y <- nudelta_normalizedcounts_pm
names(y)[names(y) == "control_pm1"] <- "NuDelta_control_pm1"
names(y)[names(y) == "control_pm2"] <- "NuDelta_control_pm2"
names(y)[names(y) == "control_pm3"] <- "NuDelta_control_pm3"
names(y)[names(y) == "control_pm4"] <- "NuDelta_control_pm4"
names(y)[names(y) == "chilling_pm1"] <- "NuDelta_chilling_pm1"
names(y)[names(y) == "chilling_pm2"] <- "NuDelta_chilling_pm2"
names(y)[names(y) == "chilling_pm3"] <- "NuDelta_chilling_pm3"
names(y)[names(y) == "chilling_pm4"] <- "NuDelta_chilling_pm4"
y %>% rownames_to_column(var = "geneID") -> y
xy <- inner_join(x, y)

###########################################################
### load cold upregulated gene list----
# load nospike DEGs
nospike_coldup <- read.delim("output/DEGs_DESeq2_nospike_ctrlcoldpm_cold.txt", header = FALSE)
nospike_colduplist <- nospike_coldup$V1
nudelta_coldup <- read.delim("output/DEGs_DESeq2_nudelta_ctrlcoldpm_cold.txt", header = FALSE)
nudelta_colduplist <- nudelta_coldup$V1

# make a venn diagram
gplots::venn(list(Median_of_Ratio = nospike_colduplist,
                  NuDelta = nudelta_colduplist))

# genes that DE in both method
coldup_both <- intersect(nospike_colduplist, nudelta_colduplist) #genes = 2611
# get DEGs only in nudelta only
coldup_nudelta_only <- dplyr::setdiff(nudelta_colduplist, nospike_colduplist) #genes = 867
# get DEGs only in median of ratio only
coldup_nospike_only <- dplyr::setdiff(nospike_colduplist, nudelta_colduplist) #genes = 0

### make a heatmap of DEGs in both method ----
coldup_both_xy <- xy[xy$geneID %in% coldup_both, ]
coldup_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> coldup_both_xy2
my_sample_col <- data.frame(condition = rep(c("control_PM", "chilling_PM", 
                                              "control_PM", "chilling_PM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(coldup_both_xy2)
my_colour = list(
  condition = c(control_PM = "#F8766D", chilling_PM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))

# save plot
#png(filename="pheatmap_DEG_ctrlcoldpm_coldup_both.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(coldup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE,
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour, fontsize = 20)
#dev.off()

### make a heatmap of DEGs in nudelta only ----
coldup_nudelta_only_xy <- xy[xy$geneID %in% coldup_nudelta_only, ]
coldup_nudelta_only_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> coldup_nudelta_only_xy2
my_sample_col <- data.frame(condition = rep(c("control_PM", "chilling_PM", 
                                              "control_PM", "chilling_PM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(coldup_nudelta_only_xy2)
my_colour = list(
  condition = c(control_PM = "#F8766D", chilling_PM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))

# save plot
#png(filename="pheatmap_DEG_ctrlcoldpm_coldup_nudelta.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(coldup_nudelta_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour, fontsize = 20)
#dev.off()


###########################################################
### load control upregulated gene list----
# load nospike DEGs
nospike_ctrlup <- read.delim("output/DEGs_DESeq2_nospike_ctrlcoldpm_ctrl.txt", header = FALSE)
nospike_ctrluplist <- nospike_ctrlup$V1
nudelta_ctrlup <- read.delim("output/DEGs_DESeq2_nudelta_ctrlcoldpm_ctrl.txt", header = FALSE)
nudelta_ctrluplist <- nudelta_ctrlup$V1

# make a venn diagram
gplots::venn(list(Median_of_Ratio = nospike_ctrluplist,
                  NuDelta = nudelta_ctrluplist))

# genes that DE in both method
ctrlup_both <- intersect(nospike_ctrluplist, nudelta_ctrluplist) #genes = 1607
# get DEGs only in nudelta only
ctrlup_nudelta_only <- dplyr::setdiff(nudelta_ctrluplist, nospike_ctrluplist) #genes = 0
# get DEGs only in median of ratio only
ctrlup_nospike_only <- dplyr::setdiff(nospike_ctrluplist, nudelta_ctrluplist) #genes = 702

### make a heatmap of DEGs in both method ----
ctrlup_both_xy <- xy[xy$geneID %in% ctrlup_both, ]
ctrlup_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> ctrlup_both_xy2
my_sample_col <- data.frame(condition = rep(c("control_PM", "chilling_PM", 
                                              "control_PM", "chilling_PM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(ctrlup_both_xy2)
my_colour = list(
  condition = c(control_PM = "#F8766D", chilling_PM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))

# save plot
#png(filename="pheatmap_DEG_ctrlcoldpm_ctrlup_both.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(ctrlup_both_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour, fontsize = 20)
#dev.off()

### make a heatmap of DEGs in nospike only ----
ctrlup_nospike_only_xy <- xy[xy$geneID %in% ctrlup_nospike_only, ]
ctrlup_nospike_only_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> ctrlup_nospike_only_xy2
my_sample_col <- data.frame(condition = rep(c("control_PM", "chilling_PM", 
                                              "control_PM", "chilling_PM"), c(4,4,4,4)),
                            method = rep(c("Median of Ratio", "Spike-in"), c(8,8)))
rownames(my_sample_col) <- colnames(ctrlup_nospike_only_xy2)
my_colour = list(
  condition = c(control_PM = "#F8766D", chilling_PM = "#00BFC4"),
  method = c(`Median of Ratio` = "#00BE67", `Spike-in` = "#C77CFF"))

# save plot
#png(filename="pheatmap_DEG_ctrlcoldpm_ctrlup_nospike.png", width=9, height=7, units="in", res=300)
pheatmap(as.matrix(ctrlup_nospike_only_xy2), cluster_rows = TRUE, cluster_cols = FALSE, 
         show_rownames = FALSE, border_color = NA, annotation_col = my_sample_col,
         show_colnames = FALSE, scale="row", annotation_colors = my_colour, fontsize = 20)
#dev.off()

### total DEGs ----
# make a venn diagram to compare DE genes
colors <- brewer.pal(4, "Pastel1")
colors <- c( "#B3CDE3", "#FBB4AE", "#B3CDE3", "#FBB4AE")
venn.diagram(list(Median_of_Ratio_chilling_up = nospike_colduplist,
                  Median_of_Ratio_chilling_down = nospike_ctrluplist,
                  SpikeIn_chilling_up = nudelta_colduplist,
                  SpikeIn_chilling_down = nudelta_ctrluplist), 
             filename = "./venn_DESeq2_2norm_ctrlcoldpm.png", 
             fill = colors, lwd = 2, lty = 'blank', cex = 2, cat.cex = 0.3, margin = 0.5,
             main = "Control_PM vs Chilling_PM", disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(320, 40, 320, 40), main.cex = 2)

# compare DEGs only in cold
color2 <- c("#FBB4AE", "#B3CDE3")
venn.diagram(list(`Median of Ratio` = nospike_colduplist,
                  `Spike-in` = nudelta_colduplist), 
             filename = "./venn_DESeq2_2norm_ctrlcoldpm_coldup.png", 
             fill = color2, lwd = 2, lty = 'blank', cex = 1.8, cat.cex = 1.5, margin = 0.5,
             main = "Control_PM vs Chilling_PM", sub = "chilling up-regulated genes",
             disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(300,0), main.cex = 2, sub.cex = 1.5)

# compare DEGs only in ctrl
color2 <- c("#FBB4AE", "#B3CDE3")
venn.diagram(list(`Median of Ratio` = nospike_ctrluplist,
                  `Spike-in` = nudelta_ctrluplist), 
             filename = "./venn_DESeq2_2norm_ctrlcoldpm_ctrlup.png", 
             fill = color2, lwd = 2, lty = 'blank', cex = 1.8, cat.cex = 1.5, margin = 0.5,
             main = "Control_PM vs Chilling_PM", sub = "chilling down-regulated genes",
             disable.logging=TRUE, main.pos = c(0.5,1.05),
             cat.pos = c(350,300), main.cex = 2, sub.cex = 1.5)






