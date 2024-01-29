##### scripts for making plots of selected DEGs from DESeq2 analysis ----
##### Effect of chilling stress on gene expression under AM and PM ----

library(DESeq2)
library(tidyverse)
setwd("~/sorghum_spikein")

###########################################################
##### DESeq2: chilling_am VS chilling_pm ---------------------------------------------
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

##############################################################################
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

### DEGs in nospike method only => in nudelta
DEG_coldampm_onlynospike_nudeltalogFC <- nudelta_logFC[rownames(nudelta_logFC) %in% DEG_coldampm_nospike, ]
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
nudelta_notpass <- intersect(nudelta_lowpadj$geneID, nudelta_lowFC$geneID) #genes = 1654
nudelta_passpval <- setdiff(nudelta_lowFC$geneID, nudelta_lowpadj$geneID) #genes = 84
nudelta_passlogFC <- setdiff(nudelta_lowpadj$geneID, nudelta_lowFC$geneID) #genes = 274

### DEGs in nudelta method only => in nospike
DEG_coldampm_onlynudelta_nospikelogFC <- nospike_logFC[rownames(nospike_logFC) %in% DEG_coldampm_nudelta, ]
# get no. of genes with p-value > 0.05, |logFC| <= 0.5
onlynudelta_nospikelogFC <- as.data.frame(DEG_coldampm_onlynudelta_nospikelogFC)
onlynudelta_nospikelogFC %>%
  rownames_to_column(var = "geneID") %>%
  filter(padj >= 0.05) -> nospike_lowpadj #dim 5919 x 6
onlynudelta_nospikelogFC %>%
  rownames_to_column(var = "geneID") %>%
  filter(log2FoldChange <= 0.5) -> nospike_lowFC #dim 5265 x 6
# make a venn diagram
gplots::venn(list(low_padj = nospike_lowpadj$geneID,
                  low_logFC = nospike_lowFC$geneID))

nospike_notpass <- intersect(nospike_lowpadj$geneID, nospike_lowFC$geneID) #genes = 5042
nospike_passpval <- setdiff(nospike_lowFC$geneID, nospike_lowpadj$geneID) #genes = 223
nospike_passlogFC <- setdiff(nospike_lowpadj$geneID, nospike_lowFC$geneID) #genes = 877

##### check gene annotation ----
Sb_annotation <- read.delim("input/Sbicolor_454_v3.1.1.annotation_info.txt", header = TRUE)
Sb_annotation %>%
  dplyr::select(locusName, Best.hit.arabi.name, arabi.symbol, arabi.defline) %>%
  distinct() -> Sb_annotation1
DEG_coldampm_both_anno <- Sb_annotation1[Sb_annotation1$locusName %in% DEG_coldampm_both,]
DEG_coldampm_nospike_anno <- Sb_annotation1[Sb_annotation1$locusName %in% DEG_coldampm_nospike,]
DEG_coldampm_nudelta_anno <- Sb_annotation1[Sb_annotation1$locusName %in% DEG_coldampm_nudelta,]

nudelta_notpass_anno <- Sb_annotation1[Sb_annotation1$locusName %in% nudelta_notpass,]
nudelta_passpval_anno <- Sb_annotation1[Sb_annotation1$locusName %in% nudelta_passpval,]
nudelta_passlogFC_anno <- Sb_annotation1[Sb_annotation1$locusName %in% nudelta_passlogFC,]

nospike_notpass_anno <- Sb_annotation1[Sb_annotation1$locusName %in% nospike_notpass,]
nospike_passpval_anno <- Sb_annotation1[Sb_annotation1$locusName %in% nospike_passpval,]
nospike_passlogFC_anno <- Sb_annotation1[Sb_annotation1$locusName %in% nospike_passlogFC,]

# circadian gene list
circ_genes <- data.frame("ID" = c("Sobic.007G047400", "Sobic.004G216700", "Sobic.004G281800", 
                                  "Sobic.010G275700", "Sobic.001G352400", "Sobic.002G275100", 
                                  "Sobic.005G044400", "Sobic.002G193000", "Sobic.005G145300", 
                                  "Sobic.003G040900", "Sobic.002G247200", "Sobic.007G210801"),
                         "name" = c("SbLHY", "SbTOC1", "SbRVE6", "SbRVE2", "SbLNK1", "SbPRR5a", 
                                    "SbPRR5b", "SbEFL4", "SbFKF1", "SbGI", "JMJ30a", "JMJ30b"))
circ_genes_list <- circ_genes$ID
circ_DEG_coldampm_both <- circ_genes[circ_genes$ID %in% DEG_coldampm_both,] # 11 genes


########################################################
### DEGs in both method ----
coldampm_both_xy <- xy[xy$geneID %in% DEG_coldampm_both, ]
coldampm_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> coldampm_both_xy2
my_sample_col <- data.frame(method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)),
                            condition = rep(c("chilling_AM", "chilling_PM", 
                                              "chilling_AM", "chilling_PM"), c(4,4,4,4)))
rownames(my_sample_col) <- colnames(coldampm_both_xy2)
my_sample_col %>%
  rownames_to_column(var = "variable") %>%
  mutate(samples = rep(c("replicate_1", "replicate_2", "replicate_3", "replicate_4"), 4)) -> my_sample_col2
coldampm_both_xy_melt <- reshape2::melt(coldampm_both_xy)
coldampm_both_xy_melt2 <- left_join(coldampm_both_xy_melt, my_sample_col2)

# Sobic.007G047400 LHY
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.007G047400") -> LHY
p <- ggplot(LHY, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  labs(title = "Sobic.007G047400",
       subtitle = expression(paste(italic("LHY"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_coldampm_circadian_LHY.png", plot = p,
       path = "~/kanjana/R_analysis/BTx_TOD_nodsp", scale = 1,
       width = 7, height = 4, units = c("in"), dpi = 300)

# Sobic.001G352400 LNK1
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.001G352400") -> LNK1
p <- ggplot(LNK1, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  labs(title = "Sobic.001G352400",
       subtitle = expression(paste(italic("LNK1"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_coldampm_circadian_LNK1.png", plot = p,
       path = "~/kanjana/R_analysis/BTx_TOD_nodsp", scale = 1,
       width = 7, height = 4, units = c("in"), dpi = 300)

# Sobic.004G281800 RVE6
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.004G281800") -> RVE6
p <- ggplot(RVE6, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  labs(title = "Sobic.004G281800",
       subtitle = expression(paste(italic("RVE6"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_coldampm_circadian_RVE6.png", plot = p,
       path = "~/kanjana/R_analysis/BTx_TOD_nodsp", scale = 1,
       width = 7, height = 4, units = c("in"), dpi = 300)

# Sobic.010G275700 RVE2
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.010G275700") -> RVE2
p <- ggplot(RVE2, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  labs(title = "Sobic.010G275700",
       subtitle = expression(paste(italic("RVE2"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_coldampm_circadian_RVE2.png", plot = p,
       path = "~/kanjana/R_analysis/BTx_TOD_nodsp", scale = 1,
       width = 7, height = 4, units = c("in"), dpi = 300)

# Sobic.004G216700 TOC1
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.004G216700") -> TOC1
p <- ggplot(TOC1, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  labs(title = "Sobic.004G216700",
       subtitle = expression(paste(italic("TOC1"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_coldampm_circadian_TOC1.png", plot = p,
       path = "~/kanjana/R_analysis/BTx_TOD_nodsp", scale = 1,
       width = 7, height = 4, units = c("in"), dpi = 300)

# Sobic.002G193000 ELF4
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.002G193000") -> ELF4
p <- ggplot(ELF4, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  labs(title = "Sobic.002G193000",
       subtitle = expression(paste(italic("ELF4"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_coldampm_circadian_ELF4.png", plot = p,
       path = "~/kanjana/R_analysis/BTx_TOD_nodsp", scale = 1,
       width = 7, height = 4, units = c("in"), dpi = 300)

# Sobic.005G145300 FKF
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.005G145300") -> FKF
p <- ggplot(FKF, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  labs(title = "Sobic.005G145300",
       subtitle = expression(paste(italic("FKF"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_coldampm_circadian_FKF.png", plot = p,
       path = "~/kanjana/R_analysis/BTx_TOD_nodsp", scale = 1,
       width = 7, height = 4, units = c("in"), dpi = 300)

# Sobic.009G036400 PIF3
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.009G036400") -> PIF3
p <- ggplot(PIF3, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  labs(title = "Sobic.009G036400",
       subtitle = expression(paste(italic("Phytochrome-interacting transcription factor 3 (PIF3)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_coldampm_PIF3.png", plot = p,
       path = "~/kanjana/R_analysis/BTx_TOD_nodsp", scale = 1,
       width = 7, height = 4, units = c("in"), dpi = 300)

# Sobic.001G480400 COR27
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.001G480400") -> COR27
p <- ggplot(COR27, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  labs(title = "Sobic.001G480400",
       subtitle = expression(paste(italic("Cold regulated gene 27 (COR27))"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_coldampm_COR27.png", plot = p,
       path = "~/kanjana/R_analysis/BTx_TOD_nodsp", scale = 1,
       width = 7, height = 4, units = c("in"), dpi = 300)

########################################################
### DEGs in nospike method only ----
DEG_coldampm_nospike_xy <- xy[xy$geneID %in% DEG_coldampm_nospike, ]
DEG_coldampm_nospike_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> DEG_coldampm_nospike_xy2
my_sample_col <- data.frame(method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)),
                            condition = rep(c("chilling_AM", "chilling_PM", 
                                              "chilling_AM", "chilling_PM"), c(4,4,4,4)))
rownames(my_sample_col) <- colnames(DEG_coldampm_nospike_xy2)
my_sample_col %>%
  rownames_to_column(var = "variable") %>%
  mutate(samples = rep(c("replicate_1", "replicate_2", "replicate_3", "replicate_4"), 4)) -> my_sample_col2
DEG_coldampm_nospike_xy_melt <- reshape2::melt(DEG_coldampm_nospike_xy)
DEG_coldampm_nospike_xy_melt2 <- left_join(DEG_coldampm_nospike_xy_melt, my_sample_col2)

## in nudelta, not pass p-value & logFC cutoff
nudelta_lowpadj %>% dplyr::rename(locusName = geneID) -> nudelta_lowpadj2
z <- left_join(nudelta_notpass_anno, nudelta_lowpadj2)

# Sobic.003G073000 IWS1
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.003G073000") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.003G073000 ", italic("(IWS1)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_IWS1.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.009G257200 CKL2
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.009G257200") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.009G257200 ", italic("(CKL2)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_CKL2.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G537300 REV
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.001G537300") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G537300 ", italic("(REV)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_REV.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.010G085200 MED6
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.010G085200") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G085200 ", italic("(MED6)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_MED6.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.007G186300 AKR2
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.007G186300") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.007G186300 ", italic("(AKR2)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_AKR2.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

## in nudelta, pass p-value cutoff but not logFC cutoff
nudelta_lowFC %>% dplyr::rename(locusName = geneID) -> nudelta_lowFC2
z2 <- left_join(nudelta_passpval_anno, nudelta_lowFC2)

# Sobic.009G160000 bHLH105
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.009G160000") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.009G160000 ", italic("(bHLH105)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_bHLH105.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.004G050000 UBQ10
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.004G050000") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.004G050000 ", italic("(UBQ10)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_UBQ10.png", scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.007G134000 GRF7
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.007G134000") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.007G134000 ", italic("(GRF7)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_GRF7.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

## pass logFC cutoff but not p-value cutoff
z3 <- left_join(nudelta_passlogFC_anno, nudelta_lowpadj2)

# Sobic.001G111366 phyB
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.001G111366") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE)+
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G030200 ", italic("(phyB)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_phyB.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.003G018700 TCP4
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.003G018700") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.003G018700 ", italic("(TCP4)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_TCP4.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.007G133000 FAB1B
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.007G133000") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.007G133000 ", italic("(FAB1B)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_FAB1B.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G273800 ACYB-2
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.001G273800") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G273800 ", italic("(ACYB-2)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_ACYB2.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.007G044100 (NRPB1)
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.007G044100") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.007G044100 ", italic("(NRPB1)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_NRPB1.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.009G101400 (DREB2A) 
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.009G101400") -> DREB2A
p <- ggplot(DREB2A, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.009G101400 ", italic("(DREB2A)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_DREB2A.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

########################################################
### DEGs in nudelta method only ----
DEG_coldampm_nudelta_xy <- xy[xy$geneID %in% DEG_coldampm_nudelta, ]
DEG_coldampm_nudelta_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> DEG_coldampm_nudelta_xy2
my_sample_col <- data.frame(method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)),
                            condition = rep(c("chilling_AM", "chilling_PM", 
                                              "chilling_AM", "chilling_PM"), c(4,4,4,4)))
rownames(my_sample_col) <- colnames(DEG_coldampm_nudelta_xy2)
my_sample_col %>%
  rownames_to_column(var = "variable") %>%
  mutate(samples = rep(c("replicate_1", "replicate_2", "replicate_3", "replicate_4"), 4)) -> my_sample_col2
DEG_coldampm_nudelta_xy_melt <- reshape2::melt(DEG_coldampm_nudelta_xy)
DEG_coldampm_nudelta_xy_melt2 <- left_join(DEG_coldampm_nudelta_xy_melt, my_sample_col2)

## in nospike, not pass p-value & logFC cutoffs
nospike_lowpadj %>% dplyr::rename(locusName = geneID) -> nospike_lowpadj2
a <- left_join(nospike_notpass_anno, nospike_lowpadj2)

# Sobic.003G093800 EIF3A
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.003G093800") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE)+
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.003G093800 ", italic("(EIF3A)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_EIF3A.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.003G257600 TAF11
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.003G257600") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.003G257600 ", italic("(TAF11)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_TAF11.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.003G048600 GRX4
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.003G048600") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.003G048600 ", italic("(GRX4)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_GRX4.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.006G267300 MED9
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.006G267300") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.006G267300 ", italic("(MED9)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_MED9.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.002G010300 CKA1
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.002G010300") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.002G010300 ", italic("(CKA1)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_CKA1.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

## in nospike, pass p-value cutoff but not logFC cutoff
nospike_lowFC %>% dplyr::rename(locusName = geneID) -> nospike_lowFC2
a2 <- left_join(nospike_passpval_anno, nospike_lowFC2)

# Sobic.001G501800 REF6
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.001G501800") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G501800 ", italic("(REF6)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_REF6.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.002G145000 NRPB3
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.002G145000") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.002G145000 ", italic("(NRPB3)"))),  
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_NRPB3.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.004G012800 NRPB4
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.004G012800") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.004G012800 ", italic("(NRPB4)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_NRPB4.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

## in nospike, pass logFC cutoff but not p-value cutoff
a3 <- left_join(nospike_passlogFC_anno, nospike_lowpadj2)

# Sobic.001G023900 GATA16
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.001G023900") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G115800 ", italic("(GATA16)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_GATA16.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G314900 BLH4
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.001G314900") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G314900 ", italic("(BLH4)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_BLH4.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.004G037800 ARF7
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.004G037800") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.004G037800 ", italic("(ARF7)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_ARF7.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.002G092400 (AM-upregulated gene that unique to NuDelta)
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.002G092400") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = "Sobic.002G092400",
       subtitle = expression(paste(italic("60S acidic ribosomal protein family"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_Sobic.002G092400.png", plot = p, scale = 1,
       width = 5, height = 3, units = c("in"), dpi = 300)

# Sobic.002G269400 CBF3
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.002G269400") -> CBF3
p <- ggplot(CBF3, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.002G269400 ", italic("(DREB1A/CBF3)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_coldampm_CBF3.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

