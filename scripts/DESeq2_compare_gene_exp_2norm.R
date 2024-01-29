##### scripts for making plots of selected DEGs from DESeq2 analysis ----
##### Effect of time of day on gene expression under control and chilling stress ----

library(DESeq2)
library(tidyverse)
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

##############################################################################
### load AM upregulated gene list----
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

### DEGs in nospike method only => in nudelta
DEG_ctrlampm_onlynospike_nudeltalogFC <- nudelta_logFC[rownames(nudelta_logFC) %in% DEG_ctrlampm_nospike, ]
# get no. of genes with p-value > 0.05, |logFC| <= 0.5
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
nudelta_notpass <- intersect(nudelta_lowpadj$geneID, nudelta_lowFC$geneID) #genes = 904
nudelta_passpval <- setdiff(nudelta_lowFC$geneID, nudelta_lowpadj$geneID) #genes = 148
nudelta_passlogFC <- setdiff(nudelta_lowpadj$geneID, nudelta_lowFC$geneID) #genes = 345

### DEGs in nudelta method only => in nospike
DEG_ctrlampm_onlynudelta_nospikelogFC <- nospike_logFC[rownames(nospike_logFC) %in% DEG_ctrlampm_nudelta, ]
# get no. of genes with p-value > 0.05, |logFC| <= 0.5
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
nospike_notpass <- intersect(nospike_lowpadj$geneID, nospike_lowFC$geneID) #genes = 1475
nospike_passpval <- setdiff(nospike_lowFC$geneID, nospike_lowpadj$geneID) #genes = 189
nospike_passlogFC <- setdiff(nospike_lowpadj$geneID, nospike_lowFC$geneID) #genes = 716

##### check gene annotation ----
Sb_annotation <- read.delim("input/Sbicolor_454_v3.1.1.annotation_info.txt", header = TRUE)
Sb_annotation %>%
  dplyr::select(locusName, Best.hit.arabi.name, arabi.symbol, arabi.defline) %>%
  distinct() -> Sb_annotation1
DEG_ctrlampm_both_anno <- Sb_annotation1[Sb_annotation1$locusName %in% DEG_ctrlampm_both,]
DEG_ctrlampm_nospike_anno <- Sb_annotation1[Sb_annotation1$locusName %in% DEG_ctrlampm_nospike,]
DEG_ctrlampm_nudelta_anno <- Sb_annotation1[Sb_annotation1$locusName %in% DEG_ctrlampm_nudelta,]

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
circ_DEG_ctrlampm_both <- circ_genes[circ_genes$ID %in% DEG_ctrlampm_both,] # 11 genes

########################################################
### DEGs in both method ----
ctrlampm_both_xy <- xy[xy$geneID %in% DEG_ctrlampm_both, ]
ctrlampm_both_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> ctrlampm_both_xy2
my_sample_col <- data.frame(method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)),
                            condition = rep(c("control_AM", "control_PM", 
                                              "control_AM", "control_PM"), c(4,4,4,4)))
rownames(my_sample_col) <- colnames(ctrlampm_both_xy2)
my_sample_col %>%
  rownames_to_column(var = "variable") %>%
  mutate(samples = rep(c("replicate_1", "replicate_2", "replicate_3", "replicate_4"), 4)) -> my_sample_col2
ctrlampm_both_xy_melt <- reshape2::melt(ctrlampm_both_xy)
ctrlampm_both_xy_melt2 <- left_join(ctrlampm_both_xy_melt, my_sample_col2)

# Sobic.007G047400 LHY
ctrlampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.007G047400") -> LHY
p <- ggplot(LHY, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = after_stat(y), ymin = after_stat(y)), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.007G047400 ", italic("(LHY)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_ctrlampm_circadian_LHY.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G352400 LNK1
ctrlampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.001G352400") -> LNK1
p <- ggplot(LNK1, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G352400 ", italic("(LNK1)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_ctrlampm_circadian_LNK1.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.004G281800 RVE6
ctrlampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.004G281800") -> RVE6
p <- ggplot(RVE6, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.004G281800 ", italic("(RVE6)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_ctrlampm_circadian_RVE6.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.010G275700 RVE2
ctrlampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.010G275700") -> RVE2
p <- ggplot(RVE2, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G275700 ", italic("(RVE2)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_ctrlampm_circadian_RVE2.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.004G216700 TOC1
ctrlampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.004G216700") -> TOC1
p <- ggplot(TOC1, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.004G216700 ", italic("(TOC1)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_ctrlampm_circadian_TOC1.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.002G193000 ELF4
ctrlampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.002G193000") -> ELF4
p <- ggplot(ELF4, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.002G193000 ", italic("(ELF4)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_ctrlampm_circadian_ELF4.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.005G145300 FKF
ctrlampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.005G145300") -> FKF
p <- ggplot(FKF, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.005G145300 ", italic("(FKF1)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
# save plot 
ggsave("geneExp_ctrlampm_circadian_FKF.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

##### candidate genes ----
# Sobic.006G127500 NRPB5A
ctrlampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.006G127500") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.006G127500 ", italic("(NRPB5A)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_ctrlampm_NRPB5A.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.008G128900 NRPB9B
ctrlampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.008G128900") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.008G128900 ", italic("(NRPB9B)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_ctrlampm_NRPB9B.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)


########################################################
### DEGs in nospike method only ----
DEG_ctrlampm_nospike_xy <- xy[xy$geneID %in% DEG_ctrlampm_nospike, ]
DEG_ctrlampm_nospike_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> DEG_ctrlampm_nospike_xy2
my_sample_col <- data.frame(method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)),
                            condition = rep(c("control_AM", "control_PM", 
                                              "control_AM", "control_PM"), c(4,4,4,4)))
rownames(my_sample_col) <- colnames(DEG_ctrlampm_nospike_xy2)
my_sample_col %>%
  rownames_to_column(var = "variable") %>%
  mutate(samples = rep(c("replicate_1", "replicate_2", "replicate_3", "replicate_4"), 4)) -> my_sample_col2
DEG_ctrlampm_nospike_xy_melt <- reshape2::melt(DEG_ctrlampm_nospike_xy)
DEG_ctrlampm_nospike_xy_melt2 <- left_join(DEG_ctrlampm_nospike_xy_melt, my_sample_col2)

## in nudelta, not pass p-value & logFC cutoff
nudelta_lowpadj %>% dplyr::rename(locusName = geneID) -> nudelta_lowpadj2
z <- left_join(nudelta_notpass_anno, nudelta_lowpadj2)

# Sobic.001G146000 tubulin beta 8
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.001G146000") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G146000 ", italic("(TUB8)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_tub_b8.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G411400 PRR7
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.001G411400") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G411400 ", italic("(PRR7)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_PRR7.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.007G155900 ABF3
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.007G155900") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.007G155900 ", italic("(ABF3)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_ABF3.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.003G305000 TCP5
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.003G305000") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.003G305000 ", italic("(TCP5)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_TCP5.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.009G030100 bZIP53
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.009G030100") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.009G030100 ", italic("(bZIP53)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_bZIP53.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

## in nudelta, pass p-value cutoff but not logFC cutoff
nudelta_lowFC %>% dplyr::rename(locusName = geneID) -> nudelta_lowFC2
z2 <- left_join(nudelta_passpval_anno, nudelta_lowFC2)

# Sobic.009G005900 ACT7
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.009G005900") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.009G005900 ", italic("(ACT7)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_ACT7.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.010G236300 auxin response factor 16
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.010G236300") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G236300 ", italic("(ARF16)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_ARF16.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.008G015200 GRL
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.008G015200") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.008G015200 ", italic("(GRL)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_GRL.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G039700 bZIP19
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.010G267500") -> bZIP19
p <- ggplot(bZIP19, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G039700 ", italic("(bZIP19)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_bZIP19.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

## pass logFC cutoff but not p-value cutoff
z3 <- left_join(nudelta_passlogFC_anno, nudelta_lowpadj2)

# Sobic.010G030200 HAT14
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.010G030200") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G030200 ", italic("(HAT14)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_HAT14.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.003G018700 TCP4
DEG_ctrlampm_nospike_xy_melt2 %>%
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
ggsave("geneExp_ctrlampm_TCP4.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.010G199000 bZIP43
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.010G199000") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G199000 ", italic("(bZIP43)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_bZIP43.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.010G077200 GRF5
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.010G077200") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G077200 ", italic("(GRF5)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_GRF5.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

########################################################
### DEGs in nudelta method only ----
DEG_ctrlampm_nudelta_xy <- xy[xy$geneID %in% DEG_ctrlampm_nudelta, ]
DEG_ctrlampm_nudelta_xy %>%
  remove_rownames() %>%
  column_to_rownames(var = "geneID") -> DEG_ctrlampm_nudelta_xy2
my_sample_col <- data.frame(method = rep(c("Median_of_Ratio", "NuDelta"), c(8,8)),
                            condition = rep(c("control_AM", "control_PM", 
                                              "control_AM", "control_PM"), c(4,4,4,4)))
rownames(my_sample_col) <- colnames(DEG_ctrlampm_nudelta_xy2)
my_sample_col %>%
  rownames_to_column(var = "variable") %>%
  mutate(samples = rep(c("replicate_1", "replicate_2", "replicate_3", "replicate_4"), 4)) -> my_sample_col2
DEG_ctrlampm_nudelta_xy_melt <- reshape2::melt(DEG_ctrlampm_nudelta_xy)
DEG_ctrlampm_nudelta_xy_melt2 <- left_join(DEG_ctrlampm_nudelta_xy_melt, my_sample_col2)

## in nospike, not pass p-value & logFC cutoffs
nospike_lowpadj %>% dplyr::rename(locusName = geneID) -> nospike_lowpadj2
a <- left_join(nospike_notpass_anno, nospike_lowpadj2)

# Sobic.001G162100 WRKY21
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.001G162100") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G162100 ", italic("(WRKY21)"))),  
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_WRKY21.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.010G209400 TAF5
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.010G209400") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G209400 ", italic("(TAF5)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_TAF5.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.003G003800 ARF11
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.003G003800") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.003G003800 ", italic("(ARF11)"))),  
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_ARF11.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.007G103100 MED14
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.007G103100") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.007G103100 ", italic("(MED14)"))),  
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_MED14.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G387800 TOM40
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.001G387800") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G387800 ", italic("(TOM40)"))),  
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_TOM40.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.004G019200 NRPA2
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.004G019200") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.004G019200 ", italic("(NRPA2)"))),  
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_NRPA2.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

## in nospike, pass p-value cutoff but not logFC cutoff
# Sobic.002G354900 SCL5
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.002G354900") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.002G354900 ", italic("(SCL5)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_SCL5.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G014500 NRPB8B
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.001G014500") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G014500 ", italic("(NRPB8B)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_NRPB8B.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.010G251100 EIF4A1
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.010G251100") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G251100 ", italic("(EIF4A1)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_EIF4A1.png", plot = p,scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

## in nospike, pass logFC cutoff but not p-value cutoff
# Sobic.010G115800 COL2
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.010G115800") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.010G115800 ", italic("(COL2)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_COL2.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.004G037800 ARF7
DEG_ctrlampm_nudelta_xy_melt2 %>%
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
ggsave("geneExp_ctrlampm_ARF7.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.006G142400 GTF2H2
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.006G142400") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.006G142400 ", italic("(GTF2H2)"))), 
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_GTF2H2.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.003G293400 STT3A
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.003G293400") -> test
p <- ggplot(test, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.003G293400 ", italic("(STT3A)"))),  
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
# save plot 
ggsave("geneExp_ctrlampm_STT3A.png", plot = p, scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)
