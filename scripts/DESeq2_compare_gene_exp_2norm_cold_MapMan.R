##### scripts for making plots of DEGs from MapMan analysis ----
##### Effect of chilling stress on gene expression in AM and PM ----

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

##### candidate genes from MapMan ----
# Sobic.003G045400 NRPB7
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.003G045400") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.003G045400 ", italic("(NRPB7)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_coldampm_NRPB7.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.006G127500 RPB5A
coldampm_both_xy_melt2 %>%
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
ggsave("MapMan_geneExp_coldampm_RPB5A.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.005G030600 RPB5B
coldampm_both_xy_melt2 %>%
  filter(geneID == "Sobic.005G030600") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.005G030600 ", italic("(NRPE5)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_coldampm_RPB5B.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)


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

# Sobic.007G044100 NRPB1
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.007G044100") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
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
ggsave("MapMan_geneExp_coldampm_NRPB1.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G155000 NRPB2
DEG_coldampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.001G155000") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G155000 ", italic("(NRPB2)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_coldampm_NRPB2.png", plot = last_plot(), scale = 1,
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

# Sobic.004G012800 NRPB4
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.004G012800") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
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
ggsave("MapMan_geneExp_coldampm_NRPB4.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G013500 NRPB8B
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.001G013500") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G013500 ", italic("(NRPB8B)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_coldampm_NRPB8B.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.007G068800 NRPB12
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.007G068800") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.007G068800 ", italic("(NRPB12)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_coldampm_NRPB12.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.008G058300
DEG_coldampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.008G058300") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.008G058300 ", italic("(NRPB10)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_coldampm_RPB10.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)
