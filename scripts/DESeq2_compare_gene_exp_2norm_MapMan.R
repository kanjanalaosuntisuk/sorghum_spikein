##### scripts for making plots of DEGs from MapMan analysis ----
##### Effect of Time of Day on gene expression under control and chilling stress ----

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

##### candidate genes from MapMan ----
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
        #legend.text = element_text(size=15),
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
########################################################
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

##### candidate genes from MapMan
# Sobic.006G209500
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.006G209500") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.006G209500 ", italic("(TFIIE2)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_ctrlampm_TFIIE2.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.003G442601
DEG_ctrlampm_nospike_xy_melt2 %>%
  filter(geneID == "Sobic.003G442601") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.003G442601 ", italic("(TFIIA2)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_ctrlampm_TFIIA.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)


########################################################
### DEGs in nudelta method only ----
########################################################
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

##### candidate genes from MapMan
# Sobic.002G145000 NRPB3
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.002G145000") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
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
ggsave("MapMan_geneExp_ctrlampm_NRPB3.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G179000 NRPB6A
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.001G179000") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
  geom_point(show.legend = FALSE, 
             position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.15)) +
  stat_summary(fun = mean, geom = "errorbar", position = position_dodge(width = 0.3),
               aes(ymax = ..y.., ymin = ..y..), width = 0.25, show.legend = FALSE) +
  scale_x_discrete(labels=c("Median_of_Ratio" = "Median of Ratio", "NuDelta" = "Spike-in")) +
  labs(title = expression(paste("Sobic.001G179000 ", italic("(NRPB6A)"))),
       y = "log2(normalized count)", x = "") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title = element_text(size=12), axis.text = element_text(size=12))
ggsave("MapMan_geneExp_ctrlampm_NRPB6A.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)

# Sobic.001G014500  NRPB8B
DEG_ctrlampm_nudelta_xy_melt2 %>%
  filter(geneID == "Sobic.001G014500") -> df_forplot
ggplot(df_forplot, aes(x = method, y = value, color = condition)) + 
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
ggsave("MapMan_geneExp_ctrlampm_NRPB8B.png", plot = last_plot(), scale = 1,
       width = 3.5, height = 2, units = c("in"), dpi = 300)
