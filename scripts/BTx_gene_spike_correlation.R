##### scripts for correlation between gene reads and spike-in reads -----

# clear workspace
rm(list = ls())
# load packages
library(tidyverse)
library(RColorBrewer)
library(edgeR)
library(EDASeq)
library(ggpubr)
setwd("~/sorghum_spikein")

####################################################
##### spike counts ---------------------------------
# load spike-in read count table
BTxspike <- readRDS("RData/BTxspike_new.RData")
# filter out non-detected spikes;  keep reads > 0 in at least 1 samples
spike_filter <- apply(BTxspike, 1, function(x) length(x[x>0])>=1) 
spike_filtered <-  BTxspike[spike_filter,] 
spikes <- rownames(spike_filtered)

# calculate total spike reads
readsumspike <- colSums(BTxspike)
readsumspikedf <- as.data.frame(readsumspike)
readsumspikedf %>%
  rownames_to_column(var = "sample") %>%
  mutate(condition = rep(c("control_am", "control_pm", "chilling_am", "chilling_pm"), each=4)) %>%
  mutate(set = rep(c("Total"), 16)) %>%
  dplyr::rename(readsum = readsumspike) -> readsumspikedf1
readsumspikedf1 %>%
  dplyr::select(sample, readsum) %>%
  spread(key = sample, value = readsum) -> readsumspikedf2
rownames(readsumspikedf2) <- c("total_spike")

##### gene read count data ------------------------
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

### combine spike-in and gene counts
genefiltered_df <- as.data.frame(genefiltered)
BTx_spike_df <- bind_rows(readsumspikedf2, genefiltered_df)

## Pearson correlation test
test_df <- BTx_spike_df
test_dft <- t(test_df)
test_dft_mt <- as.data.frame(test_dft)
cor_res <- cor(test_dft_mt[-1], test_dft_mt$total_spike)
cor_res_df <- as.data.frame(cor_res)
names(cor_res_df)[names(cor_res_df) == "V1"] <- "pearson_corr"
cor_res_df %>%
  arrange(desc(pearson_corr)) %>%
  rownames_to_column(var = "GeneID") -> cor_res_df2
# export table
write.table(cor_res_df2, file = "output/BTx_gene_spike_correlation.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

##### check gene annotation ----
Sb_annotation <- read.delim("input/Sbicolor_454_v3.1.1.annotation_info.txt", header = TRUE)
Sb_annotation %>%
  dplyr::select(locusName, Best.hit.arabi.name, arabi.symbol, arabi.defline) %>%
  distinct() -> Sb_annotation1
names(Sb_annotation1)[names(Sb_annotation1) == "locusName"] <- "GeneID"
cor_res_df2_anno <- left_join(cor_res_df2, Sb_annotation1)
# export table
write.table(cor_res_df2_anno, file = "output/BTx_gene_spike_correlation_annotation.txt", sep = "\t", 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

# select genes with highest corr
BTx_spike_df %>%
  rownames_to_column(var = "GeneID") %>%
  filter(GeneID == cor_res_df2$GeneID[1] | GeneID == "total_spike") %>%
  column_to_rownames(var = "GeneID") -> corr_gene1
corr_gene2 <- t(corr_gene1)
corr_gene3 <- as.data.frame(corr_gene2)
# save correlation plot
#png(filename="BTxgenespike_pearson_cor_Sobic001G006700.png", width=5, height=3, units="in", res=300)
ggscatter(corr_gene3, x = "total_spike", y = "Sobic.001G006700", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 5,
          xlab = "Total spike-in reads", ylab = "Filtered gene reads",
          title = "Sobic.001G006700")
#dev.off()

# select genes with lowest corr
BTx_spike_df %>%
  rownames_to_column(var = "GeneID") %>%
  filter(GeneID == cor_res_df2$GeneID[17651] | GeneID == "total_spike") %>%
  column_to_rownames(var = "GeneID") -> corr_geneLast
corr_geneLast2 <- t(corr_geneLast)
corr_geneLast3 <- as.data.frame(corr_geneLast2)
#png(filename="BTxgenespike_pearson_cor_Sobic007G153800.png", width=5, height=3, units="in", res=300)
ggscatter(corr_geneLast3, x = "total_spike", y = "Sobic.007G153800", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 5,
          xlab = "Total spike-in reads", ylab = "Filtered gene reads",
          title = "Sobic.007G153800")
#dev.off()
