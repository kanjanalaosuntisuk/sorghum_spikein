##### scripts for making spike-in read count table -----

# clear workspace
rm(list = ls())
# load packages
library(tidyverse)
library(RColorBrewer)
library(edgeR)
library(EDASeq)
library(ggsignif)
library(ggpubr)
#setwd("~/sorghum_spikein")

#################################################################
### load spike-in count table
BTxspike <- readRDS(file = "RData/BTxspike_new.RData")

###############################################
# filter out non-detected spikes;  keep reads > 0 in at least 1 samples (removes spike-ins that show 0 in all samples)
spike_filter <- apply(BTxspike, 1, function(x) length(x[x>0])>=1) 
spike_filtered <-  BTxspike[spike_filter,] 
spikes <- rownames(spike_filtered)

##### RLE and PCA plot ---------------------------
x <- as.factor(rep(c("controlam", "controlpm", "chillingam", "chillingpm"), each = 4)) 
x <- factor(x, levels = c("controlam", "controlpm", "chillingam", "chillingpm"))
set <- newSeqExpressionSet(as.matrix(spike_filtered), 
                           phenoData = data.frame(x, row.names=colnames(spike_filtered))) 
colors <- brewer.pal(4, "Set2") 
samplenames <- colnames(spike_filtered)

# make side-by-side RLE and PCA plot
#png(filename="BTxspike_RLE_PCAplot.png", width=7.5, height=3, units="in", res=300)
par(mfrow = c(1, 2))
par(mar = c(5, 4, 1, 1))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], 
        ylab = "Relative log expression", xaxt = "n")
axis(1, at=1:16, labels = samplenames, las = 2, cex.axis = 0.7)
plotPCA(set, col=colors[x], cex=0.75)
#dev.off()

###############################################
##### ploting spikes read counts and concentration ------------------------
### detected vs undetected spike-in
# spike-in read counts
BTxspike
logBTxspike <- log2(BTxspike + 1)
logBTxspike %>%
  rownames_to_column(var = "spike") -> logBTxspike2
logBTxspike2_rshp <- reshape2::melt(logBTxspike2)
treat_controlam <- c("control_am1", "control_am2", "control_am3", "control_am4")
treat_controlpm <- c("control_pm1", "control_pm2", "control_pm3", "control_pm4")
treat_chillingam <- c("chilling_am1", "chilling_am2", "chilling_am3", "chilling_am4")
treat_chillingpm <- c("chilling_pm1", "chilling_pm2", "chilling_pm3", "chilling_pm4")
logBTxspike2_rshp %>%
  mutate(treatment = factor(case_when(variable %in% treat_controlam ~ "control_am",
                                      variable %in% treat_controlpm ~ "control_pm",
                                      variable %in% treat_chillingam ~ "chilling_am",
                                      variable %in% treat_chillingpm ~ "chilling_pm",
                                      TRUE ~ NA_character_))) -> logBTxspike2_rshp2
spike_filtered_list <- rownames(spike_filtered)
logBTxspike2_rshp2$found <- ifelse(logBTxspike2_rshp2$spike %in% spike_filtered_list, "detected", "undetected")

# spike-in concentration
spike_conc <- read.delim("input/spike_conc_BTx.txt", header = TRUE)
str(spike_conc)
names(spike_conc)[names(spike_conc) == "spikeID"] <- "spike"
spike_conc$spike <- as.character(spike_conc$spike)
spike_conc %>%
  dplyr::select(spike, attomol_ul, group) %>%
  mutate(logattomol_ul = log2(attomol_ul)) -> spike_conc_logconc

# join logBTxspike2_rshp2 and spike_conc_logconc
logBTxspike_spikeconc <- inner_join(logBTxspike2_rshp2, spike_conc_logconc, by= "spike")
# make a plot of ERCC spikes
logBTxspike_spikeconcERCC <- logBTxspike_spikeconc[grep("^ERCC", logBTxspike_spikeconc$spike),]
ref_line <- log2(2) 
ggplot(logBTxspike_spikeconcERCC, aes(x = logattomol_ul, y = value, color = found, shape = treatment)) +
  geom_point(stat = "identity") +
  geom_vline(xintercept = ref_line, color = "red") +
  labs(color = "Detection", shape = "") +
  xlab("log2(ERCC concentration (attomoles/μL))") +
  ylab("log2(Average observed read counts)") +
  ggtitle("Detected vs undetected spike-in reads") +
  theme_minimal()

### average read counts from each treatment -------------------------
BTxspike_cp <- BTxspike
BTxspike_cp$control_am = rowMeans(BTxspike_cp[,c(1,2,3,4)])
BTxspike_cp$control_pm = rowMeans(BTxspike_cp[,c(5,6,7,8)])
BTxspike_cp$chilling_am = rowMeans(BTxspike_cp[,c(9,10,11,12)])
BTxspike_cp$chilling_pm = rowMeans(BTxspike_cp[,c(13,14,15,16)])
logBTxspike_Avg <- log2(BTxspike_cp + 1)
logBTxspike_Avg %>%
  rownames_to_column(var = "spike") %>%
  dplyr::select(spike, control_am, control_pm, chilling_am, chilling_pm) -> logBTxspike_Avg2
logBTxspike_Avg2rshp <- reshape2::melt(logBTxspike_Avg2)
logBTxspike_Avg2rshp$found <- ifelse(logBTxspike_Avg2rshp$spike %in% spike_filtered_list, "detected", "undetected")

# join logBTxspike2_rshp2 and spike_conc_logconc
logBTxspikeAvg_spikeconc <- inner_join(logBTxspike_Avg2rshp, spike_conc_logconc, by= "spike")

# make a plot of ERCC spikes
logBTxspikeAvg_spikeconcERCC <- logBTxspikeAvg_spikeconc[grep("^ERCC", logBTxspikeAvg_spikeconc$spike),]
ref_line <- log2(2) 
detectedERCC <-
  ggplot(logBTxspikeAvg_spikeconcERCC, aes(x = logattomol_ul, y = value, color = found, shape = variable)) +
  geom_point(stat = "identity") +
  geom_vline(xintercept = ref_line, color = "red") +
  xlab("log2(ERCC concentration (attomoles/μL))") +
  ylab("log2(Average observed read counts)") +
  ggtitle("Detected vs undetected ERCC reads") +
  labs(shape="Condition", color = "") +
  theme_minimal() +
  theme(plot.title = element_text(size = 14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
#ggsave("detectedERCC.png", plot = detectedERCC, 
#       scale = 1, width = 7, height = 4, units = c("in"), dpi = 300)


###############################################
##### total read counts of non-filtered spike-ins ----
# calculate total spike reads
readsumspike <- colSums(BTxspike)
readsumspikedf <- as.data.frame(readsumspike)
readsumspikedf %>%
  rownames_to_column(var = "sample") %>%
  mutate(condition = rep(c("control_am", "control_pm", "chilling_am", "chilling_pm"), each=4)) %>%
  mutate(set = rep(c("Total"), 16)) %>%
  dplyr::rename(readsum = readsumspike) -> readsumspikedf1

readsumspikedf1$condition <- as.factor(readsumspikedf1$condition)
readsumspikedf1$condition <- factor(readsumspikedf1$condition, 
                                levels = c("control_am", "control_pm", "chilling_am", "chilling_pm"))
readsumspikedf1$sample <- factor(readsumspikedf1$sample, 
                             levels = readsumspikedf1$sample[order(readsumspikedf1$condition)])


###############################################
##### spike and gene read counts -----
##### gene read counts (after removing low expressed genes)
genefiltered <- readRDS(file = "RData/genefiltered.RData")
# combine spike and gene read counts
BTxreadcount <- rbind(genefiltered, BTxspike)
# calculate total read counts
BTxreadsum <- colSums(BTxreadcount)
BTxreadsumdf <- as.data.frame(BTxreadsum)
BTxreadsumdf %>%
  rownames_to_column(var = "sample") %>%
  mutate(condition = rep(c("control_am", "control_pm", "chilling_am", "chilling_pm"), each=4)) -> BTxreadsumdf1
# combine total read and spike read dataframe
BTxreadsum_spike <- inner_join(BTxreadsumdf1, readsumspikedf1)
BTxreadsum_spike %>%
  mutate(spikepertotal = (readsum/BTxreadsum)*100) -> BTxreadsum_spike1
BTxreadsum_spike1$condition <- factor(BTxreadsum_spike1$condition, 
                                      levels = c("control_am", "control_pm", "chilling_am", "chilling_pm"))
BTxreadsum_spike1$sample <- factor(BTxreadsum_spike1$sample, levels = c("control_am1", "control_am2", "control_am3", "control_am4",
                                                                        "control_pm1", "control_pm2", "control_pm3", "control_pm4",
                                                                        "chilling_am1", "chilling_am2", "chilling_am3", "chilling_am4",
                                                                        "chilling_pm1", "chilling_pm2", "chilling_pm3", "chilling_pm4"))

##### one-way ANOVA ----
## Compute the analysis of variance
res.aov <- aov(spikepertotal ~ condition, data = BTxreadsum_spike1)
# Summary of the analysis
summary(res.aov)
# p-value = 0.336 > 0.05 -> no difference between groups

## check homogeneity of variances
library(car)
leveneTest(spikepertotal ~ condition, data = BTxreadsum_spike1) 
# p-value = 0.5305 > 0.05 -> data follow homogeneity

## check normality
# Run Shapiro-Wilk test
aov_residuals <- residuals(object = res.aov)
shapiro.test(x = aov_residuals)
# p-value = 0.1717 > 0.05 -> data follow normality

### make a boxplot of spike-in proportion
colors <- brewer.pal(4, "Set2") 
BTxspike_spikeratio <-
  ggplot(BTxreadsum_spike1, aes(x = condition, y = spikepertotal)) +
  geom_boxplot(fill = colors) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  geom_signif(y_position = 0.08, xmin = 1, xmax = 4, 
              annotation = c("p-value = 0.336"), tip_length = 0.1, textsize = 4.5) +
  ylim(0, 0.08) +
  xlab("") +
  ylab("Spike-in proportion (%)") +
  theme_minimal() +
  theme(axis.text=element_text(size=12),
        axis.title = element_text(size = 14))
# save plot 
#ggsave("BTxspike_spikeratio_boxplot_new.png", plot = BTxspike_spikeratio,
#       scale = 1, width = 7, height = 4, units = c("in"), dpi = 300)


####################################################################
### correlation between gene read counts and spike-in read counts
genefiltered_df <- as.data.frame(genefiltered)
genefiltered_sum <- colSums(genefiltered_df)
genefiltered_sumdf <- as.data.frame(genefiltered_sum)
genefiltered_sumdf %>%
  rownames_to_column(var = "sample") %>%
  mutate(condition = rep(c("control_am", "control_pm", "chilling_am", "chilling_pm"), each=4)) -> genefiltered_sumdf1
gene_spike_sum <- full_join(genefiltered_sumdf1, readsumspikedf1)

## Preleminary test to check the test assumptions
# Shapiro-Wilk normality test 
shapiro.test(gene_spike_sum$genefiltered_sum) #p-value = 0.3111
shapiro.test(gene_spike_sum$readsum) #p-value = 0.1083
# look at the normality plot 
ggqqplot(gene_spike_sum$genefiltered_sum, ylab = "Total filtered gene reads")
ggqqplot(gene_spike_sum$readsum, ylab = "Total spike-in reads")

## Pearson correlation test
res <- cor.test(gene_spike_sum$readsum, gene_spike_sum$genefiltered_sum, 
                method = "pearson")
res 
#cor = 0.9142788 
#p-value = 7.253e-07 -> significantly correlated 

# save correlation plot
#png(filename="BTxgenespike_pearson_cor.png", width=7, height=4, units="in", res=300)
ggscatter(gene_spike_sum, x = "readsum", y = "genefiltered_sum", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 5,
          xlab = "Total spike-in reads", ylab = "Total filtered gene reads")
#dev.off()

