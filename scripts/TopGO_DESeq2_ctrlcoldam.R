##### scripts for TopGO analysis of DE genes from DESeq2 ----
##### Effect of chilling stress on gene expression in AM

library(tidyverse)
library(topGO)
library(DESeq2)
library(pheatmap)
#library(GOplot)
setwd("~/sorghum_spikein")

####################################################################
##### TopGO analysis ---------------------------
### load SlimGO for topGO ----
Sb_annotation_good <- readMappings(file = "input/SlimGO_Sobic_forTopGO.txt", IDsep = ",")
# load BTx filtered gene list
genefilteredList <- readRDS("RData/genefilteredList_forTopGO.RData")

####################################################################
##### controlam VS chillingam ---------------------------------------
##### control_am upregulated genes
### DESeq2 nospike normalization ---------------------
nospike_ctrl <- read.delim("output/DEGs_DESeq2_nospike_ctrlcoldam_ctrl.txt", header = FALSE, sep = "")
nospike_ctrllist <- nospike_ctrl$V1
nospike_ctrl_geneList <- factor(as.integer(genefilteredList %in% nospike_ctrllist))
names(nospike_ctrl_geneList) <- genefilteredList
### make a topGOdata object for BP
nospike_ctrl_GOdata <- new("topGOdata", description = "DESeq2 nospike AM up", 
                         ontology = "BP", allGenes = nospike_ctrl_geneList,
                         annot = annFUN.gene2GO, gene2GO = Sb_annotation_good)
# run Fisher's exact test
resultFisher_BP <- runTest(nospike_ctrl_GOdata, algorithm = "weight01", statistic = "fisher")
# get a result table
nospike_ctrl_allRes <- GenTable(nospike_ctrl_GOdata, weight01Fisher = resultFisher_BP,
                                     orderBy = "weight01Fisher", topNodes = 100, numChar=1000)
View(nospike_ctrl_allRes)

# get genes in GO terms which weight01Fisher < 0.01 in BP
gene.universe <- genes(nospike_ctrl_GOdata) 
sig.genes <- sigGenes(nospike_ctrl_GOdata)
nospike_ctrl_allRes$genes <- sapply(nospike_ctrl_allRes$GO.ID, function(x)
{
  genes <- genesInTerm(nospike_ctrl_GOdata, x)
  genes[[1]][genes[[1]] %in% sig.genes] # myGenes is the queried gene list
})
nospike_ctrl_allRes2 <- subset(nospike_ctrl_allRes, 
                             as.numeric(nospike_ctrl_allRes[,"weight01Fisher"])<0.01 | grepl("e-", nospike_ctrl_allRes[,"weight01Fisher"]))
nospike_ctrl_allRes2$genes <-vapply(nospike_ctrl_allRes2$genes, paste, collapse = ",", character(1L))
nospike_ctrl_allRes2 %>%
  mutate(adj_pval = as.numeric(weight01Fisher)) %>%
  dplyr::rename(ID = GO.ID, Genes = genes) %>%
  mutate(category = as.factor(rep(c("BP"), nrow(nospike_ctrl_allRes2)))) %>%
  dplyr::select(category, ID, Term, Significant, Genes, adj_pval) -> nospike_ctrl_allRes3
nospike_ctrl_allRes3 %>%
  mutate(logadj_pval = -log(adj_pval)) -> nospike_ctrl_allRes4
# save RData
saveRDS(nospike_ctrl_allRes4, file = "RData/TopGO_ctrlcoldam_nospike.RData")


#####################################################
### DESeq2 nudelta normalization ---------------------
nudelta_ctrl <- read.delim("output/DEGs_DESeq2_nudelta_ctrlcoldam_ctrl.txt", header = FALSE, sep = "")
nudelta_ctrllist <- nudelta_ctrl$V1
nudelta_ctrl_geneList <- factor(as.integer(genefilteredList %in% nudelta_ctrllist))
names(nudelta_ctrl_geneList) <- genefilteredList
### make a topGOdata object for BP
nudelta_ctrl_GOdata <- new("topGOdata", description = "DESeq2 nospike ctrl up", 
                         ontology = "BP", allGenes = nudelta_ctrl_geneList,
                         annot = annFUN.gene2GO, gene2GO = Sb_annotation_good)
# run Fisher's exact test
resultFisher_BP <- runTest(nudelta_ctrl_GOdata, algorithm = "weight01", statistic = "fisher")
# get a result table
nudelta_ctrl_allRes <- GenTable(nudelta_ctrl_GOdata, weight01Fisher = resultFisher_BP,
                              orderBy = "weight01Fisher", topNodes = 50, numChar=1000)
View(nudelta_ctrl_allRes)

# get genes in GO terms which weight01Fisher < 0.01 in BP
gene.universe <- genes(nudelta_ctrl_GOdata) 
sig.genes <- sigGenes(nudelta_ctrl_GOdata)
nudelta_ctrl_allRes$genes <- sapply(nudelta_ctrl_allRes$GO.ID, function(x)
{
  genes <- genesInTerm(nudelta_ctrl_GOdata, x)
  genes[[1]][genes[[1]] %in% sig.genes] # myGenes is the queried gene list
})
nudelta_ctrl_allRes2 <- subset(nudelta_ctrl_allRes, 
                             as.numeric(nudelta_ctrl_allRes[,"weight01Fisher"])<0.01 | grepl("e-", nudelta_ctrl_allRes[,"weight01Fisher"]))
nudelta_ctrl_allRes2$genes <-vapply(nudelta_ctrl_allRes2$genes, paste, collapse = ",", character(1L))
nudelta_ctrl_allRes2 %>%
  mutate(adj_pval = as.numeric(weight01Fisher)) %>%
  dplyr::rename(ID = GO.ID, Genes = genes) %>%
  mutate(category = as.factor(rep(c("BP"), nrow(nudelta_ctrl_allRes2)))) %>%
  dplyr::select(category, ID, Term, Significant, Genes, adj_pval) -> nudelta_ctrl_allRes3
nudelta_ctrl_allRes3 %>%
  mutate(logadj_pval = -log(adj_pval)) -> nudelta_ctrl_allRes4
# save RData
saveRDS(nudelta_ctrl_allRes4, file = "RData/TopGO_ctrlcoldam_nudelta.RData")


#####################################################
### DESeq2 spikesf normalization ---------------------
spikesf_ctrl <- read.delim("output/DEGs_DESeq2_spikesf_ctrlcoldam_ctrl.txt", header = FALSE, sep = "")
spikesf_ctrllist <- spikesf_ctrl$V1
spikesf_ctrl_geneList <- factor(as.integer(genefilteredList %in% spikesf_ctrllist))
names(spikesf_ctrl_geneList) <- genefilteredList
### make a topGOdata object for BP
spikesf_ctrl_GOdata <- new("topGOdata", description = "DESeq2 nospike ctrl up", 
                         ontology = "BP", allGenes = spikesf_ctrl_geneList,
                         annot = annFUN.gene2GO, gene2GO = Sb_annotation_good)
# run Fisher's exact test
resultFisher_BP <- runTest(spikesf_ctrl_GOdata, algorithm = "weight01", statistic = "fisher")
# get a result table
spikesf_ctrl_allRes <- GenTable(spikesf_ctrl_GOdata, weight01Fisher = resultFisher_BP,
                              orderBy = "weight01Fisher", topNodes = 30)
View(spikesf_ctrl_allRes)

# get genes in GO terms which weight01Fisher < 0.01 in BP
gene.universe <- genes(spikesf_ctrl_GOdata) 
sig.genes <- sigGenes(spikesf_ctrl_GOdata)
spikesf_ctrl_allRes$genes <- sapply(spikesf_ctrl_allRes$GO.ID, function(x)
{
  genes <- genesInTerm(spikesf_ctrl_GOdata, x)
  genes[[1]][genes[[1]] %in% sig.genes] # myGenes is the queried gene list
})
spikesf_ctrl_allRes2 <- subset(spikesf_ctrl_allRes, 
                             as.numeric(spikesf_ctrl_allRes[,"weight01Fisher"])<0.01 | grepl("e-", spikesf_ctrl_allRes[,"weight01Fisher"]))
spikesf_ctrl_allRes2$genes <-vapply(spikesf_ctrl_allRes2$genes, paste, collapse = ",", character(1L))
spikesf_ctrl_allRes2 %>%
  mutate(adj_pval = as.numeric(weight01Fisher)) %>%
  dplyr::rename(ID = GO.ID, Genes = genes) %>%
  mutate(category = as.factor(rep(c("BP"), nrow(spikesf_ctrl_allRes2)))) %>%
  dplyr::select(category, ID, Term, Significant, Genes, adj_pval) -> spikesf_ctrl_allRes3
spikesf_ctrl_allRes3 %>%
  mutate(logadj_pval = -log(adj_pval)) -> spikesf_ctrl_allRes4


#####################################################
### DESeq2 RUV normalization ---------------------
RUV_ctrl <- read.delim("output/DEGs_DESeq2_RUV_ctrlcoldam_ctrl.txt", header = FALSE, sep = "")
RUV_ctrllist <- RUV_ctrl$V1
RUV_ctrl_geneList <- factor(as.integer(genefilteredList %in% RUV_ctrllist))
names(RUV_ctrl_geneList) <- genefilteredList
### make a topGOdata object for BP
RUV_ctrl_GOdata <- new("topGOdata", description = "DESeq2 nospike ctrl up", 
                         ontology = "BP", allGenes = RUV_ctrl_geneList,
                         annot = annFUN.gene2GO, gene2GO = Sb_annotation_good)
# run Fisher's exact test
resultFisher_BP <- runTest(RUV_ctrl_GOdata, algorithm = "weight01", statistic = "fisher")
# get a result table
RUV_ctrl_allRes <- GenTable(RUV_ctrl_GOdata, weight01Fisher = resultFisher_BP,
                              orderBy = "weight01Fisher", topNodes = 30)
View(RUV_ctrl_allRes)

# get genes in GO terms which weight01Fisher < 0.05 in BP
gene.universe <- genes(RUV_ctrl_GOdata) 
sig.genes <- sigGenes(RUV_ctrl_GOdata)
RUV_ctrl_allRes$genes <- sapply(RUV_ctrl_allRes$GO.ID, function(x)
{
  genes <- genesInTerm(RUV_ctrl_GOdata, x)
  genes[[1]][genes[[1]] %in% sig.genes] # myGenes is the queried gene list
})
RUV_ctrl_allRes2 <- subset(RUV_ctrl_allRes, 
                             as.numeric(RUV_ctrl_allRes[,"weight01Fisher"])<0.01 | grepl("e-", RUV_ctrl_allRes[,"weight01Fisher"]))
RUV_ctrl_allRes2$genes <-vapply(RUV_ctrl_allRes2$genes, paste, collapse = ",", character(1L))
RUV_ctrl_allRes2 %>%
  mutate(adj_pval = as.numeric(weight01Fisher)) %>%
  dplyr::rename(ID = GO.ID, Genes = genes) %>%
  mutate(category = as.factor(rep(c("BP"), nrow(RUV_ctrl_allRes2)))) %>%
  dplyr::select(category, ID, Term, Significant, Genes, adj_pval) -> RUV_ctrl_allRes3
RUV_ctrl_allRes3 %>%
  mutate(logadj_pval = -log(adj_pval)) -> RUV_ctrl_allRes4


############################################################
##### combine GO terms from four normalization methods ----
## combine logadj_pval
nospike_ctrl_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(MedianofRatio = logadj_pval) -> nospike_ctrl_allRes5
nudelta_ctrl_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(NuDelta = logadj_pval) -> nudelta_ctrl_allRes5
spikesf_ctrl_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(SpikeSizeFactor = logadj_pval) -> spikesf_ctrl_allRes5
RUV_ctrl_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(RUV = logadj_pval) -> RUV_ctrl_allRes5
ctrl_up_GO <- full_join(nospike_ctrl_allRes5, nudelta_ctrl_allRes5)
ctrl_up_GO <- full_join(ctrl_up_GO, spikesf_ctrl_allRes5)
ctrl_up_GO <- full_join(ctrl_up_GO, RUV_ctrl_allRes5)
ctrl_up_GO %>%
  dplyr::select(-ID) %>%
  column_to_rownames(var = "Term") -> ctrl_up_GO

## combine no. of genes
nospike_ctrl_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(MedianofRatio = Significant) -> nospike_ctrl_allRes6
nudelta_ctrl_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(NuDelta = Significant) -> nudelta_ctrl_allRes6
spikesf_ctrl_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(SpikeSizeFactor = Significant) -> spikesf_ctrl_allRes6
RUV_ctrl_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(RUV = Significant) -> RUV_ctrl_allRes6
ctrl_up_GO_genes <- full_join(nospike_ctrl_allRes6, nudelta_ctrl_allRes6)
ctrl_up_GO_genes <- full_join(ctrl_up_GO_genes, spikesf_ctrl_allRes6)
ctrl_up_GO_genes <- full_join(ctrl_up_GO_genes, RUV_ctrl_allRes6)
ctrl_up_GO_genes %>%
  dplyr::select(-ID) %>%
  column_to_rownames(var = "Term") -> ctrl_up_GO_genes

# make a heatmap
p <- pheatmap(as.matrix(ctrl_up_GO), cluster_row = FALSE, cluster_cols = FALSE,
         display_numbers = as.matrix(ctrl_up_GO_genes), fontsize = 10,
         annotation_legend = TRUE)
# save plot
png("TopGO_pheatmap_DESeq2_ctrlcoldam_ctrl.png", width = 5, height = 7, units = "in", res = 300)
p
dev.off()

##################################################
##### compare only median of ratio vs nudelta ----
# load Rdata
nospike_ctrl_allRes4 <- readRDS(file = "RData/TopGO_ctrlcoldam_nospike.RData")
nudelta_ctrl_allRes4 <- readRDS(file = "RData/TopGO_ctrlcoldam_nudelta.RData")

## combine logadj_pval
nospike_ctrl_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(MedianofRatio = logadj_pval) -> nospike_ctrl_allRes5
nudelta_ctrl_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(NuDelta = logadj_pval) -> nudelta_ctrl_allRes5
ctrl_up_GO_2norm <- full_join(nudelta_ctrl_allRes5, nospike_ctrl_allRes5)
ctrl_up_GO_2norm %>%
  dplyr::select(-ID) %>%
  column_to_rownames(var = "Term") -> ctrl_up_GO_2norm
ctrl_up_GO_2norm2 <- ctrl_up_GO_2norm
names(ctrl_up_GO_2norm2)[names(ctrl_up_GO_2norm2) == "NuDelta"] <- "SpikeIn"
names(ctrl_up_GO_2norm2)[names(ctrl_up_GO_2norm2) == "MedianofRatio"] <- "Median_of_Ratio"

## combine no. of genes
nospike_ctrl_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(MedianofRatio = Significant) -> nospike_ctrl_allRes6
nudelta_ctrl_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(NuDelta = Significant) -> nudelta_ctrl_allRes6
ctrl_up_GO_genes_2norm <- full_join(nudelta_ctrl_allRes6, nospike_ctrl_allRes6)
ctrl_up_GO_genes_2norm %>%
  dplyr::select(-ID) %>%
  column_to_rownames(var = "Term") -> ctrl_up_GO_genes_2norm
ctrl_up_GO_genes_2norm2 <- ctrl_up_GO_genes_2norm
names(ctrl_up_GO_genes_2norm2)[names(ctrl_up_GO_genes_2norm2) == "NuDelta"] <- "SpikeIn"
names(ctrl_up_GO_genes_2norm2)[names(ctrl_up_GO_genes_2norm2) == "MedianofRatio"] <- "Median_of_Ratio"

# make a heatmap
p <- pheatmap(as.matrix(ctrl_up_GO_2norm2), cluster_row = FALSE, cluster_cols = FALSE,
              display_numbers = as.matrix(ctrl_up_GO_genes_2norm2), fontsize = 12,
              annotation_legend = TRUE)
# save plot
png("TopGO_pheatmap_DESeq2_ctrlcoldam_ctrl_2norm.png", width = 5, height = 7, units = "in", res = 300)
p
dev.off()

####################################################################
####################################################################
##### chilling_am upregulated genes
### DESeq2 nospike normalization ---------------------
nospike_cold <- read.delim("output/DEGs_DESeq2_nospike_ctrlcoldam_cold.txt", header = FALSE, sep = "")
nospike_coldlist <- nospike_cold$V1
nospike_cold_geneList <- factor(as.integer(genefilteredList %in% nospike_coldlist))
names(nospike_cold_geneList) <- genefilteredList
### make a topGOdata object for BP
nospike_cold_GOdata <- new("topGOdata", description = "DESeq2 nospike cold up", 
                         ontology = "BP", allGenes = nospike_cold_geneList,
                         annot = annFUN.gene2GO, gene2GO = Sb_annotation_good)
# run Fisher's exact test
resultFisher_BP <- runTest(nospike_cold_GOdata, algorithm = "weight01", statistic = "fisher")
# get a result table
nospike_cold_allRes <- GenTable(nospike_cold_GOdata, weight01Fisher = resultFisher_BP,
                              orderBy = "weight01Fisher", topNodes = 50, numChar=1000)
View(nospike_cold_allRes)

# get genes in GO terms which weight01Fisher < 0.01 in BP
gene.universe <- genes(nospike_cold_GOdata) 
sig.genes <- sigGenes(nospike_cold_GOdata)
nospike_cold_allRes$genes <- sapply(nospike_cold_allRes$GO.ID, function(x)
{
  genes <- genesInTerm(nospike_cold_GOdata, x)
  genes[[1]][genes[[1]] %in% sig.genes] # myGenes is the queried gene list
})
nospike_cold_allRes2 <- subset(nospike_cold_allRes, 
                             as.numeric(nospike_cold_allRes[,"weight01Fisher"])<0.01 | grepl("e-", nospike_cold_allRes[,"weight01Fisher"]))
nospike_cold_allRes2$genes <-vapply(nospike_cold_allRes2$genes, paste, collapse = ",", character(1L))
nospike_cold_allRes2 %>%
  mutate(adj_pval = as.numeric(weight01Fisher)) %>%
  dplyr::rename(ID = GO.ID, Genes = genes) %>%
  mutate(category = as.factor(rep(c("BP"), nrow(nospike_cold_allRes2)))) %>%
  dplyr::select(category, ID, Term, Significant, Genes, adj_pval) -> nospike_cold_allRes3
nospike_cold_allRes3 %>%
  mutate(logadj_pval = -log(adj_pval)) -> nospike_cold_allRes4
# save RData
saveRDS(nospike_cold_allRes4, file = "RData/TopGO_ctrlcoldam_nospike_cold.RData")


####################################################################
### DESeq2 nudelta normalization ---------------------
nudelta_cold <- read.delim("output/DEGs_DESeq2_nudelta_ctrlcoldam_cold.txt", header = FALSE, sep = "")
nudelta_coldlist <- nudelta_cold$V1
nudelta_cold_geneList <- factor(as.integer(genefilteredList %in% nudelta_coldlist))
names(nudelta_cold_geneList) <- genefilteredList

### make a topGOdata object for BP
nudelta_cold_GOdata <- new("topGOdata", description = "DESeq2 nudelta cold up", 
                         ontology = "BP", allGenes = nudelta_cold_geneList,
                         annot = annFUN.gene2GO, gene2GO = Sb_annotation_good)
# run Fisher's exact test
resultFisher_BP <- runTest(nudelta_cold_GOdata, algorithm = "weight01", statistic = "fisher")
# get a result table
nudelta_cold_allRes <- GenTable(nudelta_cold_GOdata, weight01Fisher = resultFisher_BP,
                              orderBy = "weight01Fisher", topNodes = 80, numChar=1000)
View(nudelta_cold_allRes)

# get genes in GO terms which weight01Fisher < 0.01 in BP
gene.universe <- genes(nudelta_cold_GOdata) 
sig.genes <- sigGenes(nudelta_cold_GOdata)
nudelta_cold_allRes$genes <- sapply(nudelta_cold_allRes$GO.ID, function(x)
{
  genes <- genesInTerm(nudelta_cold_GOdata, x)
  genes[[1]][genes[[1]] %in% sig.genes] # myGenes is the queried gene list
})
nudelta_cold_allRes2 <- subset(nudelta_cold_allRes, 
                             as.numeric(nudelta_cold_allRes[,"weight01Fisher"])<0.01 | grepl("e-", nudelta_cold_allRes[,"weight01Fisher"]))
nudelta_cold_allRes2$genes <-vapply(nudelta_cold_allRes2$genes, paste, collapse = ",", character(1L))
nudelta_cold_allRes2 %>%
  mutate(adj_pval = as.numeric(weight01Fisher)) %>%
  dplyr::rename(ID = GO.ID, Genes = genes) %>%
  mutate(category = as.factor(rep(c("BP"), nrow(nudelta_cold_allRes2)))) %>%
  dplyr::select(category, ID, Term, Significant, Genes, adj_pval) -> nudelta_cold_allRes3
nudelta_cold_allRes3 %>%
  mutate(logadj_pval = -log(adj_pval)) -> nudelta_cold_allRes4
# save RData
saveRDS(nudelta_cold_allRes4, file = "RData/TopGO_ctrlcoldam_nudelta_cold.RData")


####################################################################
### DESeq2 spikesf normalization ---------------------
spikesf_cold <- read.delim("output/DEGs_DESeq2_spikesf_ctrlcoldam_cold.txt", header = FALSE, sep = "")
spikesf_coldlist <- spikesf_cold$V1
spikesf_cold_geneList <- factor(as.integer(genefilteredList %in% spikesf_coldlist))
names(spikesf_cold_geneList) <- genefilteredList

### make a topGOdata object for BP
spikesf_cold_GOdata <- new("topGOdata", description = "DESeq2 nudelta cold up", 
                         ontology = "BP", allGenes = spikesf_cold_geneList,
                         annot = annFUN.gene2GO, gene2GO = Sb_annotation_good)
# run Fisher's exact test
resultFisher_BP <- runTest(spikesf_cold_GOdata, algorithm = "weight01", statistic = "fisher")
# get a result table
spikesf_cold_allRes <- GenTable(spikesf_cold_GOdata, weight01Fisher = resultFisher_BP,
                              orderBy = "weight01Fisher", topNodes = 50)
View(spikesf_cold_allRes)

# get genes in GO terms which weight01Fisher < 0.001 in BP
gene.universe <- genes(spikesf_cold_GOdata) 
sig.genes <- sigGenes(spikesf_cold_GOdata)
spikesf_cold_allRes$genes <- sapply(spikesf_cold_allRes$GO.ID, function(x)
{
  genes <- genesInTerm(spikesf_cold_GOdata, x)
  genes[[1]][genes[[1]] %in% sig.genes] # myGenes is the queried gene list
})
spikesf_cold_allRes2 <- subset(spikesf_cold_allRes, 
                             as.numeric(spikesf_cold_allRes[,"weight01Fisher"])<0.01 | grepl("e-", spikesf_cold_allRes[,"weight01Fisher"]))
spikesf_cold_allRes2$genes <-vapply(spikesf_cold_allRes2$genes, paste, collapse = ",", character(1L))
spikesf_cold_allRes2 %>%
  mutate(adj_pval = as.numeric(weight01Fisher)) %>%
  dplyr::rename(ID = GO.ID, Genes = genes) %>%
  mutate(category = as.factor(rep(c("BP"), nrow(spikesf_cold_allRes2)))) %>%
  dplyr::select(category, ID, Term, Significant, Genes, adj_pval) -> spikesf_cold_allRes3
spikesf_cold_allRes3 %>%
  mutate(logadj_pval = -log(adj_pval)) -> spikesf_cold_allRes4


####################################################################
### DESeq2 RUV normalization ---------------------
RUV_cold <- read.delim("output/DEGs_DESeq2_RUV_ctrlcoldam_cold.txt", header = FALSE, sep = "")
RUV_coldlist <- RUV_cold$V1
RUV_cold_geneList <- factor(as.integer(genefilteredList %in% RUV_coldlist))
names(RUV_cold_geneList) <- genefilteredList
### make a topGOdata object for BP
RUV_cold_GOdata <- new("topGOdata", description = "DESeq2 RUV cold up", 
                         ontology = "BP", allGenes = RUV_cold_geneList,
                         annot = annFUN.gene2GO, gene2GO = Sb_annotation_good)
# run Fisher's exact test
resultFisher_BP <- runTest(RUV_cold_GOdata, algorithm = "weight01", statistic = "fisher")
# get a result table
RUV_cold_allRes <- GenTable(RUV_cold_GOdata, weight01Fisher = resultFisher_BP,
                              orderBy = "weight01Fisher", topNodes = 50)
View(RUV_cold_allRes)

# get genes in GO terms which weight01Fisher < 0.001 in BP
gene.universe <- genes(RUV_cold_GOdata) 
sig.genes <- sigGenes(RUV_cold_GOdata)
RUV_cold_allRes$genes <- sapply(RUV_cold_allRes$GO.ID, function(x)
{
  genes <- genesInTerm(RUV_cold_GOdata, x)
  genes[[1]][genes[[1]] %in% sig.genes] # myGenes is the queried gene list
})
RUV_cold_allRes2 <- subset(RUV_cold_allRes, 
                             as.numeric(RUV_cold_allRes[,"weight01Fisher"])<0.01 | grepl("e-", RUV_cold_allRes[,"weight01Fisher"]))
RUV_cold_allRes2$genes <-vapply(RUV_cold_allRes2$genes, paste, collapse = ",", character(1L))
RUV_cold_allRes2 %>%
  mutate(adj_pval = as.numeric(weight01Fisher)) %>%
  dplyr::rename(ID = GO.ID, Genes = genes) %>%
  mutate(category = as.factor(rep(c("BP"), nrow(RUV_cold_allRes2)))) %>%
  dplyr::select(category, ID, Term, Significant, Genes, adj_pval) -> RUV_cold_allRes3
RUV_cold_allRes3 %>%
  mutate(logadj_pval = -log(adj_pval)) -> RUV_cold_allRes4


############################################################
##### combine GO terms from four normalization methods ----
## combine logadj_pval
nospike_cold_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(MedianofRatio = logadj_pval) -> nospike_cold_allRes5
nudelta_cold_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(NuDelta = logadj_pval) -> nudelta_cold_allRes5
spikesf_cold_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(SpikeSizeFactor = logadj_pval) -> spikesf_cold_allRes5
RUV_cold_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(RUV = logadj_pval) -> RUV_cold_allRes5
cold_up_GO <- full_join(nospike_cold_allRes5, nudelta_cold_allRes5)
cold_up_GO <- full_join(cold_up_GO, spikesf_cold_allRes5)
cold_up_GO <- full_join(cold_up_GO, RUV_cold_allRes5)
cold_up_GO %>%
  dplyr::select(-ID) %>%
  column_to_rownames(var = "Term") -> cold_up_GO

## combine no. of genes
nospike_cold_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(MedianofRatio = Significant) -> nospike_cold_allRes6
nudelta_cold_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(NuDelta = Significant) -> nudelta_cold_allRes6
spikesf_cold_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(SpikeSizeFactor = Significant) -> spikesf_cold_allRes6
RUV_cold_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(RUV = Significant) -> RUV_cold_allRes6
cold_up_GO_genes <- full_join(nospike_cold_allRes6, nudelta_cold_allRes6)
cold_up_GO_genes <- full_join(cold_up_GO_genes, spikesf_cold_allRes6)
cold_up_GO_genes <- full_join(cold_up_GO_genes, RUV_cold_allRes6)
cold_up_GO_genes %>%
  dplyr::select(-ID) %>%
  column_to_rownames(var = "Term") -> cold_up_GO_genes

# make a heatmap
p <- pheatmap(as.matrix(cold_up_GO), cluster_row = FALSE, cluster_cols = FALSE,
              display_numbers = as.matrix(cold_up_GO_genes), fontsize = 12,
              annotation_legend = TRUE)

# save plot
png("TopGO_pheatmap_DESeq2_ctrlcoldam_cold.png", width = 6, height = 8, units = "in", res = 300)
p
dev.off()

###############################################
##### compare only median of ratio vs nudelta ----
# load RData
nospike_cold_allRes4 <- readRDS(file = "RData/TopGO_ctrlcoldam_nospike_cold.RData")
nudelta_cold_allRes4 <- readRDS(file = "RData/TopGO_ctrlcoldam_nudelta_cold.RData")

## combine logadj_pval
nospike_cold_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(MedianofRatio = logadj_pval) -> nospike_cold_allRes5
nudelta_cold_allRes4 %>%
  dplyr::select(ID, Term, logadj_pval) %>%
  dplyr::rename(NuDelta = logadj_pval) -> nudelta_cold_allRes5
cold_up_GO_2norm <- full_join(nudelta_cold_allRes5, nospike_cold_allRes5)
cold_up_GO_2norm %>%
  dplyr::select(-ID) %>%
  column_to_rownames(var = "Term") -> cold_up_GO_2norm 
cold_up_GO_2norm2 <- cold_up_GO_2norm
names(cold_up_GO_2norm2)[names(cold_up_GO_2norm2) == "NuDelta"] <- "SpikeIn"
names(cold_up_GO_2norm2)[names(cold_up_GO_2norm2) == "MedianofRatio"] <- "Median_of_Ratio"

## combine no. of genes
nospike_cold_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(MedianofRatio = Significant) -> nospike_cold_allRes6
nudelta_cold_allRes4 %>%
  dplyr::select(ID, Term, Significant) %>%
  dplyr::rename(NuDelta = Significant) -> nudelta_cold_allRes6
cold_up_GO_genes_2norm <- full_join(nudelta_cold_allRes6, nospike_cold_allRes6)
cold_up_GO_genes_2norm %>%
  dplyr::select(-ID) %>%
  column_to_rownames(var = "Term") -> cold_up_GO_genes_2norm
cold_up_GO_genes_2norm2 <- cold_up_GO_genes_2norm
names(cold_up_GO_genes_2norm2)[names(cold_up_GO_genes_2norm2) == "NuDelta"] <- "SpikeIn"
names(cold_up_GO_genes_2norm2)[names(cold_up_GO_genes_2norm2) == "MedianofRatio"] <- "Median_of_Ratio"

# make a heatmap
p <- pheatmap(as.matrix(cold_up_GO_2norm2), cluster_row = FALSE, cluster_cols = FALSE,
              display_numbers = as.matrix(cold_up_GO_genes_2norm2), fontsize = 12,
              annotation_legend = TRUE)
# save plot
png("TopGO_pheatmap_DESeq2_ctrlcoldam_cold_2norm.png", width = 8, height = 12, units = "in", res = 300)
p
dev.off()


