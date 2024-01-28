##### script for making a plot to compare number of DEGs from different normalization methods

# clear workspace
rm(list = ls())
# load package
library(readxl)
library(tidyverse)
library(lemon)
library(ggpubr)
library(ggtext)
library(glue)
# set working directory
setwd("~/sorghum_spikein")

#################################################################
##### import data --------------------------------------------------------------
Sb_DEG <- read_excel("Input/sorghum_number of DEGs_new.xlsx")
Sb_DEG <- as.data.frame(Sb_DEG)
Sb_DEG$norm_methods <- as.factor(Sb_DEG$norm_methods)
Sb_DEG$norm_methods  <- factor(Sb_DEG$norm_methods, 
                               levels = c("DESeq2_MedianofRatio", "DESeq2_RUV",
                                          "DESeq2_NuDelta", "DESeq2_SpikeSizeFactor", 
                                          "EdgeR_TMM", "EdgeR_RUV", "EdgeR_log2spike"))

#################################################################
### select controlam_pm data ----
Sb_DEG %>%
  filter(condition == "controlam_pm") %>%
  filter(norm_methods != "EdgeR_Nudelta") %>%
  dplyr::select(DE_tool, norm_methods, AM_up, PM_up) %>%
  gather(key = TOD, value = number, -DE_tool, -norm_methods) -> Sb_DEG_control

ggplot(Sb_DEG_control, aes(x = norm_methods, y = number, fill = TOD)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  xlab("DE Tools and Normalization Methods") +
  ylab("Number of DEGs") +
  labs(fill = "") +
  theme_minimal()

### make a back-to-back barchart
Sb_DEG_control$norm_methods <- factor(Sb_DEG_control$norm_methods, 
                                      levels = rev(levels(Sb_DEG_control$norm_methods)))
Sb_DEG_control$StyledClass <- c("<i>DESeq2 Median of Ratio</i>", "<b><i>DESeq2 Spike-in (Nu*Delta)</i></b>",
                                "<b><i>DESeq2 Spike-in size factor</i><b>",  "<i>DESeq2 RUV</i>",
                                "EdgeR TMM", "<b>EdgeR log2 spike-in</b>", "EdgeR RUV")
Sb_DEG_control$StyledClass <- fct_reorder(Sb_DEG_control$StyledClass, 
                                          as.integer(Sb_DEG_control$norm_methods))
Sb_DEG_control_plot <- 
  ggplot(data = Sb_DEG_control, 
       mapping = aes(x = ifelse(test = TOD == "AM_up", yes = -number, no = number), 
                       y = StyledClass, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  labs(x = "Number of DEGs", y = "Normalization methods",
       title = "Control_AM vs Control_PM", fill = "") +
  theme_minimal() +
  theme(axis.text.y = ggtext::element_markdown(),
        plot.title = element_text(size = 14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
# save plot
ggsave(filename = "control_ampm_DEG.jpeg", plot = Sb_DEG_control_plot,
       width = 6, height = 4, dpi = 300, units = "in", device='png')

### plot only DESeq2 Data ----
Sb_DEG_control %>%
  filter(DE_tool == "DESeq2") -> Sb_DEG_control_DESeq2
Sb_DEG_control_DESeq2$norm_methods <- gsub("DESeq2_", "", as.character(Sb_DEG_control_DESeq2$norm_methods))
Sb_DEG_control_DESeq2$norm_methods <- factor(Sb_DEG_control_DESeq2$norm_methods, 
                                             levels = c("RUV", "SpikeSizeFactor", "NuDelta", "MedianofRatio"))

ggplot(data = Sb_DEG_control_DESeq2, 
       mapping = aes(x = ifelse(test = TOD == "AM_up", yes = -number, no = number), 
                     y = norm_methods, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  labs(x = "Number of DEGs", y = "Normalization methods",
       title = "Control_AM vs Control_PM") +
  scale_y_discrete(breaks = c("RUV", "SpikeSizeFactor", "NuDelta", "MedianofRatio"), 
                   labels = c("RUV", "Spike-in size factor", "Nu*Delta", "Median of ratio")) +
  theme_minimal()
# save plot
ggsave(filename = "control_ampm_DEG_DESeq2.jpeg", plot = last_plot(),
       width = 6, height = 4, dpi = 300, units = "in", device='png')

### plot only DESeq2: median of ratio and nudelta ----
Sb_DEG_control_DESeq2 %>%
  filter(norm_methods == "MedianofRatio" | norm_methods == "NuDelta") -> Sb_DEG_control_DESeq2_n

ggplot(data = Sb_DEG_control_DESeq2_n, 
       mapping = aes(x = ifelse(test = TOD == "AM_up", yes = -number, no = number), 
                     y = norm_methods, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  labs(x = "Number of DEGs", y = "Normalization methods",
       title = "Control_AM vs Control_PM", fill = "") +
  scale_y_discrete(breaks = c("NuDelta", "MedianofRatio"), 
                   labels = c("Spike-in", "Median of ratio")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
# save plot
ggsave(filename = "control_ampm_DEG_DESeq2_n.jpeg", plot = last_plot(),
       width = 5, height = 2, dpi = 300, units = "in", device='png')

#################################################################
### select coldam_pm data ----
Sb_DEG %>%
  filter(condition == "coldam_pm") %>%
  filter(norm_methods != "EdgeR_Nudelta") %>%
  dplyr::select(DE_tool, norm_methods, AM_up, PM_up) %>%
  gather(key = TOD, value = number, -DE_tool, -norm_methods) -> Sb_DEG_cold

ggplot(Sb_DEG_cold, aes(x = norm_methods, y = number, fill = TOD)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  xlab("DE Tools and Normalization Methods") +
  ylab("Number of DEGs") +
  labs(fill = "") +
  theme_minimal()

### make a back-to-back barchart -----
Sb_DEG_cold$norm_methods <- factor(Sb_DEG_cold$norm_methods, 
                                      levels = rev(levels(Sb_DEG_cold$norm_methods)))
Sb_DEG_cold$StyledClass <- c("<i>DESeq2 Median of Ratio</i>", "<b><i>DESeq2 Spike-in (Nu*Delta)</i></b>",
                             "<b><i>DESeq2 Spike-in size factor</i><b>",  "<i>DESeq2 RUV</i>",
                             "EdgeR TMM", "<b>EdgeR log2 spike-in</b>", "EdgeR RUV")
Sb_DEG_cold$StyledClass <- fct_reorder(Sb_DEG_cold$StyledClass, 
                                          as.integer(Sb_DEG_cold$norm_methods))
Sb_DEG_cold_plot <- 
  ggplot(data = Sb_DEG_cold, 
       mapping = aes(x = ifelse(test = TOD == "AM_up", yes = -number, no = number), 
                     y = StyledClass, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  labs(x = "Number of DEGs", y = "Normalization methods", 
       title = "Chilling_AM vs Chilling_PM", fill = "") +
  theme_minimal() +
  theme(axis.text.y = ggtext::element_markdown(),
        plot.title = element_text(size = 14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
# save plot
ggsave(filename = "chilling_ampm_DEG.jpeg", plot = Sb_DEG_cold_plot,
       width = 6, height = 4, dpi = 300, units = "in", device='png')

### plot only DESeq2 Data ----
Sb_DEG_cold %>%
  filter(DE_tool == "DESeq2") -> Sb_DEG_cold_DESeq2
Sb_DEG_cold_DESeq2$norm_methods <- gsub("DESeq2_", "", as.character(Sb_DEG_cold_DESeq2$norm_methods))
Sb_DEG_cold_DESeq2$norm_methods <- factor(Sb_DEG_cold_DESeq2$norm_methods, 
                                          levels = c("RUV", "SpikeSizeFactor", "NuDelta", 
                                                     "MedianofRatio"))

ggplot(data = Sb_DEG_cold_DESeq2, 
       mapping = aes(x = ifelse(test = TOD == "AM_up", yes = -number, no = number), 
                     y = norm_methods, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  labs(x = "Number of DEGs", y = "Normalization methods",
       title = "Chilling_AM vs Chilling_PM") +
  scale_y_discrete(breaks = c("RUV", "SpikeSizeFactor", "NuDelta", "MedianofRatio"), 
                   labels = c("RUV", "Spike-in size factor", "Nu*Delta", "Median of ratio")) +
  theme_minimal()
# save plot
ggsave(filename = "chilling_ampm_DEG_DESeq2.jpeg", plot = last_plot(),
       width = 6, height = 4, dpi = 300, units = "in", device='png')

### plot only DESeq2: median of ratio and nudelta ----
Sb_DEG_cold_DESeq2 %>%
  filter(norm_methods == "MedianofRatio" | norm_methods == "NuDelta") -> Sb_DEG_cold_DESeq2_n

ggplot(data = Sb_DEG_cold_DESeq2_n, 
       mapping = aes(x = ifelse(test = TOD == "AM_up", yes = -number, no = number), 
                     y = norm_methods, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  labs(x = "Number of DEGs", y = "Normalization methods",
       title = "Chilling_AM vs Chilling_PM", fill = "") +
  scale_y_discrete(breaks = c("NuDelta", "MedianofRatio"), 
                   labels = c("Spike-in", "Median of ratio")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
# save plot
ggsave(filename = "cold_ampm_DEG_DESeq2_n.jpeg", plot = last_plot(),
       width = 5, height = 2, dpi = 300, units = "in", device='png')

###################################################
### combine figure 
Sb_DEG_plot <- ggarrange(Sb_DEG_control_plot, Sb_DEG_cold_plot, 
                           ncol = 2, common.legend = TRUE, legend = "bottom")
# save plot
ggsave(filename = "number_DEGs_plot.png", plot = Sb_DEG_plot, width = 12, height = 5, 
       dpi = 300, units = "in", device='png')

#################################################################
### select control_cold_am data ----
Sb_DEG %>%
  filter(condition == "control_cold_am") %>%
  filter(norm_methods != "EdgeR_Nudelta") %>%
  dplyr::select(DE_tool, norm_methods, AM_up, PM_up) %>%
  gather(key = TOD, value = number, -DE_tool, -norm_methods) -> Sb_DEG_ctrlcold_am
Sb_DEG_ctrlcold_am$TOD[Sb_DEG_ctrlcold_am$TOD == "AM_up"] <- "control_up"
Sb_DEG_ctrlcold_am$TOD[Sb_DEG_ctrlcold_am$TOD == "PM_up"] <- "chilling_up"

ggplot(Sb_DEG_ctrlcold_am, aes(x = norm_methods, y = number, fill = TOD)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  xlab("DE Tools and Normalization Methods") +
  ylab("Number of DEGs") +
  labs(fill = "") +
  theme_minimal()

### make a back-to-back barchart ----
Sb_DEG_ctrlcold_am$norm_methods <- factor(Sb_DEG_ctrlcold_am$norm_methods, 
                                          levels = rev(levels(Sb_DEG_ctrlcold_am$norm_methods)))
Sb_DEG_ctrlcold_am$StyledClass <- c("<i>DESeq2 Median of Ratio</i>", "<b><i>DESeq2 Spike-in (Nu*Delta)</i></b>",
                                    "<b><i>DESeq2 Spike-in size factor</i><b>",  "<i>DESeq2 RUV</i>",
                                    "EdgeR TMM", "<b>EdgeR log2 spike-in</b>", "EdgeR RUV")
Sb_DEG_ctrlcold_am$StyledClass <- fct_reorder(Sb_DEG_ctrlcold_am$StyledClass, 
                                          as.integer(Sb_DEG_ctrlcold_am$norm_methods))

Sb_DEG_ctrlcold_am_plot <- 
  ggplot(data = Sb_DEG_ctrlcold_am, 
         mapping = aes(x = ifelse(test = TOD == "control_up", yes = -number, no = number), 
                       y = StyledClass, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"),
                      breaks = c("control_up", "chilling_up"),
                      labels = c("chilling_down", "chilling_up")) +
  labs(x = "Number of DEGs", y = "Normalization methods", 
       title = "Control_AM vs Chilling_AM", fill = "") +
  theme_minimal() +
  theme(axis.text.y = ggtext::element_markdown(),
        plot.title = element_text(size = 14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
# save plot
ggsave(filename = "ctrl_cold_am_DEG.jpeg", plot = Sb_DEG_ctrlcold_am_plot,
       width = 6, height = 4, dpi = 300, units = "in", device='png')

### plot only DESeq2 Data ----
Sb_DEG_ctrlcold_am %>%
  filter(DE_tool == "DESeq2") -> Sb_DEG_ctrlcold_am_DESeq2
Sb_DEG_ctrlcold_am_DESeq2$norm_methods <- gsub("DESeq2_", "", 
                                               as.character(Sb_DEG_ctrlcold_am_DESeq2$norm_methods))
Sb_DEG_ctrlcold_am_DESeq2$norm_methods <- factor(Sb_DEG_ctrlcold_am_DESeq2$norm_methods, 
                                                 levels = c("RUV", "SpikeSizeFactor", "NuDelta", 
                                                            "MedianofRatio"))

ggplot(data = Sb_DEG_ctrlcold_am_DESeq2, 
       mapping = aes(x = ifelse(test = TOD == "control_up", yes = -number, no = number), 
                     y = norm_methods, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"),
                    breaks = c("control_up", "chilling_up"),
                    labels = c("chilling_down", "chilling_up")) +
  labs(x = "Number of DEGs", y = "Normalization methods",
       title = "Control_AM vs Chilling_AM") +
  scale_y_discrete(breaks = c("RUV", "SpikeSizeFactor", "NuDelta", "MedianofRatio"), 
                   labels = c("RUV", "Spike-in size factor", "Nu*Delta", "Median of ratio")) +
  theme_minimal() +
  theme(legend.title = element_blank())
# save plot
ggsave(filename = "ctrl_cold_am_DEG_DESeq2.jpeg", plot = last_plot(),
       width = 6, height = 4, dpi = 300, units = "in", device='png')

### plot only DESeq2: median of ratio and nudelta ----
Sb_DEG_ctrlcold_am_DESeq2 %>%
  filter(norm_methods == "MedianofRatio" | norm_methods == "NuDelta") -> Sb_DEG_ctrlcold_am_DESeq2_n

ggplot(data = Sb_DEG_ctrlcold_am_DESeq2_n, 
       mapping = aes(x = ifelse(test = TOD == "control_up", yes = -number, no = number), 
                     y = norm_methods, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"),
                    breaks = c("control_up", "chilling_up"),
                    labels = c("chilling_down", "chilling_up")) +
  labs(x = "Number of DEGs", y = "Normalization methods",
       title = "Control_AM vs Chilling_AM", fill = "") +
  scale_y_discrete(breaks = c("NuDelta", "MedianofRatio"), 
                   labels = c("Spike-in", "Median of ratio")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
# save plot
ggsave(filename = "ctrl_cold_am_DEG_DESeq2_n.jpeg", plot = last_plot(),
       width = 5, height = 2, dpi = 300, units = "in", device='png')

#################################################################
### select control_cold_pm data ----
Sb_DEG %>%
  filter(condition == "control_cold_pm") %>%
  filter(norm_methods != "EdgeR_Nudelta") %>%
  dplyr::select(DE_tool, norm_methods, AM_up, PM_up) %>%
  gather(key = TOD, value = number, -DE_tool, -norm_methods) -> Sb_DEG_ctrlcold_pm
Sb_DEG_ctrlcold_pm$TOD[Sb_DEG_ctrlcold_pm$TOD == "AM_up"] <- "control_up"
Sb_DEG_ctrlcold_pm$TOD[Sb_DEG_ctrlcold_pm$TOD == "PM_up"] <- "chilling_up"

ggplot(Sb_DEG_ctrlcold_pm, aes(x = norm_methods, y = number, fill = TOD)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  xlab("DE Tools and Normalization Methods") +
  ylab("Number of DEGs") +
  labs(fill = "") +
  theme_minimal()

### make a back-to-back barchart ----
Sb_DEG_ctrlcold_pm$norm_methods <- factor(Sb_DEG_ctrlcold_pm$norm_methods, 
                                          levels = rev(levels(Sb_DEG_ctrlcold_pm$norm_methods)))
Sb_DEG_ctrlcold_pm$StyledClass <- c("<i>DESeq2 Median of Ratio</i>", "<b><i>DESeq2 Spike-in (Nu*Delta)</i></b>",
                                    "<b><i>DESeq2 Spike-in size factor</i><b>",  "<i>DESeq2 RUV</i>",
                                    "EdgeR TMM", "<b>EdgeR log2 spike-in</b>", "EdgeR RUV")
Sb_DEG_ctrlcold_pm$StyledClass <- fct_reorder(Sb_DEG_ctrlcold_pm$StyledClass, 
                                          as.integer(Sb_DEG_ctrlcold_pm$norm_methods))
Sb_DEG_ctrlcold_pm_plot <- 
  ggplot(data = Sb_DEG_ctrlcold_pm, 
         mapping = aes(x = ifelse(test = TOD == "control_up", yes = -number, no = number), 
                       y = StyledClass, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"),
                    breaks = c("control_up", "chilling_up"),
                    labels = c("chilling_down", "chilling_up")) +
  labs(x = "Number of DEGs", y = "Normalization methods", 
       title = "Control_PM vs Chilling_PM", fill = "") +
  theme_minimal() +
  theme(axis.text.y = ggtext::element_markdown(),
        plot.title = element_text(size = 14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
# save plot
ggsave(filename = "ctrl_cold_pm_DEG.jpeg", plot = Sb_DEG_ctrlcold_pm_plot,
       width = 6, height = 4, dpi = 300, units = "in", device='png')

### plot only DESeq2 Data ----
Sb_DEG_ctrlcold_pm %>%
  filter(DE_tool == "DESeq2") -> Sb_DEG_ctrlcold_pm_DESeq2
Sb_DEG_ctrlcold_pm_DESeq2$norm_methods <- gsub("DESeq2_", "", 
                                               as.character(Sb_DEG_ctrlcold_pm_DESeq2$norm_methods))
Sb_DEG_ctrlcold_pm_DESeq2$norm_methods <- factor(Sb_DEG_ctrlcold_pm_DESeq2$norm_methods, 
                                                 levels = c("RUV", "SpikeSizeFactor", "NuDelta", 
                                                            "MedianofRatio"))
ggplot(data = Sb_DEG_ctrlcold_pm_DESeq2, 
       mapping = aes(x = ifelse(test = TOD == "control_up", yes = -number, no = number), 
                     y = norm_methods, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"),
                    breaks = c("control_up", "chilling_up"),
                    labels = c("chilling_down", "chilling_up")) +
  labs(x = "Number of DEGs", y = "Normalization methods",
       title = "Control_PM vs Chilling_PM") +
  scale_y_discrete(breaks = c("RUV", "SpikeSizeFactor", "NuDelta", "MedianofRatio"), 
                   labels = c("RUV", "Spike-in size factor", "Nu*Delta", "Median of ratio")) +
  theme_minimal() +
  theme(legend.title = element_blank())
# save plot
ggsave(filename = "ctrl_cold_pm_DEG_DESeq2.jpeg", plot = last_plot(),
       width = 6, height = 4, dpi = 300, units = "in", device='png')

### plot only DESeq2: median of ratio and nudelta ----
Sb_DEG_ctrlcold_pm_DESeq2 %>%
  filter(norm_methods == "MedianofRatio" | norm_methods == "NuDelta") -> Sb_DEG_ctrlcold_pm_DESeq2_n

ggplot(data = Sb_DEG_ctrlcold_pm_DESeq2_n, 
       mapping = aes(x = ifelse(test = TOD == "control_up", yes = -number, no = number), 
                     y = norm_methods, fill = TOD, label = number)) +
  geom_col() +
  geom_text(size = 4, position = position_stack(vjust = 0.5)) +
  scale_x_symmetric(labels = abs) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"),
                    breaks = c("control_up", "chilling_up"),
                    labels = c("chilling_down", "chilling_up")) +
  labs(x = "Number of DEGs", y = "Normalization methods",
       title = "Control_PM vs Chilling_PM", fill = "") +
  scale_y_discrete(breaks = c("NuDelta", "MedianofRatio"), 
                   labels = c("Spike-in", "Median of ratio")) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
# save plot
ggsave(filename = "ctrl_cold_pm_DEG_DESeq2_n.jpeg", plot = last_plot(),
       width = 5, height = 2, dpi = 300, units = "in", device='png')

###################################################
### combine figure 
Sb_DEG_plot2 <- ggarrange(Sb_DEG_ctrlcold_am_plot, Sb_DEG_ctrlcold_pm_plot, 
                         ncol = 2, common.legend = TRUE, legend = "bottom")
# save plot
ggsave(filename = "number_DEGs_plot2.png", plot = Sb_DEG_plot2, width = 12, height = 5, 
       dpi = 300, units = "in", device='png')

