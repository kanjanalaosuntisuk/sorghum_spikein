# sorghum_spikein

This is a repository of data and scripts used in the manuscript "A normalization method that controls for total RNA abundance affects the identification of differentially expressed genes, revealing bias toward morning-expressed responses".

Data used in the analysis are in the 'RData' folder, and R scripts are in the 'scripts' folder.
**Note: please do not source these codes directly as many of the commands are not wrapped in functions**

**Scripts**
- BTx_gene_count_table.R: to tidy gene count table
- BTx_spike_count_table.R: to tidy and analyze spike-in count table
- BTx_DESeq2_MedianofRatio.R: to perform traditional normalization in DESeq2
- BTx_DESeq2_spike_calibration.R: to perform spike-in normalization using a method from Athanasiadou _et al._ (2019) _PLoS Computational Biology_
- BTx_DESeq2_RUV.R: to perform spike-in normalization using a method from Risso _et al._ (2014) _Nature Biotechnology_
- BTx_DESeq2_spike_sizefactor.R: to perform spike-in normalization using a method from Brennecke _et al._ (2013) _Nature Methods_
- BTx_number_DEGs_compare.R: to make plots comparing the numbers of DEGs from different normalization methods
- BTx_compare_DEGs_2norm_TOD.R: to make heat maps and Venn diagrams comparing AM-vs-PM DEGs between traditional and spike-in normalization methods
- BTx_compare_DEGs_2norm_cold.R: to make heat maps and Venn diagrams comparing control-vs-chilling DEGs between traditional and spike-in normalization methods
- DESeq2_compare_gene_exp_2norm_MapMan.R: to make dot plots of selected AM-vs-PM genes from MapMan analysis
- DESeq2_compare_gene_exp_2norm_cold_MapMan.R: to make dot plots of selected control-vs-chilling DEGs from MapMan analysis
- DESeq2_compare_gene_exp_2norm.R: to make dot plots of selected genes to demonstrate differences in normalized expression levels between normalization methods
- DESeq2_compare_gene_exp_2norm_cold.R: to make dot plots of selected genes to demonstrate differences in normalized expression levels between normalization methods
- TopGO_DESeq2_ctrlcoldam.R: to perform GO analysis and make heat maps of significant GO terms using AM control-vs-chilling DEGs
- TopGO_DESeq2_ctrlcoldpm.R: to perform GO analysis and make heat maps of significant GO terms using PM control-vs-chilling DEGs

