# mb2021_phyloseq
Different methods for analyzing microbiome data in a phyloseq data format

Here, data is in phyloseq format - "Marissa_Osyter" (Oyster spelled wrong - I should fix that one day lol)

Description of each script: 

1. PCoA_MDS script - creating pseq and pseq.rel objects, generating composition plots and PCOA plots. Still need to learn how to run PERMANOVA tests (if possible when data is in phyloseq format).
2. Plotting abundance phyloseq script - For plotting top # of microbiota and/or plotting specific microbiota. Code unable to filter out samples. Code is unable to generate boxplots (phyloseq objects cannot be coerced into data frame). Code can plot relative or total abundance.
3. mb2021 network analysis - Attempting to generate code for network analysis using phyloseq microbiome data.
4. Correlations analysis heatmaps - Running Spearman's or Pearson's correlations with phyloseq microbiome data and visualizing with heatmaps.
5. mb2021 3-way ANOVA - Using raw data files (NOT phyloseq data), can test 3-way interactions (in this case, the 3 factors are Treatment (Control, HS, LS), Genetics (F1-4), and Age (day 1 - spat). Tested with Flavobactereceae, Rhodobactereceae, and Vibrionaceae. 
