#Microbiotaprocesses package
#source: https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html#1-load-required-packages


Marissa_MU42022_rarefied_20231016 <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rarefied_20231016.rds")

pseq <-  Marissa_MU42022_rarefied_20231016

#install ----
install.packages("MicrobiotaProcess")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("MicrobiotaProcess")

install.packages("coin")
install.packages("reshape2")
install.packages("ggnewscale")


#load
suppressPackageStartupMessages({
  library(MicrobiotaProcess) # an R package for analysis, visualization and biomarker discovery of Microbiome.
  library(phyloseq) # Handling and analysis of high-throughput microbiome census data.
  library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics.
  library(tidyverse) # Easily Install and Load the 'Tidyverse'.
  library(vegan) # Community Ecology Package.
  library(coin) # Conditional Inference Procedures in a Permutation Test Framework.
  library(reshape2) # Flexibly Reshape Data: A Reboot of the Reshape Package.
  library(ggnewscale) # Multiple Fill and Colour Scales in 'ggplot2'.
})


#biomarker discovery ----

# for the kruskal_test and wilcox_test
library(coin)
library(MicrobiotaProcess)
# Since the effect size was calculated by randomly re-sampling, 
# the seed should be set for reproducibly results.
set.seed(1024)
deres <- diff_analysis(obj = pseq, classgroup = "",
                       mlfun = "lda",
                       filtermod = "pvalue",
                       firstcomfun = "kruskal_test",
                       firstalpha = 0.05,
                       strictmod = TRUE,
                       secondcomfun = "wilcox_test",
                       subclmin = 3,
                       subclwilc = TRUE,
                       secondalpha = 0.01,
                       lda=3)
deres

diffres <- diff_analysis(pseq, classgroup="Treatment", taxda = "", 
                         mlfun="lda", filtermod="fdr",
                         firstcomfun = "kruskal.test",
                         firstalpha=0.05, strictmod=TRUE,
                         secondcomfun = "wilcox.test",
                         subclmin=3, subclwilc=TRUE,
                         secondalpha=0.01, ldascore=3)


deres <- diff_analysis(obj = pseq, 
                      classgroup = "Treatment",
                      mlfun = "mda",
                      filtermod = "pvalue",
                      firstcomfun = "kruskal_test",
                      firstalpha = 0.05,
                      strictmod = TRUE,
                      secondcomfun = "wilcox_test",
                      subclmin = 3,
                      subclwilc = TRUE,
                      secondalpha = 0.01,
                      mda = 3)



#> The original data: 689 features and 41 samples
#> The sample data: 1 variables and 41 samples
#> The taxda contained 1432 by 7 rank
#> after first test (kruskal_test) number of feature (pvalue<=0.05):71
#> after second test (wilcox_test) number of significantly discriminative feature:28
#> after lda, Number of discriminative features: 22 (certain taxonomy classification:15; uncertain taxonomy classication: 7)

#visualize

diffbox <- ggdiffbox(obj=deres, box_notch=FALSE, 
                     colorlist=c("#00AED7", "#FD9347"), l_xlabtext="relative abundance")
diffbox

ggdifftaxbar(obj=deres, xtextsize=1.5, 
             output="IBD_biomarkder_barplot",
             coloslist=c("#00AED7", "#FD9347"))

es_p <- ggeffectsize(obj=deres, 
                     lineheight=0.1,
                     linewidth=0.3) + 
  scale_color_manual(values=c("#00AED7", 
                              "#FD9347")) 

es_p