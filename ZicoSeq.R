#Zicoseq differential abundace test
#see paper: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-022-01320-0
#Supposed to control for FDR but keep higher power (unlike ANCOMBC and ALDEx2)
#See package info: https://cran.r-project.org/web/packages/GUniFrac/vignettes/ZicoSeq.html

#load
install.packages("GUniFrac")
library(GUniFrac)

#pseq object - data

#may want to use CSS normalized 
pseq <- PB2023_spat_not_rarefied_CSSnormalized_Jan2025
pseq <- PB2023_spat_filtered_not_rarefied
pseq <- subset_samples(pseq, !Treatment %in% c("James", "Continuous Probiotics"))

#Quality check
any(taxa_sums(pseq) == 0)

#if true

pseq_filtered <- prune_taxa(taxa_sums(pseq) > 0, pseq)
any(taxa_sums(pseq_filtered) == 0)

pseq <- pseq_filtered

OTU <- pseq@otu_table
Meta <- pseq@sam_data

comm <- t(OTU)

meta.dat = sample_data(pseq)
meta = as.matrix.data.frame(meta.dat)
meta.dat = as.data.frame(meta)
str(meta.dat)

meta.dat$Treatment <- as.factor(meta.dat$Treatment)  # Convert to factor

ZicoSeq.obj <- ZicoSeq(meta.dat = meta.dat, feature.dat = comm, 
                       grp.name = 'Treatment', feature.dat.type = "count",
                       # Filter to remove rare taxa
                       prev.filter = 0.2, mean.abund.filter = 0,  
                       max.abund.filter = 0.002, min.prop = 0, 
                       # Winsorization to replace outliers
                       is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                       # Posterior sampling 
                       is.post.sample = FALSE, post.sample.no = 25, 
                       # Use the square-root transformation
                       link.func = list(function (x) x^0.5), stats.combine.func = max,
                       # Permutation-based multiple testing correction
                       perm.no = 99,  strata = NULL, 
                       # Reference-based multiple stage normalization
                       ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                       # Family-wise error rate control
                       is.fwer = FALSE, verbose = TRUE, return.feature.dat = TRUE)

#For sample size less than 40, posterior sampling will not be used!
# 53  features are filtered!
#   The data has  17  samples and  130  features will be tested!
#   On average,  1  outlier counts will be replaced for each feature!

names(ZicoSeq.obj)

# Raw p-values
raw_p_values <- ZicoSeq.obj$p.raw

# Adjusted p-values for FDR
adj_p_values_fdr <- ZicoSeq.obj$p.adj.fdr

# Adjusted p-values for FWER
adj_p_values_fwer <- ZicoSeq.obj$p.adj.fwer


ZicoSeq.plot(ZicoSeq.obj, pvalue.type = 'p.adj.fdr', cutoff = 0.1, text.size = 10,
             out.dir = NULL, width = 10, height = 6)
