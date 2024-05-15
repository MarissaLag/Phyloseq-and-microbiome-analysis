
#selecting pca axis using supervised method

library("phyloseq")

library("ggplot2")

library("plyr")



theme_set(theme_classic())

#making scree plots ----
#scree plots = show you different PCOa axes - can choose different axes sometimes based on data


pseq_fam <- microbiome::aggregate_rare(pseq,level = "Family", detection = 50/100, prevalence = 70/100)

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

pseq.core <- microbiome::transform(pseq.core, "compositional")

set.seed(4235421)

ord <- ordinate(pseq, "MDS", "bray")


eigenvalues <- ord$values$Eigenvalues
eigenvalues <- eigenvalues^2

#select axes for plot

plot_ordination(pseq, ord, type = "samples", axes = 3:4,
                color = "Treatment", shape = "Age", label = NULL, title = NULL,
                justDF = FALSE) + theme_classic() + geom_point(size = 4)



#for other po modification such as plotting taxa information see below
#source: https://joey711.github.io/phyloseq/plot_ordination-examples.html

#Remove OTUs that do not show appear more than 5 times in more than half the samples

wh0 = genefilter_sample(pseq, filterfun_sample(function(x) x > 5), A=0.5*nsamples(pseq))

GP1 = prune_taxa(wh0, pseq)

#transform to even sample depth

GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

#Keep only the most abundant five families.

phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Family"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Family"] %in% top5phyla), GP1)

#38 taxa remaining

#Plot by taxa ----
#Let’s start by plotting just the OTUs, and shading the points by Phylum

GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 = plot_ordination(GP1, GP.ord, type="taxa", color="Family", title="taxa")
print(p1)

#facet by phylum

p1 + facet_wrap(~Family, 3)

#plot just samples ----

p2 = plot_ordination(GP1, GP.ord, type="samples", color="Treatment", shape="Age") 
p2 + geom_polygon(aes(fill=Treatment)) + geom_point(size=4) + ggtitle("samples")

#biplot ----

p3 = plot_ordination(GP1, GP.ord, type="biplot", color="Treatment", shape="Family", title="biplot")
# Some stuff to modify the automatic shape scale
GP1.shape.names = get_taxa_unique(GP1, "Family")
GP1.shape <- 15:(15 + length(GP1.shape.names) - 1)
names(GP1.shape) <- GP1.shape.names
GP1.shape["samples"] <- 16
p3 + scale_shape_manual(values=GP1.shape)

p4 = plot_ordination(GP1, GP.ord, type="split", color="Family", shape="Age", label="Treatment", title="split") 
p4

#to override code and make factor (treatment) black
gg_color_hue <- function(n){
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}
color.names <- levels(p4$data$Phylum)
p4cols <- gg_color_hue(length(color.names))
names(p4cols) <- color.names
p4cols["samples"] <- "black"
p4 + scale_color_manual(values=p4cols)

#use loop to run through different ordination methods - note: need a phy_tree within phyloseq object

dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, pseq, dist){
  ordi = ordinate(pseq, method=i, distance=dist)
  plot_ordination(pseq, ordi, "samples", color="Treatment")
}, GP1, dist)

#Use the ordinate function to simultaneously perform weightd UniFrac and then perform a Principal Coordinate Analysis on that distance matrix (first line). Next pass that data and the ordination results to plot_ordination to create the ggplot2 output graphic with default ggplot2 settings.

ordu = ordinate(GP1, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(GP1, ordu, color="Treatment", shape="Age")



p = p + geom_point(size=7, alpha=0.75)
p = p + scale_colour_brewer(type="qual", palette="Set1")
