###top phylum/family/genera within microbiome data
##bar plots and violin plots

#Quick lib ----

library("devtools")
library(phyloseq)
library(microbiome)

#Load data ----


pseq <-  Marissa_MU42022_rare
pseq <- Marissa_MU42022_rarefied_20231016


#Create objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree


#convert to compositional data


pseq.rel <- microbiome::transform(pseq, "compositional")

pseq.fam.rel <- microbiome::transform(pseq_fam, "compositional")

pseq.phy.rel <- microbiome::transform(pseq_phy, "compositional")

pseq.gen.rel <- microbiome::transform(pseq_gen, "compositional")

pseq.core <- microbiome::transform(pseq.core, "compositional")


pseq_fam <- microbiome::aggregate_rare(pseq, level = "Family", detection = 50/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(pseq, level = "Phylum", detection = 50/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(pseq, level = "Genus", detection = 50/100, prevalence = 70/100)

pseq.core <- core(pseq.fam.rel, detection = .1/100, prevalence = 90/100)

#source for below: https://rstudio-pubs-static.s3.amazonaws.com/330760_8bba830836324bf6b100d4e76f49e3d2.html#other_visualizations 


#Bar plots ----

#PHYLUM LEVEL ----

##STARTING AT PHYLUM LEVEL - TOTAL ABUNDANCE

plot_bar(pseq.rel, fill="Phylum", facet_grid = "Age") 

#Sort the Phyla by abundance and pick the top 5

top5P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Phylum"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 10 Phyla

top5P = subset_taxa(pseq, Phylum %in% names(top5P.names))

#Plot

plot_bar(top5P, x="Age", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(top5P, x="Age", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

##code to change x-axis labels

plot <- plot_bar(top5P, x = "Age", fill = "Phylum") +
  geom_boxplot(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack")

plot + scale_x_discrete(labels = c("Day 1", "Day 18", "Day 3", "Spat"))


#Plot

plot_bar(top5P, x="Treatment", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")


plot_bar(top5P, x="Treatment", fill="Phylum", facet_grid = ~Phylum) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae"))

plot_bar(top5P, x="Treatment", fill="Phylum", facet_grid = ~Phylum) + 
  geom_boxplot(aes(color=Phylum, fill=Phylum), 
  stat="identity", 
  position="stack") 
+ theme_bw() 
+ theme(panel.grid = element_blank()) 
+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
+ theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)   
)



#Phylum relative abundance ----

###PHYLUM LEVEL - instead of using pseq (total abundance) plot with pseq.fam.rel (relative abundance at Family level)

top5P.names = sort(tapply(taxa_sums(pseq.phy.rel), tax_table(pseq.phy.rel)[, "Phylum"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 10 Phyla

top5P = subset_taxa(pseq.phy.rel, Phylum %in% names(top5P.names))


  labs(y = "Relative Abundance")

  plot_bar(top5P, x="Treatment", fill="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
  
  
  plot_bar(top5P, x="Treatment", fill="Phylum", facet_grid = ~Phylum) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12)   
  ) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Relative abundance", x = " ")


#FAMILY LEVEL ----

#Sort the Phyla by abundance and pick the top 5

top5P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Family"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 5 Family

top5P = subset_taxa(pseq, Family %in% names(top5P.names))

#Plot

plot_bar(top5P, x="Treatment", fill="Family", facet_grid = ~Age) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 12)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")

plot_bar(top5P, x="Treatment", fill="Phylum", facet_grid = ~Phylum) + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 11)   
)plot_

#FAMILY LEVEL - plot by bar by treatment and plot relative abundance

top5P.names = sort(tapply(taxa_sums(pseq.fam.rel), tax_table(pseq.fam.rel)[, "Family"], sum), TRUE)[1:5]

#Cut down the physeq data to only the top 5 Phyla

top5P = subset_taxa(pseq.fam.rel, Family %in% names(top5P.names))

plot_bar(top5P, x="Treatment", fill="Family", facet_grid = ~Age) + geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 12)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Relative abundance", x = " ")



#GENUS LEVEL----
#- sort by treatment

top15P.names = sort(tapply(taxa_sums(pseq), tax_table(pseq)[, "Genus"], sum), TRUE)[1:15]

#Cut down the physeq data to only the top 15 Genus

top15P = subset_taxa(pseq, Genus %in% names(top15P.names))

#Plot


plot_bar(top15P, x="Treatment", fill="Genus", facet_grid = ~Genus) + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")

#plot specific family

family_to_plot <- "Rhodobacteraceae"
pseq_filtered <- subset_taxa(pseq, Family == family_to_plot)

plot_bar(pseq_filtered, x="Treatment", fill="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ", title = "Rhodobacteraceae") +
  theme(plot.title = element_text(hjust = 0.5))


##PLOT SPECIFIC GENUS 

###Flavobacterium ----

# Filter data for only "Flavobacterium" genus
genus_to_plot <- "Flavobacterium"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

# Create the bar plot for "Flavobacterium" abundance between different "Tank_treatments"
plot <- plot_bar(pseq_filtered, x = "Treatment", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack")

plot_bar(pseq_filtered, x="Treatment", fill="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")


###Celeribacter ----

genus_to_plot <- "Celeribacter"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

# Create the bar plot for "Flavobacterium" abundance between different "Tank_treatments"
plot <- plot_bar(pseq_filtered, x = "Treatment", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack")

plot_bar(pseq_filtered, x="Treatment", fill="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")

###Phaeobacter ----

genus_to_plot <- "Phaeobacter"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

# Create the bar plot for "Flavobacterium" abundance between different "Tank_treatments"
plot <- plot_bar(pseq_filtered, x = "Treatment", fill = "Genus") +
  geom_bar(aes(color = Genus, fill = Genus), stat = "identity", position = "stack")

plot_bar(pseq_filtered, x="Treatment", fill="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")

###Vibrio ----

# Filter data for only "Vibrio" genus
genus_to_plot <- "Vibrio"
pseq_filtered <- subset_taxa(pseq, Genus == genus_to_plot)

plot_bar(pseq_filtered, x="Treatment", fill="Genus") + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme_bw() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
  axis.text.y = element_text(size = 10),
  strip.text = element_text(size = 10)   
) + scale_x_discrete(labels = c("Control", "High temperature", "Probiotics", "Probiotics + HT", "Algae")) + labs(y = "Abundance", x = " ")


#Plot ASVs ----

#In MU42022: ASV3, 88, and 201 signif different in PB/PB+heat from control at spat stage (inval results)
pseq <-  Marissa_MU42022_rarefied_20231016
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% c("Spat"))
pseq.rel <- microbiome::transform(pseq, "compositional")

psmelt <- psmelt(pseq.rel)

psmelt_filtered <- psmelt %>%
  filter(OTU %in% c("ASV3", "ASV88", "ASV201"))

library(ggplot2)

# Assuming psmelt_filtered is your filtered data frame with abundance values


# Define a colorblind-friendly palette
cud_palette <- c("#009E73", "#E69F00", "#8A79B9","#CC79A7", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Create a box plot
ggplot(psmelt_filtered, aes(x = Treatment, y = Abundance, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = cud_palette) +
  labs(x = "ASV", y = "Abundance") +
  ggtitle("Abundance of Selected ASVs") +
  facet_grid(~OTU) +
  theme(axis.text.x = element_blank()) +
  geom_boxplot(outlier.shape = NA)

  
#theme(axis.text.x = element_text(size=11, angle=45, hjust=1, face = "bold")) +



# Define the desired order of ASVs
asv_order <- c("ASV3", "ASV88", "ASV201")

# Convert the OTU column to a factor with specified levels and order
psmelt_filtered$OTU <- factor(psmelt_filtered$OTU, levels = asv_order)

# Create a box plot
ggplot(psmelt_filtered, aes(x = Treatment, y = Abundance, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +  # Remove outliers
  scale_fill_manual(values = cud_palette) +  # Set the fill colors using the CUD palette
  labs(x = "ASV", y = "Abundance") +
  ggtitle("Abundance of Selected ASVs") +
  facet_grid(~OTU) +
  theme(axis.text.x = element_blank())

View(psmelt_filtered)

#rename ASVs to highest matches on NCBI

asv_names <- c("ASV3" = "Phaeobacter", "ASV88" = "Lokatenella", "ASV201" = "Sulfitobacter")

# Add a new column with the ASV names
psmelt_filtered <- psmelt_filtered %>%
  mutate(ASV_Name = case_when(
    OTU == "ASV3" ~ "Phaeobacter",
    OTU == "ASV88" ~ "Lokatenella",
    OTU == "ASV201" ~ "Sulfitobacter",
    TRUE ~ OTU  # Keep the original OTU value if it doesn't match any of the specified ASVs
  ))


asv_order <- c("Phaeobacter", "Lokatenella", "Sulfitobacter")

# Convert the OTU column to a factor with specified levels and order
psmelt_filtered$ASV_Name <- factor(psmelt_filtered$ASV_Name, levels = asv_order)


ggplot(psmelt_filtered, aes(x = Treatment, y = Abundance, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = cud_palette) +
  labs(x = "", y = "Relative Abundance") +
  ggtitle("") +
  facet_grid(~ASV_Name) +
  theme(axis.text.x = element_blank(),
          panel.border = element_rect(color = "black", fill = NA))
