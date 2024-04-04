#Abundance plots using psmelt
#Written by Marissa WL (2024-01-18)
#workflow: Load data (pseq) -> filter data (if you want) -> manipulate data (e.g., top families, relative composition, etc.) -> psmelt -> plot


library("devtools")
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggplot2)
library(dplyr)

#Load data ----

Marissa_MU42022_rarefied_20231016 <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rarefied_20231016.rds")

pseq<- Marissa_mb2021_filtered_20240203

pseq <-  Marissa_MU42022_rarefied_20231016

#set theme ----

theme.marissa <- function() {
  theme_classic(base_size = 14) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold"))
}

theme_set(theme.marissa())

#Sample selection ----

#Remove F4 (MU42022)

pseq <- subset_samples(pseq, !Genetics %in% c("4"))

#Remove day 3 (only 1 sample remaining) for mb2021 project

pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))

#Remove algae
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#Spat only
pseq <- subset_samples(pseq, !Age %in% c("Day 01", "Day 03", "Day 06", "Day 15"))

#Day 1 only
pseq <- subset_samples(pseq, Age %in% c("1 dpf"))






#Relative composition ----

#convert to compositional data

pseq.rel <- microbiome::transform(pseq, "compositional")

pseq_fam <- microbiome::aggregate_rare(pseq.rel, level = "Family", detection = 50/100, prevalence = 50/100)
#pseq_fam <- microbiome::aggregate_rare(pseq.rel, level = "Family", detection = 1/100, prevalence = 70/100)

pseq_phy <- microbiome::aggregate_rare(pseq.rel, level = "Phylum", detection = 50/100, prevalence = 50/100)
#pseq_phy <- microbiome::aggregate_rare(pseq.rel, level = "Phylum", detection = 1/100, prevalence = 70/100)

pseq_gen <- microbiome::aggregate_rare(pseq.rel, level = "Genus", detection = 50/100, prevalence = 70/100)
#pseq_phy <- microbiome::aggregate_rare(pseq.rel, level = "Genus", detection = 1/100, prevalence = 70/100)

pseq.core <- core(pseq, detection = .1/100, prevalence = 90/100)


pseq_group <- pseq.rel %>%
  aggregate_taxa(level = "Phylum") %>%
  transform(transform = "compositional") %>%
  psmelt() %>%
  group_by(Age)


pseq_group <- pseq.rel %>%
  aggregate_taxa(level = "Family") %>%
  transform(transform = "clr") %>%
  psmelt() %>%
  group_by(Age)


View(pseq_group)

ggplot(pseq_group, aes(fill=Phylum, y=Abundance, x=Age)) + 
  geom_bar(position="stack", stat="identity")


#Top families relative abundance ----

top5F.names = sort(tapply(taxa_sums(pseq.rel), tax_table(pseq.rel)[, "Family"], sum), TRUE)[1:10]

top5F = subset_taxa(pseq.rel, Family %in% names(top5F.names))

#Top Genus relative abundance ----

top10G.names = sort(tapply(taxa_sums(pseq.gen.rel), tax_table(pseq.gen.rel)[, "Genus"], sum), TRUE)[1:10]

top10G = subset_taxa(pseq.gen.rel, Genus %in% names(top10G.names))

#psmelt ----

#all data

top5F <- psmelt(top5F) 

top10G <- psmelt(top10G)

top10 <- psmelt(top10)

top10C <- psmelt(top10C)

pseq <- psmelt(pseq.rel)

relative_fam <- psmelt(pseq_fam)

relative_phy <- psmelt(pseq.phy.rel)

core <- psmelt(pseq.core)




#Manipulating data ----
#If you plot rel abundances straight away, factors with more samples (e.g., between Age groups)
#will look like false higher abundances. To correct for this, average samples.

#look at sample numbers

sample_counts <- top5F %>%
  group_by(Age, Treatment) %>%
  summarise(Num_Samples = n())

print(sample_counts)


#uneven sample numbers - average relative abundances by groups
#1st create new column that contains Treat*Age information 

sample_data <- pseq.rel@sam_data

# Create a new column combining Age and Treatment ----
sample_data$Age_Treatment <- paste(sample_data$Age, sample_data$Treatment, sep = "_")
sample_data(pseq.rel) <- sample_data
View(pseq.rel@sam_data)

pseq <- psmelt(pseq.rel)
View(pseq)

#Now average relative abundance data by column "Age_Treatment"
Avg_abundance <- top5F %>%
  group_by(Treatment, Age, Family) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )

Avg_abundance <- pseq %>%
  group_by(Age_Treatment, Phylum) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )
View(Avg_abundance)

#composition boxplots ----

ggplot(Avg_abundance, aes(fill=Family, y=Avg_Abundance, x=Treatment)) + 
  geom_bar(position="dodge", stat="identity") + scale_fill_ipsum()

ggplot(Avg_abundance, aes(fill=Phylum, y=Avg_Abundance, x=Age_Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Spectral") +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1))

ggplot(Avg_abundance, aes(fill=Family, y=Avg_Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  #scale_fill_brewer(palette = "Paired") +
  labs(title = "All time-points", x = "", y = "Average relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Age) +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1))
  

ggplot(top10G, aes(fill=Genus, y=Abundance, x=Treatment)) + 
  geom_bar(position="dodge", stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "All time-points", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5))

#if need more colours for colour theme
new_palette <- c(brewer.pal(12, "Paired"), "orange")

ggplot(core, aes(fill=Family, y=Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5) +
          facet_grid(~pseq$Age))

ggplot(core, aes(fill = Family, y = Abundance, x = Treatment)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Age) +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 11)) +
  theme(axis.title.y = element_text(size = 12)) +  # Adjust the size of y-axis title here
  theme(legend.text = element_text(size = 11)) +
  theme(legend.title = element_text(size = 12))

  

ggplot(relative_fam, aes(fill=Family, y=Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "All time-points", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Age) +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) +
  theme(axis.text.y = element_text(size=11)) +
  theme(legend.text = element_text(size = 11)) +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1))


ggplot(top10, aes(fill=Order, y=Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "All time-points", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Age) +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) +
  theme(axis.text.y = element_text(size=11)) +
  theme(legend.text = element_text(size = 11))


ggplot(top5F, aes(fill=Family, y=Abundance, x=Treatment)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(title = "All time-points", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~Age) +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) +
  theme(axis.text.y = element_text(size=11)) +
  theme(legend.text = element_text(size = 11))



ggplot(Avg_abundance, aes(fill=Family, y=Avg_Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_ipsum() +
  labs(title = "Larval samples", x = "", y = "Average abundance") +
  theme(plot.title = element_text(hjust = 0.5)) 


#Vibrio subset ----

nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

#color blind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 = c("#CC6666", "#9999CC", "#66CC99")

pseq.rel <- microbiome::transform(pseq, "compositional")

pseq.rel.vibrio <- subset_taxa(pseq.rel, Family=="Vibrionaceae")

pseq.rel.vibrio <- psmelt(pseq.rel.vibrio)

View(pseq.rel.vibrio)

#Looks like high vibrio abundance comes from 1 tank (F2H1) in high sal treatment
#need to check survival data to see if this family (F2) had higher survival

filtered_data <- pseq.rel.vibrio %>% 
filter(Library_Name != "F2H1")

ggplot(pseq.rel.vibrio, aes(fill=Treatment, y=Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "All time-points - F2H1 removed", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) + 
  facet_wrap(~Age) + 
  #scale_fill_hue(c=45, l=80) +
  scale_fill_brewer(palette = "Dark2") 

ggplot(test) + geom_bar(aes(x=a, y=b, fill=c), colour="black", stat="identity")


ggplot(pseq.rel.vibrio, aes(fill = Genus, y = Abundance, x = Treatment)) + 
  geom_boxplot(position = position_dodge(width = 0.8)) +  # Use geom_boxplot for box plots
  labs(title = "All time-points", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1)) + 
  scale_fill_manual(values = mycolors) +
  facet_wrap(~Age)

  
ggplot(pseq.rel.vibrio, aes(fill = Genus, y = Abundance, x = Treatment)) + 
  geom_boxplot(position="stack") +  # Change from geom_bar to geom_boxplot
  labs(title = "", x = "", y = "Relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 11, angle = 45, hjust = 1)) + 
  facet_wrap(~sample_Family) +
  scale_fill_manual(values=mycolors)





Avg_abundance <- pseq.rel.vibrio %>%
  group_by(Treatment, sample_Family) %>%
  summarise(
    Avg_Abundance = mean(Abundance),
    SD_Abundance = sd(Abundance),
    .groups = 'drop'
  )

ggplot(Avg_abundance, aes(fill=Treatment, y=Avg_Abundance, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Day 1 - Vibrionaceae", x = "", y = "Avg relative abundance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(size=11, angle=45, hjust=1)) +
  scale_fill_manual(values=mycolors)


#
