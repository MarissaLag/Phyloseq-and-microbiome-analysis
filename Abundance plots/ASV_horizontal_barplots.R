#Horizontal bar plot with abundances of signif ASVs
#Here, for mb2021 project
#Only graphing ASv abundance to compare low and high salinity treatments (not control)
#Graphing ASVs that are signif different to control (i.e., NOT different between LS and HS treatments)

install.packages("dplyr")
install.packages("ggplot2")
library(dplyr)
library(ggplot2)



#First lets graph abundance of specific ASVs
#names_diff_treat_unadj_spat
#[1] "ASV17"   "ASV30"   "ASV37"   "ASV49"   "ASV216"  "ASV237"  "ASV240"  "ASV254"  "ASV271"  "ASV314"  "ASV325"  "ASV347" 
#[13] "ASV348"  "ASV383"  "ASV391"  "ASV412"  "ASV446"  "ASV507"  "ASV553"  "ASV571"  "ASV580"  "ASV595"  "ASV656"  "ASV751" 
#[25] "ASV841"  "ASV845"  "ASV900"  "ASV927"  "ASV964"  "ASV1088" "ASV1162" "ASV1164"


#names_diff_treat_unadj_1dpf
#"ASV21", "ASV34", "ASV40", "ASV42", "ASV44", "ASV48", "ASV51", "ASV60", "ASV76", "ASV79", "ASV90", "ASV91", "ASV92", "ASV100", "ASV109", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV130", "ASV131", "ASV150", "ASV163", "ASV176", "ASV181", "ASV187", "ASV189", "ASV193", "ASV195", "ASV216", "ASV218", "ASV219", "ASV237", "ASV240", "ASV244", "ASV257", "ASV262", "ASV309", "ASV313", "ASV315", "ASV319", "ASV326", "ASV333", "ASV337", "ASV343", "ASV349", "ASV351", "ASV356", "ASV375", "ASV379", "ASV385", "ASV387", "ASV394", "ASV397", "ASV401", "ASV416", "ASV417", "ASV428", "ASV431", "ASV447", "ASV450", "ASV468", "ASV484", "ASV487", "ASV499", "ASV501", "ASV510", "ASV523", "ASV524", "ASV538", "ASV556", "ASV559", "ASV563", "ASV575", "ASV579", "ASV595", "ASV608", "ASV620", "ASV641", "ASV649", "ASV652", "ASV703", "ASV738", "ASV756", "ASV758", "ASV789", "ASV798", "ASV818", "ASV819", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV874", "ASV882", "ASV889", "ASV916", "ASV918", "ASV919", "ASV929", "ASV952", "ASV955", "ASV960", "ASV982", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1030", "ASV1035", "ASV1036", "ASV1061", "ASV1070", "ASV1074", "ASV1078", "ASV1082", "ASV1086", "ASV1110", "ASV1111", "ASV1138", "ASV1146", "ASV1153", "ASV1170"

pseq <- mb2021_filtered_NOT_rarefied
pseq <- subset_samples(pseq, Age %in% "1 dpf")
#pseq <- subset_samples(pseq, !Family %in% "9")
pseq <- subset_samples(pseq, !Treatment %in% c("Low salinity"))

pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4


pseq <- microbiome::transform(pseq, "compositional")

ps <- psmelt(pseq)

#Spat ASVs (unadj)
filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV17",  "ASV30",   "ASV37",   "ASV49",   "ASV216",  "ASV237",  "ASV240",  "ASV254",  "ASV271",  "ASV314", "ASV325", "ASV347", 
                    "ASV348", "ASV383", "ASV391", "ASV412", "ASV446", "ASV507", "ASV553", "ASV571", "ASV580", "ASV595", "ASV656", "ASV751", 
                    "ASV841",  "ASV845",  "ASV900",  "ASV927",  "ASV964",  "ASV1088", "ASV1162", "ASV1164"))

filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV49",   "ASV216",  "ASV237",  "ASV240",  "ASV254",  "ASV271",  "ASV314", "ASV325", "ASV347", 
                    "ASV348", "ASV383", "ASV391", "ASV412", "ASV446", "ASV507", "ASV553", "ASV571", "ASV580", "ASV595", "ASV656", "ASV751", 
                    "ASV841",  "ASV845",  "ASV900",  "ASV927",  "ASV964",  "ASV1088", "ASV1162", "ASV1164"))

#Spat unadj LS removed - testing only signif different ASVs between COntrol and HS
filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV49", "ASV240", "ASV325", "ASV348", "ASV412", "ASV580", "ASV656", "ASV751", "ASV841", "ASV927", "ASV1088", "ASV1164"))

filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV17", "ASV37"))

#Spat unadj HS removed 
filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV49", "ASV571", "ASV595", "ASV751", "ASV927", "ASV1088", "ASV1162"))

#Day 1 ASVs (unadj)
filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV21", "ASV34", "ASV40", "ASV42", "ASV44", "ASV48", "ASV51", "ASV60", "ASV76", "ASV79", "ASV90", "ASV91", "ASV92", "ASV100")) 

#1 dpf (signif ASVs) with HS removed

filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV34", "ASV40", "ASV48", "ASV51", "ASV60", "ASV79", "ASV90", "ASV91", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV131", "ASV163", "ASV181", "ASV187", "ASV189", "ASV195", "ASV219", "ASV240", "ASV257", "ASV309", "ASV313", "ASV319", "ASV326", "ASV333", "ASV337", "ASV343", "ASV356", "ASV387", "ASV394", "ASV417", "ASV431", "ASV447", "ASV450", "ASV468", "ASV484", "ASV487", "ASV501", "ASV559", "ASV563", "ASV575", "ASV579", "ASV595", "ASV608", "ASV620", "ASV641", "ASV649", "ASV738", "ASV798", "ASV818", "ASV819", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV874", "ASV882", "ASV916", "ASV929", "ASV955", "ASV960", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1035", "ASV1036", "ASV1061", "ASV1070", "ASV1082", "ASV1086", "ASV1111", "ASV1146", "ASV1153", "ASV1170"))

#1 dpf (signif ASVs) with LS removed

filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV34", "ASV40", "ASV44")) 

filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV21", "ASV48", "ASV51", "ASV76", "ASV79", "ASV90", "ASV100", "ASV114", "ASV118", "ASV119", "ASV124")) 
                    
"ASV21", "ASV48", "ASV51", "ASV76", "ASV79", "ASV90", "ASV100", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV131", "ASV163", "ASV176", "ASV181", "ASV189", "ASV193", "ASV195", "ASV216", "ASV219", "ASV237", "ASV244", "ASV309", "ASV319", "ASV333", "ASV337", "ASV343", "ASV349", "ASV351", "ASV356", "ASV375", "ASV385", "ASV387", "ASV394", "ASV397", "ASV401", "ASV416", "ASV417", "ASV428", "ASV431", "ASV468", "ASV484", "ASV487", "ASV499", "ASV510", "ASV524", "ASV538", "ASV559", "ASV575", "ASV595", "ASV649", "ASV652", "ASV703", "ASV738", "ASV756", "ASV798", "ASV818", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV882", "ASV889", "ASV916", "ASV929", "ASV952", "ASV960", "ASV982", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1030", "ASV1035", "ASV1036", "ASV1061", "ASV1074", "ASV1078", "ASV1082", "ASV1086", "ASV1110", "ASV1111", "ASV1138", "ASV1146"


# Calculate the average abundance for each treatment group
average_abundance <- filtered_ps %>%
  group_by(Treatment, Family, OTU) %>%
  summarise(Average_Abundance = mean(Abundance)) 


# Plot
paired_palette <- brewer.pal(12, "Paired")
paired_palette <- c(paired_palette, "#ff7f00", "pink", "red", "yellow", "lightgreen")
paired_palette <- c(brewer.pal(8, "Dark2"),  # 8 more colors from the Dark2 palette
                       brewer.pal(8, "Set2"),   # 8 more colors from the Set2 palette
                       brewer.pal(8, "Accent"), # 8 more colors from the Accent palette
                       brewer.pal(6, "Set3"))

ggplot(average_abundance, aes(fill = Family, y = Average_Abundance, x = Treatment)) + 
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(values = paired_palette) +
  labs(title = "", 
       x = "", 
       y = "Average Relative Abundance", 
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),  # Adjust y-axis title font size
        legend.text = element_text(size = 12),   # Adjust legend text font size
        legend.title = element_text(size = 14),
        panel.border = element_blank()) +  # Adjust legend title font size
  facet_wrap("OTU")

ggplot(average_abundance, aes(fill = OTU, y = Average_Abundance, x = Treatment)) + 
  geom_bar(position = "stack", stat = "identity", color = "black") +
  labs(title = "", 
       x = "", 
       y = "Average Relative Abundance", 
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),  # Adjust y-axis title font size
        legend.text = element_text(size = 12),   # Adjust legend text font size
        legend.title = element_text(size = 14),
        panel.border = element_blank()) +  # Adjust legend title font size
  facet_wrap("OTU")


#horizontal bar plot

pseq <- mb2021_filtered_NOT_rarefied
pseq <- subset_samples(pseq, Age %in% "Spat")
pseq <- subset_samples(pseq, !Family %in% "9")
pseq <- microbiome::transform(pseq, "compositional")

# Subset the phyloseq object to include only high and low salinity treatments
filtered_phyloseq <- subset_samples(mb2021_filtered_NOT_rarefied, Treatment %in% c("High salinity", "Low salinity"))

#signif ASVs
selected_asvs <- c("ASV17",  "ASV30",   "ASV37",   "ASV49",   "ASV216",  "ASV237",  "ASV240",  "ASV254",  "ASV271",  "ASV314", "ASV325", "ASV347", 
                   "ASV348", "ASV383", "ASV391", "ASV412", "ASV446", "ASV507", "ASV553", "ASV571", "ASV580", "ASV595", "ASV656", "ASV751", 
                   "ASV841",  "ASV845",  "ASV900",  "ASV927",  "ASV964",  "ASV1088", "ASV1162", "ASV1164")

# Subset the phyloseq object to keep only the selected ASVs
filtered_phyloseq <- prune_taxa(selected_asvs, filtered_phyloseq)

# Melt the phyloseq object to a long data frame
melted_phyloseq <- psmelt(filtered_phyloseq)


average_abundance <- melted_phyloseq %>%
  group_by(Treatment, OTU, Age) %>%
  summarise(Average_Abundance = mean(Abundance)) %>%
mutate(Log_Average_Abundance = log(Average_Abundance))  # Adding 1 to avoid log(0)



# Convert ASV IDs to factors to control the order in the plot
melted_phyloseq$OTU <- factor(melted_phyloseq$OTU, levels = selected_asvs)
melted_phyloseq$Abundance <- ifelse(melted_phyloseq$Treatment == "low salinity", -melted_phyloseq$Abundance, melted_phyloseq$Abundance)

ggplot(average_abundance, aes(x = Log_Average_Abundance, y = OTU, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  facet_wrap(~ Treatment) +
  labs(x = "Log Relative Abundance", y = "", title = "Abundance of Signif ASVs - Spat") +
  theme(legend.position = "bottom")

