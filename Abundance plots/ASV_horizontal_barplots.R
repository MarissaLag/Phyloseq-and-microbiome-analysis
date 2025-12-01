#Horizontal bar plot with abundances of signif ASVs
#Here, for mb2021 project
#Only graphing ASv abundance to compare low and high salinity treatments (not control)
#Graphing ASVs that are signif different to control (i.e., NOT different between LS and HS treatments)

install.packages("dplyr")
install.packages("ggplot2")
library(dplyr)
library(ggplot2)
library(phyloseq)
library(RColorBrewer)

#set theme
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

#First lets graph abundance of specific ASVs
#names_diff_treat_unadj_spat
#[1] "ASV17"   "ASV30"   "ASV37"   "ASV49"   "ASV216"  "ASV237"  "ASV240"  "ASV254"  "ASV271"  "ASV314"  "ASV325"  "ASV347" 
#[13] "ASV348"  "ASV383"  "ASV391"  "ASV412"  "ASV446"  "ASV507"  "ASV553"  "ASV571"  "ASV580"  "ASV595"  "ASV656"  "ASV751" 
#[25] "ASV841"  "ASV845"  "ASV900"  "ASV927"  "ASV964"  "ASV1088" "ASV1162" "ASV1164"


#names_diff_treat_unadj_1dpf
#"ASV21", "ASV34", "ASV40", "ASV42", "ASV44", "ASV48", "ASV51", "ASV60", "ASV76", "ASV79", "ASV90", "ASV91", "ASV92", "ASV100", "ASV109", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV130", "ASV131", "ASV150", "ASV163", "ASV176", "ASV181", "ASV187", "ASV189", "ASV193", "ASV195", "ASV216", "ASV218", "ASV219", "ASV237", "ASV240", "ASV244", "ASV257", "ASV262", "ASV309", "ASV313", "ASV315", "ASV319", "ASV326", "ASV333", "ASV337", "ASV343", "ASV349", "ASV351", "ASV356", "ASV375", "ASV379", "ASV385", "ASV387", "ASV394", "ASV397", "ASV401", "ASV416", "ASV417", "ASV428", "ASV431", "ASV447", "ASV450", "ASV468", "ASV484", "ASV487", "ASV499", "ASV501", "ASV510", "ASV523", "ASV524", "ASV538", "ASV556", "ASV559", "ASV563", "ASV575", "ASV579", "ASV595", "ASV608", "ASV620", "ASV641", "ASV649", "ASV652", "ASV703", "ASV738", "ASV756", "ASV758", "ASV789", "ASV798", "ASV818", "ASV819", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV874", "ASV882", "ASV889", "ASV916", "ASV918", "ASV919", "ASV929", "ASV952", "ASV955", "ASV960", "ASV982", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1030", "ASV1035", "ASV1036", "ASV1061", "ASV1070", "ASV1074", "ASV1078", "ASV1082", "ASV1086", "ASV1110", "ASV1111", "ASV1138", "ASV1146", "ASV1153", "ASV1170"

# pseq <- mb2021_filtered_NOT_rarefied
# pseq <- subset_samples(pseq, Age %in% "1 dpf")
# #pseq <- subset_samples(pseq, !Family %in% "9")
# pseq <- subset_samples(pseq, !Treatment %in% c("Low salinity"))
# 
# pseq@sam_data$Family[pseq@sam_data$Family %in% c(9, 13)] <- 1
# pseq@sam_data$Family[pseq@sam_data$Family %in% c(10, 14)] <- 2
# pseq@sam_data$Family[pseq@sam_data$Family %in% c(11, 15)] <- 3
# pseq@sam_data$Family[pseq@sam_data$Family %in% c(12, 16)] <- 4

#MU42022 filtering
pseq <- MU42022_filtered_Oct92024
#pseq <- MU42022_filtered_NOT_rarefied
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% "Spat")
pseq <- subset_samples(pseq, !Treatment %in% "High temperature")
pseq <- microbiome::transform(pseq, "compositional")
ps <- psmelt(pseq)

#PB2023
pseq <- PB2023_spat_not_rarefied_normalized
pseq <- PB2023_spat_10X_limited_CSS
pseq <- subset_samples(pseq, !Treatment %in% c("Continuous Probiotics", "James"))
pseq <- microbiome::transform(pseq, "compositional")
ps <- psmelt(pseq)

#Spat ASVs (unadj)
# filtered_ps <- ps %>% 
#   filter(OTU %in% c("ASV17",  "ASV30",   "ASV37",   "ASV49",   "ASV216",  "ASV237",  "ASV240",  "ASV254",  "ASV271",  "ASV314", "ASV325", "ASV347", 
#                     "ASV348", "ASV383", "ASV391", "ASV412", "ASV446", "ASV507", "ASV553", "ASV571", "ASV580", "ASV595", "ASV656", "ASV751", 
#                     "ASV841",  "ASV845",  "ASV900",  "ASV927",  "ASV964",  "ASV1088", "ASV1162", "ASV1164"))
# 
# filtered_ps <- ps %>% 
#   filter(OTU %in% c("ASV49",   "ASV216",  "ASV237",  "ASV240",  "ASV254",  "ASV271",  "ASV314", "ASV325", "ASV347", 
#                     "ASV348", "ASV383", "ASV391", "ASV412", "ASV446", "ASV507", "ASV553", "ASV571", "ASV580", "ASV595", "ASV656", "ASV751", 
#                     "ASV841",  "ASV845",  "ASV900",  "ASV927",  "ASV964",  "ASV1088", "ASV1162", "ASV1164"))
# 
# #Spat unadj LS removed - testing only signif different ASVs between Control and HS
# filtered_ps <- ps %>% 
#   filter(OTU %in% c("ASV49", "ASV240", "ASV325", "ASV348", "ASV412", "ASV580", "ASV656", "ASV751", "ASV841", "ASV927", "ASV1088", "ASV1164"))
# 
# filtered_ps <- ps %>% 
#   filter(OTU %in% c("ASV17", "ASV37"))
# 
# #Spat unadj HS removed 
# filtered_ps <- ps %>% 
#   filter(OTU %in% c("ASV49", "ASV571", "ASV595", "ASV751", "ASV927", "ASV1088", "ASV1162"))
# 
# #Day 1 ASVs (unadj)
# filtered_ps <- ps %>% 
#   filter(OTU %in% c("ASV21", "ASV34", "ASV40", "ASV42", "ASV44", "ASV48", "ASV51", "ASV60", "ASV76", "ASV79", "ASV90", "ASV91", "ASV92", "ASV100")) 
# 
# #1 dpf (signif ASVs) with HS removed
# 
# filtered_ps <- ps %>% 
#   filter(OTU %in% c("ASV34", "ASV40", "ASV48", "ASV51", "ASV60", "ASV79", "ASV90", "ASV91", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV131", "ASV163", "ASV181", "ASV187", "ASV189", "ASV195", "ASV219", "ASV240", "ASV257", "ASV309", "ASV313", "ASV319", "ASV326", "ASV333", "ASV337", "ASV343", "ASV356", "ASV387", "ASV394", "ASV417", "ASV431", "ASV447", "ASV450", "ASV468", "ASV484", "ASV487", "ASV501", "ASV559", "ASV563", "ASV575", "ASV579", "ASV595", "ASV608", "ASV620", "ASV641", "ASV649", "ASV738", "ASV798", "ASV818", "ASV819", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV874", "ASV882", "ASV916", "ASV929", "ASV955", "ASV960", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1035", "ASV1036", "ASV1061", "ASV1070", "ASV1082", "ASV1086", "ASV1111", "ASV1146", "ASV1153", "ASV1170"))
# 
# #1 dpf (signif ASVs) with LS removed
# 
# filtered_ps <- ps %>% 
#   filter(OTU %in% c("ASV34", "ASV40", "ASV44")) 
# 
# filtered_ps <- ps %>% 
#   filter(OTU %in% c("ASV21", "ASV48", "ASV51", "ASV76", "ASV79", "ASV90", "ASV100", "ASV114", "ASV118", "ASV119", "ASV124")) 
#                     
# "ASV21", "ASV48", "ASV51", "ASV76", "ASV79", "ASV90", "ASV100", "ASV114", "ASV118", "ASV119", "ASV124", "ASV126", "ASV131", "ASV163", "ASV176", "ASV181", "ASV189", "ASV193", "ASV195", "ASV216", "ASV219", "ASV237", "ASV244", "ASV309", "ASV319", "ASV333", "ASV337", "ASV343", "ASV349", "ASV351", "ASV356", "ASV375", "ASV385", "ASV387", "ASV394", "ASV397", "ASV401", "ASV416", "ASV417", "ASV428", "ASV431", "ASV468", "ASV484", "ASV487", "ASV499", "ASV510", "ASV524", "ASV538", "ASV559", "ASV575", "ASV595", "ASV649", "ASV652", "ASV703", "ASV738", "ASV756", "ASV798", "ASV818", "ASV820", "ASV837", "ASV838", "ASV845", "ASV846", "ASV872", "ASV882", "ASV889", "ASV916", "ASV929", "ASV952", "ASV960", "ASV982", "ASV983", "ASV993", "ASV1005", "ASV1007", "ASV1030", "ASV1035", "ASV1036", "ASV1061", "ASV1074", "ASV1078", "ASV1082", "ASV1086", "ASV1110", "ASV1111", "ASV1138", "ASV1146"


#MU42022 spat inval results
filtered_ps <- ps %>% 
   filter(OTU %in% c("ASV88", "ASV178")) 

#PB2023 results
filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV275", "ASV231", "ASV190", "ASV227")) 

filtered_ps <- ps %>% 
  filter(OTU %in% c("ASV1211", "ASV227", "ASV231", 
                    "ASV461", "ASV468", "ASV190", "ASV777", "ASV275", "ASV236", "ASV348"))

filtered_ps$Genus[filtered_ps$OTU == "ASV227" & is.na(filtered_ps$Genus)] <- "Litoreibacter"

filtered_ps$Genus[filtered_ps$OTU == "ASV461" & is.na(filtered_ps$Genus)] <- "Unknown Saccharospirillaceae"

filtered_ps$Genus[filtered_ps$OTU == "ASV236" & is.na(filtered_ps$Genus)] <- "Unknown Methyloligellaceae"

filtered_ps$Genus[filtered_ps$OTU == "ASV348" & is.na(filtered_ps$Genus)] <- "Unknown Microtrichaceae"


# Calculate the average abundance for each treatment group
average_abundance <- filtered_ps %>%
  group_by(Treatment, OTU, Genus) %>%
  summarise(Average_Abundance = mean(Abundance),
            std_Abundance = sd(Abundance))


# Plot

paired_palette <- brewer.pal(12, "Paired")
paired_palette <- c(paired_palette, "#ff7f00", "pink", "red", "yellow", "lightgreen")
paired_palette <- c(brewer.pal(8, "Dark2"),  # 8 more colors from the Dark2 palette
                       brewer.pal(8, "Set2"),   # 8 more colors from the Set2 palette
                       brewer.pal(8, "Accent"), # 8 more colors from the Accent palette
                       brewer.pal(6, "Set3"))


#With error bars or withotu error bars
ggplot(average_abundance, aes(y = Average_Abundance, x = Treatment, fill = Genus)) + 
  geom_bar(position = "stack", stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin = Average_Abundance - std_Abundance, 
                    ymax = Average_Abundance + std_Abundance),
                width = 0.2, 
                position = position_dodge(width = 0.9)) +  # Adjusts error bar positioning
  #scale_fill_manual(values = paired_palette) +
  scale_fill_viridis_d(option = "rocket") + 
  labs(title = "", 
       x = "", 
       y = "Average Relative Abundance", 
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.border = element_blank()) + 
  facet_wrap(~OTU)


#As box plot 
#change order
# Set OTU as a factor with the specified order
filtered_ps$OTU <- factor(filtered_ps$OTU, levels = c("ASV275", "ASV777", "ASV348", "ASV236", "ASV190", 
                                                      "ASV1211", "ASV227", "ASV231", "ASV461", "ASV478"))

# Create the OTU_Genus column while maintaining the OTU order
filtered_ps <- filtered_ps %>%
  mutate(OTU_Genus = factor(paste(OTU, Genus, sep = " ; "), levels = unique(paste(OTU, Genus, sep = " ; "))))

filtered_ps <- filtered_ps %>%
  mutate(OTU_Genus = paste(OTU, Genus, sep = " ; ")) %>%
  mutate(OTU_Genus = factor(OTU_Genus, levels = paste(levels(filtered_ps$OTU), 
                                                      Genus[match(levels(filtered_ps$OTU), OTU)], 
                                                      sep = " ; ")))

filtered_ps$Genus <- factor(filtered_ps$Genus, 
                            levels = unique(filtered_ps$Genus[order(filtered_ps$OTU_Genus)]))

# Generate the plot
ggplot(filtered_ps, aes(x = Treatment, y = Abundance, fill = Family)) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Creates the boxplot without showing outliers
  geom_jitter(aes(color = Family), position = position_jitter(width = 0.2, height = 0), 
              size = 1.5, alpha = 0.8) +  # Adds jittered points for individual data
  scale_fill_viridis_d(option = "rocket") + 
  scale_color_viridis_d(option = "rocket", guide = "none") +  # Keeps jitter points colored but removes extra legend
  labs(title = "", 
       x = "", 
       y = "Average Relative Abundance", 
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        panel.border = element_blank()) + 
  facet_wrap(~OTU_Genus, scales = "free_y") +  
  scale_x_discrete(labels = c("Control" = "Control", 
                              "Killed-Probiotics" = "Killed-Bacteria Added",
                              "Probiotics" = "Bacteria Added"))


ggplot(filtered_ps, aes(x = Treatment, y = Abundance, fill = Family)) + 
  geom_violin() +  # Violin plot with scaled width and no trimming
  #geom_jitter(aes(color = Family), position = position_jitter(width = 0.2, height = 0), 
             # size = 1.5, alpha = 0.8) +  # Adds jittered points for individual data
  scale_fill_viridis_d(option = "rocket") + 
  scale_color_viridis_d(option = "rocket", guide = "none") +  # Keeps jitter points colored but removes extra legend
  labs(title = "", 
       x = "", 
       y = "Average Relative Abundance", 
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.border = element_blank()) + 
  facet_wrap(~Genus) +
  scale_x_discrete(labels = c("Control" = "Control", 
                              "Killed-Probiotics" = "Killed-Bacteria Added",
                              "Probiotics" = "Bacteria Added"))




# Update Treatment labels
filtered_ps <- filtered_ps %>%
  mutate(Treatment = recode(Treatment,
                            "Killed-Probiotics" = "Killed-Bacteria Added",
                            "Probiotics" = "Bacteria Added"))
         

filtered_ps$Treatment <- factor(filtered_ps$Treatment, 
                                levels = c("Control", "Killed-Bacteria Added", "Bacteria Added"))


ggplot(filtered_ps, aes(x = Treatment, y = Abundance, fill = Treatment)) + 
  geom_boxplot() + # Create the boxplot
   scale_fill_manual(values = c("Control" = "darkgrey", 
                                 "Killed-Bacteria Added" = "orange", 
                             "Bacteria Added" = "#00A572")) +
  labs(title = "", 
       x = "", 
       y = "Relative Abundance", 
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        panel.border = element_blank()) + 
  facet_wrap("OTU",  labeller = custom_labeller) # Create separate plots for each OTU

# Custom labeller function
# custom_labeller <- function(variable, value) {
#   if (variable == "OTU") {
#     return(c("ASV1211" = "ASV1211 - Lewinella", "ASV227" = "ASV227 - Litoreibacter", 
#              "ASV231" = "ASV231 - Stappia", "ASV461" = "ASV461 - Unknown", 
#              "ASV478" = "ASV478 - Unknown")[value])
#   }
#   return(value)
# }

custom_labeller <- function(variable, value) {
  if (variable == "OTU") {
    return(c("ASV190" = "ASV190 - Phaeobacter", "ASV227" = "ASV227 = Litoreibacter",
           "ASV231" = "ASV231 - Stappia", "ASV275" = "ASV275 - Loktanella")
           [value])
  }
  return(value)
}

# Plot
ggplot(average_abundance, aes(fill = Family, y = Average_Abundance, x = Treatment)) + 
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(values = paired_palette) +
  labs(title = "", 
       x = "", 
       y = "Average Relative Abundance", 
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),  
        legend.text = element_text(size = 12),   
        legend.title = element_text(size = 14),
        panel.border = element_blank()) +
  facet_wrap(~OTU, labeller = custom_labeller)

#horizontal bar plot ----

pseq <- mb2021_filtered_NOT_rarefied
pseq <- subset_samples(pseq, Age %in% "Spat")
pseq <- subset_samples(pseq, !Family %in% "9")
pseq <- microbiome::transform(pseq, "compositional")


pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, Age %in% "90 dpf")
pseq <- subset_samples(pseq, !Genetics %in% "4")
pseq <- subset_samples(pseq, Treatment %in% c("Control", "Probiotics", "Probiotics + HT"))
pseq <- microbiome::transform(pseq, "compositional")

# Subset the phyloseq object to include only high and low salinity treatments
# filtered_phyloseq <- subset_samples(mb2021_filtered_NOT_rarefied, Treatment %in% c("High salinity", "Low salinity"))
# 
# #signif ASVs
# selected_asvs <- c("ASV17",  "ASV30",   "ASV37",   "ASV49",   "ASV216",  "ASV237",  "ASV240",  "ASV254",  "ASV271",  "ASV314", "ASV325", "ASV347", 
#                    "ASV348", "ASV383", "ASV391", "ASV412", "ASV446", "ASV507", "ASV553", "ASV571", "ASV580", "ASV595", "ASV656", "ASV751", 
#                    "ASV841",  "ASV845",  "ASV900",  "ASV927",  "ASV964",  "ASV1088", "ASV1162", "ASV1164")

selected_asvs <- c("ASV88", "ASV178")

selected_asvs <- c("ASV11", "ASV198", "ASV201", "ASV471", "ASV613")

# Subset the phyloseq object to keep only the selected ASVs
filtered_phyloseq <- prune_taxa(selected_asvs, pseq)

# Melt the phyloseq object to a long data frame
melted_phyloseq <- psmelt(filtered_phyloseq)

#Calulcate log-fold change
average_abundance <- melted_phyloseq %>%
  group_by(Treatment, OTU, Family) %>%
  summarise(Average_Abundance = mean(Abundance), .groups = 'drop') %>%
  # Pivot to create a column for control values
  pivot_wider(names_from = Treatment, values_from = Average_Abundance, names_prefix = "Treatment_") %>%
  # Calculate the log average abundance for each treatment
  mutate(across(starts_with("Treatment_"), log, .names = "Log_{.col}")) %>%
  # Calculate the log-fold change for each treatment against the control
  mutate(across(starts_with("Log_Treatment_"), 
                ~ . - Log_Treatment_Control, 
                .names = "Log_Fold_Change_{.col}")) %>%
  # Select relevant columns for log-fold changes
  select(OTU, starts_with("Log_Fold_Change_"))

average_abundance_long <- average_abundance %>%
  pivot_longer(
    cols = starts_with("Log_Fold_Change_"),
    names_to = "Treatment",
    values_to = "Log_Fold_Change"
  ) %>%
  # Clean up Treatment names if needed
  mutate(Treatment = gsub("Log_Fold_Change_Log_Treatment_", "", Treatment))

average_abundance_long <- average_abundance_long %>%
  filter(Treatment != "Control")

#Update names to include genus
# Update OTU names
average_abundance_long$OTU <- ifelse(
  average_abundance_long$OTU == "ASV88", "ASV88 - Loktanella",
  ifelse(average_abundance_long$OTU == "ASV178", "ASV178 - Loktanella", average_abundance_long$OTU)
)


# Plotting
ggplot(average_abundance_long, aes(x = Log_Fold_Change, y = OTU, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") + # Adds black borders for clarity
  labs(
    title = "",
    x = "Log-Fold Change",
    y = ""
  ) +
  scale_fill_manual(values = c("Probiotics" = "skyblue2", "Probiotics + HT" = "#fc8d62")) + # Custom colors
  theme_bw(base_size = 14) + 
  theme(
    panel.grid = element_blank(),      
    legend.position = "bottom",        
    legend.title = element_blank(),       
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    axis.title.y = element_text(face = "bold"),   
    axis.title.x = element_text(face = "bold"),           
    axis.text = element_text(size = 14, face = "bold"),                 
    legend.text = element_text(size = 14)                
  )

