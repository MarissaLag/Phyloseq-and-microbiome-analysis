#Packages ----
install.packages("indicspecies")
install.packages("randomcoloR")

#load
library(indicspecies)
library(ggplot2)
library(dplyr)
library(microbiome)
library(phyloseq)
library(RColorBrewer)
library(randomcoloR)

#Set data ----

Marissa_MU42022_rare <- readRDS("~/GitHub/mb2021_phyloseq/Marissa_MU42022_rare.rds")

pseq <- MU42022_filtered_Oct92024

pseq <- MU42022_filtered_NOT_rarefied #579 taxa

#pseq <- Marissa_mb2021_filtered_20240203

pseq <- mb2021_filtered_NOT_rarefied_normalized #1007 taxa #use Deseq normalized data???

#Extract abundance matrix
#from the phyloseq object using phyloseq

OTU = pseq@otu_table

OTU1 = as(OTU, "matrix")
write.csv(OTU1, file="Data_fram_1.cvs",row.names=TRUE)

write.table(OTU1,file="data_table_MU42022_oct2024.csv",sep=",",dec = " ")

####Format to example data and reload below for actual test (add metadata)

#reload edited table
data_table <- read.csv("data_table_MU42022_oct2024.csv")

pc_FUN = read.csv("data_table_mb2021_unrarefied.csv", header= TRUE)

pc_FUN <- data_table

#if removing samples ----

#for mb2021 project remove remaining day 3 samples and T9 spat data
#pc_FUN <- pc_FUN[!pc_FUN$`Time-point` == "3 dpf", ]
pc_FUN <- pc_FUN[!pc_FUN$`Tank` == "9", ] 

#MU42022 filtering
pc_FUN <- pc_FUN[!pc_FUN$Genetics == "4", ] 
pc_FUN <- na.omit(pc_FUN)

#Day 1 only 
pc_FUN <- pc_FUN[pc_FUN$'Age' == "Day 15", ]
#Spat only
pc_FUN <- pc_FUN[pc_FUN$`Age` == "Spat", ]

####Test ASVs ----

#Inverse data
funi_df<- t(pc_FUN)

###make into a matrix and populate::: This tells r what is metadata and what is the actual data ... Below 5-952 are the coloumns that are the data
#Note: sum columns add up to zero so you may get an error
#matrix_F = pc_FUN[ ,6:1012] 

dim(pc_FUN)

matrix_F = pc_FUN[ ,8:586] 

View(matrix_F)

#mb2021
#matrix_F = pc_FUN[ ,6:585]

### Make the equation. Saying we want to examine specific column of metadata
#Note: has difficulty testing with more than one factor at a time
time_a_F = pc_FUN$Treatment

### Run test 
inv_F_day1 = multipatt(matrix_F, time_a_F, func = "r.g", control = how(nperm=9999))
results <- summary(inv_F_day1)

inv_F_spat = multipatt(matrix_F, time_a_F, func = "r.g", control = how(nperm=9999))
results <- summary(inv_F_spat)

#No good way to extract signif ASVs, just do it manually
# Subset the matrix for the selected ASVs from `inv_F_day1`
selected_asvs_day1 <- c("ASV531", "ASV27", "ASV17", "ASV323", "ASV172", "ASV397", "ASV85", "ASV316", "ASV87", "ASV135", 
                   "ASV30", "ASV235", "ASV105", "ASV7", "ASV18", "ASV345")

selected_asvs_spat <- c("ASV444", "ASV471", "ASV69", "ASV613", "ASV262", 
                        "ASV11", "ASV198", "ASV201", "ASV241", "ASV174", 
                        "ASV254", "ASV68", "ASV49", "ASV233", "ASV109", 
                        "ASV395", "ASV747", "ASV221", "ASV153", "ASV360", 
                        "ASV88", "ASV178")


significant_ASVs_spat <- data.frame(
  ASV = c("ASV444", "ASV471", "ASV69", "ASV613", "ASV262", 
          "ASV11", "ASV198", "ASV201", "ASV241", "ASV174", 
          "ASV254", "ASV68", "ASV49", "ASV233", "ASV109", 
          "ASV395", "ASV747", "ASV221", "ASV153", "ASV360", 
          "ASV88", "ASV178"),
  Treatment_significant = c("Control", "Probiotics", "Probiotics", 
                            "Probiotics", "Probiotics", "Probiotics", 
                            "Probiotics", "Probiotics", "Probiotics + HT", 
                            "Probiotics + HT", "Probiotics + HT", 
                            "Probiotics + HT", "Probiotics + HT", 
                            "Probiotics + HT", "Probiotics + HT", 
                            "Probiotics + HT", "Probiotics + HT", 
                            "Control + Probiotics", "Control + Probiotics", "Control + Probiotics",
                            "Probiotics + Probiotics + HT", "Probiotics + Probiotics + HT")
                        )

pseq <- psmelt(pseq)

# pseq_ASVs <- pseq %>%
#   filter(Age == "Spat", OTU %in% selected_ASVs_spat) %>%
#   group_by(Treatment,OTU,Age) %>%
#   summarise(avg_abundance = mean(Abundance),
#             std_abundance = sd(Abundance))
# 
# pseq_ASVs <- pseq_ASVs %>%
#   left_join(significant_ASVs_spat, by = c("OTU" = "ASV"))
# 
# pseq_ASVs_control <- pseq_ASVs %>%
#   filter(Treatment_significant == "Control")
# 
# pseq_ASVs_PB <- pseq_ASVs %>%
#   filter(Treatment_significant == "Probiotics")
# 
# pseq_ASVs_PBH <- pseq_ASVs %>%
#   filter(Treatment_significant == "Probiotics + HT")
# 
# pseq_ASVs_PB_PBH <- pseq_ASVs %>%
#   filter(Treatment_significant == "Probiotics + Probiotics + HT")
# 
# 
# ggplot(pseq_ASVs_PB_PBH, aes(x = OTU, y = avg_abundance, fill = Treatment)) +
#   geom_bar(stat = "identity", position = position_dodge(), ) +
#   scale_colour_manual(values = c("darkgrey", "cornflowerblue", "#3CB371")) +
#   scale_fill_manual(values = c("darkgrey", "cornflowerblue", "#3CB371")) +
#   geom_errorbar(aes(ymin = avg_abundance - std_abundance, 
#                     ymax = avg_abundance + std_abundance),
#                 color = "black", 
#                 position = position_dodge(0.9), width = 0.2) +  # Error bars
#   labs(title = "Spat",
#        x = "Operational Taxonomic Unit (OTU)",
#        y = "Average Abundance") +
#   theme_bw() +  
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())  



#Column ("condition") stating which treatment is signif ----

#Make column with ASV names
inv_F_day1$sign$Feature_ID <- row.names(inv_F_day1$sign)
View(inv_F_day1$sign)

inv_F_spat$sign$Feature_ID <- row.names(inv_F_spat$sign)
View(inv_F_spat$sign)

# Extract the relevant columns from the 'sign' data frame
#Assigns a 1 if associated with a group or zero if not

inv_F_sign_df_day1 <- data.frame(
  Feature_ID = inv_F_day1$sign$Feature_ID,
  p.value = inv_F_day1$sign$p.value,
  s.Control = inv_F_day1$sign$s.Control,
  s.High.temperature = inv_F_day1$sign$s.High.temperature,
  s.Probiotics = inv_F_day1$sign$s.Probiotics,
  s.Probiotics.+.HT = inv_F_day1$sign$s.Probiotics.+.HT,
  s.Probiotics.+.HT = inv_F_day1$sign$s.Probiotics.+.HT
)

inv_F_sign_df_spat <- data.frame(
  Feature_ID = inv_F_spat$sign$Feature_ID,
  p.value = inv_F_spat$sign$p.value,
  s.Control = inv_F_spat$sign$s.Control,
  s.High = inv_F_spat$sign$s.High,
  s.Low = inv_F_spat$sign$s.Low
)


# inv_F_sign_df_day1 <- data.frame(
#   Feature_ID = inv_F_day1$sign$Feature_ID,
#   p.value = inv_F_day1$sign$p.value,
#   s.Control = inv_F_day1$sign$s.Control,
#   s.High = inv_F_day1$sign$s.High,
#   s.Low = inv_F_day1$sign$s.Low
# )
# inv_F_sign_df_spat <- data.frame(
#   Feature_ID = inv_F_spat$sign$Feature_ID,
#   p.value = inv_F_spat$sign$p.value,
#   s.Control = inv_F_spat$sign$s.Control,
#   s.High = inv_F_spat$sign$s.High,
#   s.Low = inv_F_spat$sign$s.Low
# )


#below create a "condition" column that will state which treatment is signficant
# Initialize the "condition" column as an empty character vector
inv_F_sign_df_day1$condition <- ""
inv_F_sign_df_spat$condition <- ""

# Check each row for conditions and concatenate letters accordingly
inv_F_sign_df_day1$condition <- ifelse(inv_F_sign_df_day1$s.Control == 1, paste0(inv_F_sign_df_day1$condition, "C"), inv_F_sign_df_day1$condition)
inv_F_sign_df_day1$condition <- ifelse(inv_F_sign_df_day1$s.High == 1, paste0(inv_F_sign_df_day1$condition, "H"), inv_F_sign_df_day1$condition)
inv_F_sign_df_day1$condition <- ifelse(inv_F_sign_df_day1$s.Low == 1, paste0(inv_F_sign_df_day1$condition, "L"), inv_F_sign_df_day1$condition)

inv_F_sign_df_spat$condition <- ifelse(inv_F_sign_df_spat$s.Control == 1, paste0(inv_F_sign_df_spat$condition, "C"), inv_F_sign_df_spat$condition)
inv_F_sign_df_spat$condition <- ifelse(inv_F_sign_df_spat$s.High == 1, paste0(inv_F_sign_df_spat$condition, "H"), inv_F_sign_df_spat$condition)
inv_F_sign_df_spat$condition <- ifelse(inv_F_sign_df_spat$s.Low == 1, paste0(inv_F_sign_df_spat$condition, "L"), inv_F_sign_df_spat$condition)

View(inv_F_sign_df_day1)

#Make a list of signif ASVs ----

#Day 1 ASVs
# Filter rows where p.value is less than 0.05
significant_rows <- inv_F_day1$sign[inv_F_day1$sign$p.value < 0.05, ]

# Extract the Feature_IDs from those rows
significant_ASVs_day1 <- significant_rows$Feature_ID

significant_ASVs_day1 <- significant_ASVs_day1[!is.na(significant_ASVs_day1)]


#Spat ASVs
significant_rows <- inv_F_spat$sign[inv_F_spat$sign$p.value < 0.05, ]

significant_ASVs_spat <- significant_rows$Feature_ID

significant_ASVs_spat <- significant_ASVs_spat[!is.na(significant_ASVs_spat)]


# Find common ASVs between the time points
common_ASVs <- intersect(significant_ASVs_day1, significant_ASVs_spat)




#Load data
pseq<- mb2021_filtered_NOT_rarefied_normalized
pseq.rel <- microbiome::transform(pseq, "compositional")

pseq_day1 <- subset_samples(pseq.rel, Age %in% c("1 dpf"))
ps_day1 <- psmelt2(pseq_day1) #long format

pseq_spat <- subset_samples(pseq.rel, Age %in% c("Spat"))
ps_spat <- psmelt2(pseq_spat) #long format
View(ps_spat)



#Merge objects ---- 
#want to get column "condition" from inv_F object and add to correct ros of ps object
#And filter for signficant ASVs

ps_updated_day1 <- ps_day1 %>%
  left_join(inv_F_sign_df_day1, by = c("FeatureID" = "Feature_ID")) %>%
  filter(FeatureID %in% significant_ASVs_day1)

ps_updated_spat <- ps_spat %>%
  left_join(inv_F_sign_df_spat, by = c("FeatureID" = "Feature_ID")) %>%
  filter(FeatureID %in% significant_ASVs_spat)

#Try plotting ASVs that are only signif more abundant in HS and LS compared to control to make plot cleaner

ps_updated_day1_conditions <- ps_updated_day1 %>%
  filter(condition %in% c("H", "L"))


#Colours ----
nb.cols <- 25
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

# Define your custom colors
my_colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", 
               "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6",
               "#ff33bb", "#ffcc00", "#00ffcc", "#0033ff", "#66ff33")

n <- 50
palette <- distinctColorPalette(n)

custom_shapes <- c(16, 15, 18, 17, 19, 20)  # Choose from a list of available shapes (0-25) in ggplot2

#filter to only include a treatment and caluclate average abundance of signif ASVs
ps_filtered_treatment <- ps_updated_day1_conditions %>%
  filter(Treatment %in% c("High salinity", "Low salinity")) %>%
  group_by(FeatureID, Treatment, condition) %>%
  summarize(mean_abundance = mean(value, na.rm = TRUE))

# Extract unique FeatureID values and set them as factor levels
ps_filtered_treatment$FeatureID <- factor(ps_filtered_treatment$FeatureID, 
                                          levels = unique(ps_filtered_treatment$FeatureID))



# Create the horizontal bar plot with vertical panels

plot <- ggplot(ps_filtered_treatment, aes(x = FeatureID)) +
  geom_col(data = subset(ps_filtered_treatment, Treatment == "High salinity"), 
           aes(y = mean_abundance, fill = 'High salinity'), position = "dodge2") +
  geom_col(data = subset(ps_filtered_treatment, Treatment == "Low salinity"), 
           aes(y = mean_abundance, fill = 'Low salinity'), position = "dodge2") +
  coord_flip() +
  scale_y_continuous(breaks = seq(-20,20, by = 4),
                     labels = (c(seq(20, 0, by = -4), seq(4,20,by=4))))

plot


#For diverging bar plot
#Caculate difference between a treatment and the control

#Select ASVs that are different in High sal treatment only
ps_updated_day1_LC <- ps_updated_day1 %>%
  filter(condition %in% c("H", "C"))

#Calculate difference between High sal and control
ps_filtered_treatment <- ps_updated_day1_LC %>%
  filter(Treatment %in% c("Low salinity", "Control")) %>%
  group_by(FeatureID, Order) %>%
  summarize(mean_abundance_low_salinity = mean(value[Treatment == "Low salinity"], na.rm = TRUE),
            mean_abundance_control = mean(value[Treatment == "Control"], na.rm = TRUE)) %>%
  mutate(diff_abundance = ifelse(mean_abundance_low_salinity > mean_abundance_control, 
                                 mean_abundance_low_salinity - mean_abundance_control, 
                                 mean_abundance_low_salinity - mean_abundance_control)) %>%
  arrange(desc(diff_abundance))

ps_filtered_treatment <- ps_updated_day1_HC %>%
  filter(Treatment %in% c("High salinity", "Control")) %>%
  group_by(FeatureID, Phylum) %>%
  summarize(mean_abundance_high_salinity = mean(value[Treatment == "High salinity"], na.rm = TRUE),
            mean_abundance_control = mean(value[Treatment == "Control"], na.rm = TRUE)) %>%
         mutate(diff_abundance = mean_abundance_high_salinity - mean_abundance_control) %>%
  arrange(desc(diff_abundance))


ps_filtered_treatment <- ps_updated_day1_HC %>%
  filter(Treatment %in% c("High salinity", "Control")) %>%
  group_by(FeatureID, Phylum) %>%
  summarize(mean_abundance_high_salinity = mean(value[Treatment == "High salinity"], na.rm = TRUE),
            mean_abundance_control = mean(value[Treatment == "Control"], na.rm = TRUE)) %>%
  mutate(log_abundance_high_salinity = log(mean_abundance_high_salinity + 1),  # Adding 1 to avoid log(0)
         log_abundance_control = log(mean_abundance_control + 1),
         diff_abundance = log_abundance_high_salinity - log_abundance_control) %>%
  arrange(desc(diff_abundance))

ps_filtered_treatment <- ps_updated_day1_HC %>%
  filter(Treatment %in% c("High salinity", "Control")) %>%
  group_by(FeatureID, Class) %>%
  summarize(mean_abundance_high_salinity = mean(value[Treatment == "High salinity"], na.rm = TRUE),
            mean_abundance_control = mean(value[Treatment == "Control"], na.rm = TRUE)) %>%
  mutate(diff_abundance = mean_abundance_high_salinity - mean_abundance_control,
         sqrt_diff_abundance = sign(diff_abundance) * sqrt(abs(diff_abundance))) %>%  # Square root of the absolute difference with sign
  arrange(desc(sqrt_diff_abundance))

ps_filtered_treatment <- ps_updated_day1_HC %>%
  filter(Treatment %in% c("Low salinity", "Control")) %>%
  group_by(FeatureID, Class) %>%
  summarize(mean_abundance_low_salinity = mean(value[Treatment == "Low salinity"], na.rm = TRUE),
            mean_abundance_control = mean(value[Treatment == "Control"], na.rm = TRUE)) %>%
  mutate(diff_abundance = mean_abundance_low_salinity - mean_abundance_control,
         sqrt_diff_abundance = sign(diff_abundance) * sqrt(abs(diff_abundance))) %>%  # Square root of the absolute difference with sign
  arrange(desc(sqrt_diff_abundance))

my_colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", 
              "#b2df8a", "#fb9a99", "#fdbf6f", "brown", "navy",
               "#ff33bb", "#ffcc00", "#00ffcc", "#0033ff", "#66ff33")

my_colors <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", 
               "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6",
               "#ff33bb", "#ffcc00", "#00ffcc", "#0033ff", "#66ff33",
               "#ff99cc", "#66ccff", "#cc66ff", "#99ff66", "#ff6666",
               "#6666ff", "#ffff66", "#66ff66", "#ff66ff", "#6666ff",
               "#ff6666", "#66ff66", "#ffff66", "#66ff66", "#ff66ff",
               "#6666ff", "#ff6666", "#66ff66", "grey", "lightblue", "darkgreen", "navy", "magenta",
               "purple", "brown")


# Convert FeatureID to a factor with custom order
ps_filtered_treatment$FeatureID <- factor(ps_filtered_treatment$FeatureID, levels = ps_filtered_treatment$FeatureID)
# Plot with customized order
ggplot(ps_filtered_treatment, aes(x = FeatureID, y = sqrt_diff_abundance, fill = Class)) +
  geom_col() +
  scale_fill_manual(values = my_colors) +  # Set custom colors for the fill
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Low Salinity vs Control - Day 1",
       x = "",
       y = "Change in Relative Abundance") +
  theme_bw() +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size =9),
        axis.text.y = element_blank())



#Select ASVs ----

#unrarefied
pseq <- Marissa_MU42022_rare_nochloro
pseq <- subset_samples(pseq, !Genetics %in% "4")
pseq <- subset_samples(pseq, !Organism %in% "Algae")
pseq <- subset_samples(pseq, Age %in% "Spat")

pseq <- microbiome::transform(pseq, "compositional")

pseq <-psmelt(pseq)

#Inval ASVs different at spat stage MU42022, rarefied data
#PB
ASV_list_PB <- pseq %>%
  filter(OTU %in% c("ASV69", "ASV128")) 

ASV_list_PBH <- filter(pseq, OTU %in% c("ASV68", "ASV109", "ASV116", "ASV174", "ASV208", "ASV241", "ASV338"))

ASV_list_PB_PBH <- filter(pseq, OTU %in% c("ASV3", "ASV88", "ASV201")
                      
#Calculate means
average_abundance_PB <- ASV_list_PB %>%
  group_by(OTU, Treatment, Family, Class, Genus) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE))

#if want to group PB and PNH together
average_abundance_PB_PBH <- ASV_list_PB_PBH %>%
  # Combine treatments "Probiotics" and "Probiotics + HT" into a single category
  mutate(Treatment_Group = ifelse(Treatment %in% c("Probiotics", "Probiotics + HT"), "Probiotics_Group", Treatment)) %>%
  # Group by the new treatment category along with OTU, Family, Class, and Genus
  group_by(OTU, Treatment_Group, Family, Class, Genus) %>%
  # Calculate the mean abundance
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE))


# Pivot the data to wide format
average_abundance_wide_PB <- average_abundance_PB %>%
  pivot_wider(names_from = Treatment, values_from = mean_abundance)

# Calculate log2-fold changes
epsilon <- 1e-9  # Small constant to avoid division by zero
average_abundance_wide_PB <- average_abundance_wide_PB %>%
  mutate(
    log2FC_Probiotics_vs_Control = log2((`Probiotics` + epsilon) / (`Control` + epsilon))
    #,log2FC_Probiotics_HT_vs_Control = log2((`Probiotics + HT` + epsilon) / (`Control` + epsilon))
  )

#plot 

custom_palette <- c(  "#8DD3C7","#B3DE69","#FB8072","#FDB462","#CCEBC5", "#FFFFB3", "#BEBADA",    
                      "#80B1D3", "#B3DE69", "#8DD3C7", "#FCCDE5", 
                      "#D9D9D9", "#BC80BD", "#FFED6F")


library(ggplot2)

ggplot(average_abundance_wide_PB, aes(x = OTU, y = log2FC_Probiotics_vs_Control, fill = Family)) +
  geom_col(color = "black") +  # Add black border around bars
  scale_fill_manual(values = custom_palette) +  # Set custom colors for the fill
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "PB/PBH vs Control - Spat",
       x = "",
       y = "Log-Fold Change in Relative Abundance") +
  theme_bw() +
  coord_flip() +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.text.x = element_text(angle = 0, hjust = 1, size = 9),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 9)) 
#+ ylim(0, 2.1)

