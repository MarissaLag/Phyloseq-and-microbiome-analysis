# Here's all the info you should need for running DEICODE in windows!
# Parker K Lund
# August 10, 2023; Updated December 17, 2024

#install.packages("tibble")
library(tibble)
#install.packages("phyloseq")
library(phyloseq)

#Edited for Marissa's data December, 2024

#Open data
MU42022_filtered_Oct92024 <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/MU42022_filtered_Oct92024.rds")
phylo_obj <- MU42022_filtered_Oct92024

directory <- "/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis"
marker <- "MU42022_rare"


# Export the files you need for QIIME2/DEICODE from your phyloseq object in R (ASV count table, taxonomy table, and metadata table)

phylo_obj_asv_table <- as.data.frame(otu_table(phylo_obj))
phylo_obj_asv_table <- rownames_to_column(phylo_obj_asv_table,"#OTUID")
write.table(phylo_obj_asv_table,paste0(directory,"/","phylo_obj_",marker,"_asv_table_for_biom.txt"),sep = "\t",row.names = FALSE,quote = FALSE)
phylo_obj_tax_table <- as.data.frame(as(tax_table(phylo_obj),"matrix"))
phylo_obj_tax_table <- rownames_to_column(phylo_obj_tax_table,"#OTUID")
write.table(phylo_obj_tax_table,paste0(directory,"/","phylo_obj_",marker,"_tax_table_for_biom.txt"),sep = "\t",row.names = FALSE,quote = FALSE)
phylo_obj_metadata <- as.data.frame(as.matrix(sample_data(phylo_obj)))
phylo_obj_metadata <- rownames_to_column(phylo_obj_metadata,"#SampleID")
write.table(phylo_obj_metadata,paste0(directory,"/","phylo_obj_",marker,"_metadata_for_biom.txt"),sep = "\t",row.names = FALSE,quote = FALSE)

# TURN ON WINDOWS SUBSYSTEM LINUX #
1. Open Windows Features through the search bar in Windows 10.
2. Click the checkbox next to Windows Subsystem for Linux.
3. Restart computer. 
4. Download Ubuntu from the Windows 10 app store. 
5. Open Ubuntu.
6. Activate the ability to copy-paste between Windows 10 and Ubuntu by right-clicking in the upper-left corner, selecting Properties, and clicking Use Ctrl+Shift+C/V as Copy/Paste.

# DOWNLOAD MINICONDA #
1. From within Ubuntu, download Miniconda for Linux.
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
2. Once Miniconda is activated, close Ubuntu so we can restart WSL.
3. Open Windows Features and unclick Windows Subsystem for Linux. Restart computer.
4. Open Windows Features and re-click Windows Subsystem for Linux. Restart computer. 

# Download qiime2
wget https://data.qiime2.org/distro/core/qiime2-2023.5-py38-linux-conda.yml

# If there is an environment active, even (base), input
conda deactivate

# Create conda environment (if the connection error occurs, you may need to go into windows features and restart the windows linux subsystem)
conda env create -n qiime2-2023.5 --file qiime2-2023.5-py38-linux-conda.yml

# Activate Qiime2 environment
conda activate qiime2-2023.5

#Now install Deicode within the Qiime2 environment!
#Deicode installation info from: https://forum.qiime2.org/t/robust-aitchison-pca-beta-diversity-with-deicode/8333
conda install -c conda-forge deicode

# Convert .txt files to .biom files
# Biom conversion in Qiime2: http://biom-format.org/documentation/biom_conversion.html
# Note: -i denotes the input .txt file of your ASV table, -o denotes the output biom file.
biom convert -i /mnt/c/Users/Parker/Desktop/MBARI/Analysis/THL/COI/DeicodeFiles/phylo_obj_env_COI_asv_table_for_biom.txt -o table.from_txt_json.biom --table-type="OTU table" --to-json

# Add taxonomy table and metadata .txt to biom file
# Note: -i denotes the biom file just creates, -o denotes the output biom file, --observation-metadata-fp denotes the taxonomy table, and --sample-metadata-fp denotes the metadata about all of the samples (temperature, depth, etc.)
biom add-metadata -i table.from_txt_json.biom -o table.w_md.biom --observation-metadata-fp /mnt/c/Users/Parker/Desktop/MBARI/Analysis/THL/COI/DeicodeFiles/phylo_obj_env_COI_tax_table_for_biom.txt --sample-metadata-fp /mnt/c/Users/Parker/Desktop/MBARI/Analysis/THL/COI/DeicodeFiles/phylo_obj_env_COI_metadata_for_biom.txt

# Import into Qiime2 - Change output file.
qiime tools import \
--input-path table.w_md.biom \
--output-path master.biom.qza \
--type FeatureTable[Frequency]

# Run DEICODE - Change input table.
qiime deicode rpca \
--i-table master.biom.qza \
--p-n-components 3 \
--p-min-feature-count 20 \
--p-min-sample-count 500 \
--o-biplot ordination.qza \
--o-distance-matrix distance.qza

# Create biplot - Change metadata files (number of features is number of taxonomic levels).
qiime emperor biplot \
--i-biplot ordination.qza \
--m-sample-metadata-file /mnt/c/Users/Parker/Desktop/MBARI/Analysis/THL/COI/DeicodeFiles/phylo_obj_env_COI_metadata_for_biom.txt \
--m-feature-metadata-file /mnt/c/Users/Parker/Desktop/MBARI/Analysis/THL/COI/DeicodeFiles/phylo_obj_env_COI_tax_table_for_biom.txt \
--o-visualization biplot.qzv \
--p-number-of-features 8


# Run a permanova | permdisp | or anosim test between whatever metadata column you choose (argument --m-metadata-column [column of choosing])
# Here I ran an ANOSIM

qiime diversity beta-group-significance \
--i-distance-matrix distance.qza \
--m-metadata-file /mnt/c/Users/Parker/Desktop/MBARI/Analysis/THL/COI/DeicodeFiles/phylo_obj_env_COI_metadata_for_biom.txt \
--m-metadata-column SAMPLING_station \
--p-method anosim \
--o-visualization SAMPLING_station_significance_anosim.qzv

# You can make a numeric data column 
into categorical (i.e. month) to run a permanova/anosim in Qiime2.
# Add a new row to metadata .txt file. The row’s first cell must be #q2:types to indicate the row is a comment directive. Subsequent cells may contain the values categorical or numeric (both case-insensitive). The empty cell is also supported if you do not wish to assign a type to a column (the type will be inferred in that case). Thus, it is easy to include this comment directive without having to declare types for every column in your metadata. Qiime2 CANNOT convert to numeric data if there are spots listed as NA, it will read these as characters. 

# All of these .qza files can be viewed using the online viewer! https://view.qiime2.org/
# Just drag and drop the files into your browser. This is actually pretty cool for the biplot because you can view it in 3D space and rotate it (coloring the points by some variable).


# Now we go back to R!!

################ PLOT PRINCIPAL COMPONENTS ################

library(qiime2R)

# Import Qiime2 Results (this is the biplot)
file = paste0(directory,marker,"/DeicodeFiles/Output/","biplot.qza")
pco=read_qza(file)
#look at data
head(pco$data$ProportionExplained)
pco$data$Vectors[1:5, 1:4]

# Create "proportion explained" labels for the top three PC axes
label.PC1 <- paste("PC1: ", round(pco$data$ProportionExplained$PC1, 3)*100,"%")
label.PC1
label.PC2 <- paste("PC2: ", round(pco$data$ProportionExplained$PC2, 3)*100,"%")
label.PC2
label.PC3 <- paste("PC3: ", round(pco$data$ProportionExplained$PC3, 3)*100,"%")
label.PC3

# Join with sample data (it's written here as sample_dat, but it should be a dataframe containing the sample_data from your phyloseq object)
pco$data$Vectors <- pco$data$Vectors %>% 
  rename(., "sample_name" = "SampleID")
pcscores <- left_join(pco$data$Vectors, sample_dat, by="sample_name")

# Format loading scores. The "ASV" column should be labeled whatever you have labeled that column in your otu_table and taxa_table from your phyloseq object. 
loadings <- as.data.frame(pco$data$Species)
loadings$ASV <- loadings$FeatureID

# Join with taxa table from phyloseq object (written here as taxa_tab). Same as the previous step, "ASV" should be whatever it is labeled for your data.
loadings <- left_join(loadings, taxa_tab, by="ASV")

# Export pcscores
file = paste0(directory,marker,"/Figures/PCA/",marker,"_pcscores_Qiime2_env.csv",sep="")
print(file)
write.csv(pcscores, file)


# Plot PCA Scores

library(ggthemes) # This is for scale_fill_tableau, not mandatory.

p<- pcscores %>% 
  ggplot(aes(as.factor(season), PC1 ,color=factor(season)))+geom_boxplot(lwd=1) +
  scale_fill_tableau(palette = "Tableau 10", type = c("regular"), direction = 1) +
  labs(x="season", y=label.PC1, colour="season", subtitle="environmental samples")+ ggtitle(paste(project, marker, "PC1 boxplot"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        axis.text.x=element_text(size=8, angle = 45),
        legend.title = element_text(face = "bold", size=12),
        legend.key = element_rect(fill = FALSE, colour = FALSE),
        # legend.key.size = unit(0.1,"line")
  )+geom_jitter(alpha=.3)

# This example is only plotting the scores for PC1 by Season. You can change it to PC2 or PC3. 


################ LOADING SCORES ################

# Group Loading Scores by Taxonomic Level

loadings %>% 
  full_join(species_label) %>% #join with taxonomy
  #merge some groups:
  mutate(Genus = case_when(Genus=='g_'| Genus =='unassigned' | Genus =='unknown'~as.character('unknown'),
                           TRUE ~ as.character(Genus)))

filename = paste(directory,marker,"/Figures/PCA/",marker,'_loading_scores_env.csv', sep='')
filename
write_csv(loadings,filename)


# PLOT LOADING SCORES

# Extract the top 25 and bottom 25 loading scores for PC1 and then combine them in one dataframe. arrange(desc(PC1)) arranges them in descending numerical order, and the following mutation factors FeatureID by the current order (so that it is in descending order when plotted).

# PC1
load_top = loadings %>%
  slice_max(PC1, n=10)
load_bottom = loadings %>%
  slice_min(PC1, n=10)
load_combo = load_top %>%
  full_join(load_bottom) %>%
  arrange(desc(PC1)) %>%
  mutate(FeatureID = factor(FeatureID, unique(FeatureID)))


load = ggplot(data=load_combo, aes(PC1, FeatureID)) + geom_col(aes(fill=Class)) +
  labs(x="PC1", y="ASV", colour="Family", subtitle="all environmental samples") + 
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1) +
  ggtitle(paste(project, marker, "PC1 Loading Scores - Depth")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        axis.text.x=element_text(size=13, angle = 90),
        axis.text.y=element_text(size=13),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(face = "bold", size=12),
        legend.key = element_rect(fill = FALSE, colour = FALSE),
  )

# PC2
load_top = loadings %>%
  slice_max(PC2, n=10)
load_bottom = loadings %>%
  slice_min(PC2, n=10)
load_combo = load_top %>%
  full_join(load_bottom) %>%
  arrange(desc(PC2)) %>%
  mutate(FeatureID = factor(FeatureID, unique(FeatureID)))


load = ggplot(data=load_combo, aes(PC2, FeatureID)) + geom_col(aes(fill=Class)) +
  labs(x="PC2", y="ASV", colour="Family", subtitle="all environmental samples") + 
  scale_fill_tableau(palette = "Tableau 20", type = c("regular"), direction = 1) +
  ggtitle(paste(project, marker, "PC2 Loading Scores - Season")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        axis.text.x=element_text(size=13, angle = 90),
        axis.text.y=element_text(size=13),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        legend.title = element_text(face = "bold", size=12),
        legend.key = element_rect(fill = FALSE, colour = FALSE),
  )


################ PLOT RPCA ORDINATION ################

p<- ggplot(pcscores,aes(PC1,PC2,color=factor(season)))+geom_point(size=3, alpha=0.8) +
  labs(x=label.PC1 , y=label.PC2, colour = "Season")+ ggtitle(paste(project, marker, "PCA"))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.text = element_text(size = 11, colour = "black"),
        axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"),
        plot.subtitle=element_text(size=9),
        legend.title = element_text(face = "bold", size=12),
        legend.key = element_rect(fill = FALSE, colour = FALSE),
  )+ guides(color = guide_legend(override.aes = list(size=5, shape=15)))