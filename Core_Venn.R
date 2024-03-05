#venn diagrams for core microbiome data
#source: https://microbiome.github.io/tutorials/core_venn.html

#Packages ----

install.packages("eulerr")
install.packages("microbiomeutilities")

library(eulerr)
library(microbiome)
devtools::install_github('microsud/microbiomeutilities')
library(microbiomeutilities)

pseq<- Marissa_mb2021_filtered_20240203

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree


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

#Remove day 3 (only 1 sample remaining) for mb2021 project

pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))


#Compositional ----

pseq.rel <- microbiome::transform(pseq, "compositional")

#Get variable to analyze
Treatment <- unique(as.character(Metadata$Treatment))
print(Treatment)


#Write a for loop to go through each of the Treatments one by one and combine identified core taxa into a list

list_core <- c() # an empty object to store information

for (n in Treatment){ # for each variable n in Treatment
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, Treatment == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.90)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)

#With all timepts very few core taxa

#[1] "No. of core taxa in High salinity : 2"
#[1] "No. of core taxa in Control : 5"
#[1] "No. of core taxa in Low salinity : 4"

#$`High salinity`
#[1] "ASV1"  "ASV13"

#$Control
#[1] "ASV1"  "ASV3"  "ASV8"  "ASV13" "ASV30"

#$`Low salinity`
#[1] "ASV1"  "ASV8"  "ASV13" "ASV40"

mycols <- c(nonCRC="#d6e2e9", CRC="#cbf3f0", H="#fcf5c7") 
mycols <- c(nonCRC="lightblue", CRC="lightgreen", H="violet") 

plot(venn(list_core),
     fills = mycols)

for (name in names(list_core)) {
  filename <- paste0(name, "venn_core.txt")  # Create a filename based on the list element name
  writeLines(list_core[[name]], filename)  # Save the list element to a text file
}

#Age

Age <- unique(as.character(Metadata$Age))
print(Age)


#Write a for loop to go through each of the Treatments one by one and combine identified core taxa into a list

list_core <- c() # an empty object to store information

for (n in Age){ # for each variable n in Treatment
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, Age == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.90)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)

mycols <- c(nonCRC="#d6e2e9", CRC="#cbf3f0", H="#fcf5c7") 
mycols <- c(nonCRC="lightblue", CRC="lightgreen", H="violet") 

plot(venn(list_core),
     fills = mycols)


#only ASV 1 and 13 shared between all timepts (Rhodo and alteromonas)

#To get taxa info with ASV # combined

# use the pseq.rel object created at the begening of this tutorial. 
taxa_names(pseq.rel)[1:5]

# format names
pseq.rel.f <- format_to_besthit(pseq.rel)
# check names
taxa_names(pseq.rel.f)[1:5]

#Run the for loop to go through each of the Treatments one by one and combine identified core taxa into a list.

list_core <- c() # an empty object to store information

for (n in Age){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel.f, Age == n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.90)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

View(list_core)
print(list_core)

#Save list of ASVs+Tax info for each group

# Iterate over each list element and save it to a separate text file
for (name in names(list_core)) {
  filename <- paste0(name, "venn_core.txt")  # Create a filename based on the list element name
  writeLines(list_core[[name]], filename)  # Save the list element to a text file
}



