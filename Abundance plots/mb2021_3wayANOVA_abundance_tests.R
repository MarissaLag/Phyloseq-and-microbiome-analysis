##t-tests of SAR bacterial isolates from mb2021 - control and high salinity treatments
##sorce: https://www.statology.org/three-way-anova-in-r/

install.packages("datarium")
install.packages("rstatix")
install.packages("ggpubr")
install.packages("tidyverse")

library("datarium")
library(tidyverse)
library(ggpubr)
library(rstatix)
library("devtools")
library(phyloseq)
library(microbiome)
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggplot2)
library(ggResidpanel)

#load data
pseq<- Marissa_mb2021_filtered_20240203
#filter data (if needed)
pseq <- subset_samples(pseq, !Age %in% c("3 dpf"))
pseq <- subset_samples(pseq, Age %in% c("1 dpf"))

#convert to compositional (always recommended)
pseq <- microbiome::transform(pseq, "compositional")

#extract taxa you would like to test
#For mb2021 I want to test for differences of Flavobacteriaceae, 
#Rhodobactereae, and Vibrionaceae at family level

pseq.selected <- subset_taxa(pseq, Family %in% c("Vibrionaceae", "Flavobacteriaceae", "Rhodobacteraceae"))

#convert to data frame
pseq_psmelt <- psmelt(pseq)
View(pseq_psmelt)

set.seed(123)

#get means of groups ----
pseq_psmelt %>% sample_n_by(Treatment, size = 3)

pseq_psmelt %>%
  group_by(Family, Treatment) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE),
            sd_abundance = sd(Abundance, na.rm = TRUE))

View(pseq_psmelt)

#test each seperately
Data <- pseq_psmelt %>%
  filter(Family == "Vibrionaceae")

Data <- pseq_psmelt %>%
  filter(Family == "Flavobacteriaceae")

vibrionaceae_data <- subset(Data, Family == "Vibrionaceae", select = c("Abundance", "Treatment"))

Data <- pseq_psmelt %>%
  filter(Family == "Rhodobacteraceae")

View(Data)


#perform ANOVA ----
model <- aov(Abundance ~ Treatment*sample_Family, data=Data)
#view summary of three-way ANOVA
summary(model)


#Another way
vibrionaceae_data <- subset(Data, Family == "Vibrionaceae", select = c("Abundance", "Treatment"))
flavo_data <- subset(Data, Family == "Flavobacteriaceae", select = c("Abundance", "Treatment"))

# Perform ANOVA
anova_result <- aov(Abundance ~ Treatment, data = vibrionaceae_data)

# Summarize ANOVA results
summary(anova_result)

#non-paramtric test
#Kruskal-Wallis test
kruskal_result <- kruskal.test(Abundance ~ Treatment, data = vibrionaceae_data)

print(kruskal_result)

#perform non parametric t-test on average values

#compare HS and control only

pseq_psmelt <- psmelt(pseq)

vibrionaceae_data <- subset(pseq_psmelt, Treatment %in% c("High salinity", "Control"))
vibrionaceae_data <- subset(vibrionaceae_data, Family == "Vibrionaceae")

#Wilcoxon rank-sum test
wilcox_test_result <- wilcox.test(Abundance ~ Treatment, data = vibrionaceae_data)

print(wilcox_test_result)



#if need to change variable type
str(Data)

Data$Family_Number <- as.character(Data$Family_Number)
Data$Size <- as.numeric(Data$Size)





##rhodo is signif different with Age (p = 0.00165), and genetics (0.00133)

Df  Sum Sq Mean Sq F value  Pr(>F)   
Age                     3 1709700  569900   6.417 0.00165 **
  Treatment               2   77578   38789   0.437 0.65003   
Genetics                1 1104840 1104840  12.441 0.00133 **
  Age:Treatment           4  312374   78093   0.879 0.48758   
Age:Genetics            2  137140   68570   0.772 0.47072   
Treatment:Genetics      2  139783   69891   0.787 0.46410   
Age:Treatment:Genetics  4   74337   18584   0.209 0.93133   
Residuals              31 2753049   88808                   

Data %>%
  group_by(Age) %>%
  summarize(Rhodobacteraceae = mean(Rhodobacteraceae))

Data %>%
  group_by(Genetics) %>%
  summarize(Rhodobacteraceae = mean(Rhodobacteraceae))

#results

Age    Rhodobacteraceae
<chr>             <dbl>
  1 day 1             1241.
2 day 18             905.
3 day 3             1221 
4 spat              1396.
> Data %>%
  +     group_by(Genetics) %>%
  +     summarize(Rhodobacteraceae = mean(Rhodobacteraceae))
# A tibble: 4 × 2
Genetics Rhodobacteraceae
<int>            <dbl>
  1        1            1599.
2        2            1187.
3        3            1204.
4        4            1146.

#perform three-way ANOVA
model <- aov(Flavobacteriaceae ~ Age * Treatment * Genetics, data=Data)

#view summary of three-way ANOVA
summary(model)

##results

Df  Sum Sq Mean Sq F value   Pr(>F)    
Age                     3 5787073 1929024  45.179 1.99e-11 ***
  Treatment               2  582309  291155   6.819 0.003513 ** 
  Genetics                1  137450  137450   3.219 0.082538 .  
Age:Treatment           4 2613390  653348  15.302 5.19e-07 ***
  Age:Genetics            2  802714  401357   9.400 0.000644 ***
  Treatment:Genetics      2   60433   30216   0.708 0.500567    
Age:Treatment:Genetics  4  342547   85637   2.006 0.118181    
Residuals              31 1323604   42697               

Data %>%
  group_by(Age) %>%
  summarize(Flavobacteriaceae = mean(Flavobacteriaceae))

Data %>%
  group_by(Treatment) %>%
  summarize(Flavobacteriaceae = mean(Flavobacteriaceae))


#results
Age    Flavobacteriaceae
<chr>              <dbl>
  1 day 1              1395.
2 day 18              376 
3 day 3               725 
4 spat                503.

Treatment     Flavobacteriaceae
<chr>                     <dbl>
  1 Control                    539.
2 High salinity              788.
3 Low salinity               553.


#boxplot source: https://www.statmethods.net/graphs/boxplot.html 

boxplot(Flavobacteriaceae ~ Treatment,data=Data, main="Flavobacteriaceae",
        xlab="Treatment", ylab="Abundance", ylim = c(0, 1500))

##ylim is to set y axis limits, xlim for x axis
##names function is to set x-label names
#las function is to rotate x-labels to vertical direction (not enough space horizontally)

boxplot(Rhodobacteraceae ~ Treatment*Age,data=Data, main="Rhodobacteraceae", names = c("Control, day 1", "HS, day 1", "LS, day 1", "Control, day 18", "HS, day 18", "LS, day 18", "Control, day 3", "HS, day 3", "LS, day 3", "Control, spat", "HS, spat", "LS, spat"), xlab="", ylab="Abundance", las=2)

boxplot(Flavobacteriaceae ~ Age*Treatment,data=Data, main="Flavobacteriaceae",
        xlab="Treatment", ylab="Abundance")

boxplot(Flavobacteriaceae ~ Treatment*Age,data=Data, main="Flavobacteriaceae", names = c("Control, day 1", "HS, day 1", "LS, day 1", "Control, day 18", "HS, day 18", "LS, day 18", "Control, day 3", "HS, day 3", "LS, day 3", "Control, spat", "HS, spat", "LS, spat"), xlab="", ylab="Abundance", las=2)

***on day 1, high sal treatment has much higher abundance of flavo

##confused why rhodo does not have much differences on boxplot... looks much higher on other plot.


###now checking Vibrio abundance (on different data sheet)

setwd("~/USRA2021/mb2021/m2021_new")

#load data

Data=read.csv("Vibrio_abundance.csv")
View(Data)



set.seed(123)
Data %>% sample_n_by(Treatment, size = 3)

Data %>%
  group_by(Treatment) %>%
  get_summary_stats(Vibrionaceae, type = "mean_sd")

##results

Treatment     variable         n  mean    sd
<chr>         <fct>        <dbl> <dbl> <dbl>
  1 Control       Vibrionaceae    18  32.8  43.0
2 High salinity Vibrionaceae    17  91.1 270. 
3 Low salinity  Vibrionaceae    15  62.5  93.3


#perform three-way ANOVA
model <- aov((Vibrionaceae)) ~ Age * Treatment * Genetics, data=Data)

#view summary of three-way ANOVA
summary(model)

##results - Vibrio is signif different with Age (p = 0.00722)

Age                     3 356803  118934   4.819 0.00722 **
  Treatment               2  32492   16246   0.658 0.52485   
Genetics                1   2318    2318   0.094 0.76132   
Age:Treatment           4 174402   43601   1.767 0.16080   
Age:Genetics            2   7321    3661   0.148 0.86277   
Treatment:Genetics      2    917     459   0.019 0.98160   
Age:Treatment:Genetics  4   5860    1465   0.059 0.99313   
Residuals              31 765124   24681                   

Data %>%
  group_by(Age) %>%
  summarize(Vibrionaceae = mean(Vibrionaceae))

#results

Age    Vibrionaceae
<chr>         <dbl>
  1 day 1        213.  
2 day 18       104.  
3 day 3        289   
4 spat           4.41


#boxplot source: https://www.statmethods.net/graphs/boxplot.html 

boxplot(Vibrionaceae ~ Treatment*Age,data=Data, main="Vibrionaceae", names = c("Control, day 1", "HS, day 1", "LS, day 1", "Control, day 18", "HS, day 18", "LS, day 18", "Control, day 3", "HS, day 3", "LS, day 3", "Control, spat", "HS, spat", "LS, spat"), xlab="", ylab="Abundance", las=2)


##F2H1 appears to be large outlier = can I remove?

filtered_data <- subset(Data, Sample != "F2H1")

# Create the boxplot with the filtered data
boxplot(Vibrionaceae ~ Treatment*Age, data = filtered_data, main = "Vibrionaceae",
        names = c("Control, day 1", "HS, day 1", "LS, day 1", "Control, day 18", "HS, day 18", "LS, day 18",
                  "Control, day 3", "HS, day 3", "LS, day 3", "Control, spat", "HS, spat", "LS, spat"),
        xlab = "", ylab = "Abundance", las = 2)

#test for differences without F2H1

filtered_data <- subset(Data, Sample != "F2H1")

# Fit the model with the filtered data
model <- aov(Vibrionaceae ~ Age * Treatment * Genetics, data = filtered_data)

Df Sum Sq Mean Sq F value   Pr(>F)    
Age                     3 152646   50882 143.204  < 2e-16 ***
  Treatment               2   5347    2673   7.524  0.00225 ** 
  Genetics                1    136     136   0.382  0.54104    
Age:Treatment           4  13519    3380   9.512 4.34e-05 ***
  Age:Genetics            2   2511    1256   3.534  0.04186 *  
  Treatment:Genetics      2    197      99   0.278  0.75938    
Age:Treatment:Genetics  4   3938     984   2.770  0.04519 *  
  Residuals              30  10659     355                   






##test Rhodo abundance in spat only

Data=read.csv("Rhodo_spat.csv")
View(Data)

set.seed(123)
Data %>% sample_n_by(Treatment, size = 3)

Data %>%
  group_by(Treatment) %>%
  get_summary_stats(Rhodobacteraceae, type = "mean_sd")

##results
Treatment     variable             n  mean    sd
<chr>         <fct>            <dbl> <dbl> <dbl>
  1 Control       Rhodobacteraceae    11 1348.  385.
2 High salinity Rhodobacteraceae    11 1499.  361.
3 Low salinity  Rhodobacteraceae    10 1335.  299.

#perform three-way ANOVA
model <- aov(Rhodobacteraceae ~ Treatment * Genetics, data=Data)

#view summary of two-way ANOVA
summary(model)

#results - signif genetics
Df  Sum Sq Mean Sq F value  Pr(>F)   
Treatment           2  179058   89529   0.891 0.42238   
Genetics            1  811513  811513   8.077 0.00861 **
  Treatment:Genetics  2  163990   81995   0.816 0.45317   
Residuals          26 2612367  100476                   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

set.seed(123)
Data %>% sample_n_by(Genetics, size = 4)
Data %>%
  group_by(Genetics) %>%
  get_summary_stats(Rhodobacteraceae, type = "mean_sd")

Genetics variable             n  mean    sd
<int> <fct>            <dbl> <dbl> <dbl>
  1        1 Rhodobacteraceae     8 1790   410.
2        2 Rhodobacteraceae     8 1236.  221.
3        3 Rhodobacteraceae     8 1253   250.
4        4 Rhodobacteraceae     8 1305.  138.

boxplot(Rhodobacteraceae ~ Genetics,data=Data, main="Rhodobacteraceae", xlab="", ylab="Abundance", las=2)
boxplot(Rhodobacteraceae ~ Genetics,data=Data, main="Rhodobacteraceae", names = c("F1", "F2", "F3", "F4"), xlab="", ylab="Abundance", las=2)
boxplot(Rhodobacteraceae ~ Genetics*Treatment,data=Data, main="Rhodobacteraceae", names = c("F1, Control", "F2, Control", "F3, Control", "F4, Control", "F1, HS", "F2, HS", "F3, HS", "F4, HS", "F1, LS", "F2, LS", "F3, LS", "F4, LS"), xlab="", ylab="Abundance", las=2)

