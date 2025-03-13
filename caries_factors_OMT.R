############################################################################################
############## FACTORS ASSOCIATED WITH DONOR SELECTION FOR CARIES MODEL ##############
####################################################################################
##Written by Sonia Nath
## Date:10/02/25
#Clear existing data and graphics
rm(list=ls())
graphics.off()

#SETTING THE DIRECTORY
setwd("~/DIET_CARIES/All_data/")

#Library the packages
library(DECIPHER) # This package will help in importing, maintaining, analyzing, manipulating, and exporting a massive amount of sequences.

library(ape) # Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq

library(DESeq2) # This package will help analyze "differential expression" in the microbiota alongside phyloseq

library(ggplot2) # Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.

library(phyloseq) # The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.

library(plotly) # A package to create interactive web graphics of use in 3D plots

library(vegan) # The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.

library(philr) # This package provides functions for the analysis of compositional data 

library(tidyverse) # This package is designed to make it easy to install and load multiple 'tidyverse' packages in a single step

library(adespatial) # Tools for the multiscale spatial analysis of multivariate data

library(devtools) # Make package development easier by providing R functions that simplify and expedite common tasks

library(qiime2R) # A package for importing qiime artifacts into an R session

library(MicrobeR) # Data visualization

library(microbiome) # Data analysis and visualization

library("pander") # provide a minimal and easy tool for rendering R objects into Pandoc's markdown

library(grid) # support data visualization

library(gridExtra)  # support data visualization

library(knitr) # Provides a general-purpose tool for dynamic report generation in R using Literate Programming techniques.

library(png) # Figure download

library("ggdendro") #set of tools for dendrograms and tree plots using 'ggplot2'

library(ggpubr) # publication quality figures, based on ggplot2

library(RColorBrewer) # nice color options

library(DT) #for interactive tables

library(reshape2)

library(scales)

library(data.table)

library(Biostrings)
library(readr)

#UPLOADING THE DATASET###
metadata <-read_tsv("metadata.tsv")

####CLEANIG UP THE VARAIBLES########
#1. Energy consumption
#2. Water consumption
#3. Carbohydrate consumption
#4. Sugar consumption
#5. Caries score- NO NEEDED TO TRANFORM
#6. Salivary pH
#7. Salivary flow rate
#8. Professionally applied fluoride

# Calculate the median of the continuous variable
median_value <- median(metadata$Energy.inc.fibre)
# Dichotomize based on the median
metadata$Energy_di <- ifelse(metadata$Energy.inc.fibre > median_value, "High energy consumption", 
                             "Low energy consumption")

# Calculate the median of the continuous variable
median_value <- median(metadata$Water)
# Dichotomize based on the median
metadata$Water_di <- ifelse(metadata$Water > median_value, "High water consumption", 
                            "Low water consumption")median_value <- median(metadata$Sugars)

# Dichotomize based on the median
metadata$Sugars_di <- ifelse(metadata$Sugars > median_value, "High sugar", "Low sugar")

median_value <- median(metadata$Carbohydrate)
# Dichotomize based on the median
metadata$Carbohydrate_di <- ifelse(metadata$Carbohydrate > median_value, "High carbohydrate", "Low carbohydrate")

# Calculate the median of the continuous variable
median_value <- median(metadata$saliva_flow_rate)

# Dichotomize based on the median
metadata$Saliva_flow_di <- ifelse(metadata$saliva_flow_rate > median_value, "High saliva flow rate", 
                                  "Low saliva flow rate")

# Calculate the median of the continuous variable
median_value <- median(metadata$post.glucose_2)
# Dichotomize based on the median
metadata$post.gluocose_di <- ifelse(metadata$post.glucose_2 > median_value, "High pH", "Low pH")

metadata$fluoride <- ifelse(metadata$fluoride_prof_status== "Fluoride applied", 1,0)

write.table(metadata,
            file = "metadata_phy.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

####################################
####Creating a phyloseq object#####
###################################
physeq<-qza_to_phyloseq(
  features="table.qza",
  tree="rooted-tree.qza",
  "taxonomy.qza", 
  metadata= "metadata_phy.tsv"
)
physeq

###Summarising the phyloseq object
summarize_phyloseq(physeq) #gives the brief summary of the samples and the data
print_ps(physeq) #prints the same summary
summary(sample_sums(physeq)) #fives the min, max, median, mean for the reads
sample_names(physeq)[1:5] #overview of sample names 
rank_names(physeq)  #Rank names
sample_variables(physeq)  # the variables in the metadata

otu_table(physeq)[1:5, 1:5] ## the feature table 

tax_table(physeq)[1:5, 1:4]  # the taxonomy table

ps <- physeq
######Alpha divrsity of categorical variables########
#1. Energy consumption
#2. Water consumption
#3. Carbohydrate consumption
#4. Sugar consumption
#5. Caries score
#6. Salivary pH
#7. Salivary flow rate
#8. Professionally applied fluoride

#####    ALPHA DIVERSITY #############
#########################################
richness <- estimate_richness(ps)
head(richness)
alpha_diversity <- data.frame(sample_data(ps), richness)

write.table(richness,
            file = "richness.tsv", 
           sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(alpha_diversity,
            file = "alpha_diversity.tsv", 
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

###Check the normal distribution curve for alpha diversity
hist(richness$Shannon, main="Shannon index", xlab="")
hist(richness$Observed, main="Observed features", xlab="")
hist(richness$Chao1, main="Chao 1 index", xlab="")

# Perform Shapiro-Wilk normality test
shapiro.test(alpha_diversity$Shannon)
shapiro.test(alpha_diversity$Observed)
shapiro.test(alpha_diversity$Chao1)

symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

#####1. PLOT

####### Caries Experience #######
levels <- unique(sample_data(ps)$caries_experience_status)
comps <- combn(levels, 2, simplify = FALSE)

comps <- make_pairs(sample_data(ps)$caries_experience_status)
print(comps)

kruskal.test(richness$Shannon ~ sample_data(ps)$caries_experience_status)
pairwise.wilcox.test(richness$Shannon, sample_data(ps)$caries_experience_status, p.adj = "bonf")

t.test(Observed~ caries_experience_status, data = alpha_diversity)

p1 <- plot_richness(ps, x = "caries_experience_status", measures = c("Observed", "Shannon"), 
                    color = "caries_experience_status") +
  geom_boxplot(aes(fill = caries_experience_status), color = "black", outlier.color = "black") +  # Specifying black color for outline
  scale_fill_manual(values = c("Caries Experience" = "#0072B2", "No Caries Experience" = "#D55E00")) +  # Specify actual colors
  labs(x = "Past Caries Experience") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = caries_experience_status), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)


# Print the plot
print(p1)

####### Energy #######
kruskal.test(richness$Observed ~ sample_data(ps)$Energy_di)
kruskal.test(richness$Shannon ~ sample_data(ps)$Energy_di)
comps <- make_pairs(sample_data(ps)$Energy_di)

p2 <- plot_richness(ps, x = "Energy_di", measures = c("Observed", "Shannon"), 
                    color = "Energy_di") +
  geom_boxplot(aes(fill = Energy_di), color = "black", outlier.color = "black") +  # Specifying black color for outline
  scale_fill_manual(values = c("Low energy consumption" = "#0072B2", "High energy consumption" = "#D55E00")) +  # Specify actual colors
  labs(x = "Total energy consumption per day") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = Energy_di), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

print(p2)

####### Water consumption. #######
kruskal.test(richness$Observed ~ sample_data(ps)$Water_di)
kruskal.test(richness$Shannon ~ sample_data(ps)$Water_di)


p3 <- plot_richness(ps, x = "Water_di", measures = c("Observed", "Shannon"), 
                    color = "Water_di") +
  geom_boxplot(aes(fill = Water_di), color = "black") +  # Removed outlier.color since it's not a valid parameter for violin
  scale_fill_manual(values = c("High water consumption" = "#0072B2", "Low water consumption" = "#D55E00")) +  # Specify actual colors
  labs(x = "Water consumption in a day") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = Water_di), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, method = "wilcox.test", label = "p.signif", exact = FALSE,
                     label.y = c(5, 7, 6), symnum.args = symnum.args) +  # Wilcoxon test
  stat_compare_means(label.y = 7)  # Overall comparison

# Print the plot
print(p3)


####### Carbohydrate #######
kruskal.test(richness$Observed ~ sample_data(ps)$Carbohydrate_di)
kruskal.test(richness$Shannon ~ sample_data(ps)$Carbohydrate_di)
comps <- make_pairs(sample_data(ps)$Carbohydrate_di)

p4 <- plot_richness(ps, x = "Carbohydrate_di", measures = c("Observed", "Shannon"), 
                    color = "Carbohydrate_di") +
  geom_boxplot(aes(fill = Carbohydrate_di), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("High carbohydrate" = "#882255", "Low carbohydrate" = "#005AB5")) +  # Specify actual colors
  labs(x = "Carbohydrate consumption/day in gms") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = Carbohydrate_di), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", exact = FALSE, symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)
print(p4)


####### Sugars #######
kruskal.test(richness$Observed ~ sample_data(ps)$Sugars_di)
kruskal.test(richness$Shannon ~ sample_data(ps)$Sugars_di)
comps <- make_pairs(sample_data(ps)$Sugars_di)

p5 <- plot_richness(ps, x = "Sugars_di", measures = c("Observed", "Shannon"), 
                    color = "Sugars_di") +
  geom_boxplot(aes(fill = Sugars_di), color = "black", outlier.color = NA) +  # Specifying black color for outline
  scale_fill_manual(values = c("High sugar" ="#882255", "Low sugar" = "#005AB5")) +  # Specify actual colors
  labs(x = "Sugar consumption/day in gms") +
  theme_minimal() +  # Using a minimal theme for a clean look
  theme(strip.background = element_blank(),  # Removing grey background completely
        legend.position = "none") +
  geom_jitter(aes(color = Sugars_di), position = position_jitter(width = 0.2)) +  # Added jitter
  stat_compare_means(comparisons = comps, label.y = c(5, 7, 6), method = "wilcox.test",  
                     label = "p.signif", symnum.args = symnum.args) +
  stat_compare_means(label.y = 7)

print(p5)

##### Saliva flow rate #######
kruskal.test(richness$Observed ~ sample_data(ps)$Saliva_flow_di)
kruskal.test(richness$Shannon ~ sample_data(ps)$Saliva_flow_di)
comps <- make_pairs(sample_data(ps)$Saliva_flow_di)

p6 <- plot_richness(ps, x = "Saliva_flow_di", measures = c("Observed", "Shannon"), 
                    color = "Saliva_flow_di") +
  geom_boxplot(aes(fill = Saliva_flow_di), color = "black") + 
  scale_fill_manual(values = c("High saliva flow rate" = "#0072B2", "Low saliva flow rate" = "#D55E00")) + 
  labs(x = "Saliva flow rate") +
  theme_minimal() + 
  theme(strip.background = element_blank(),
        legend.position = "none") +
  geom_jitter(aes(color = Saliva_flow_di), position = position_jitter(width = 0.2)) +
  stat_compare_means(comparisons = comps, method = "wilcox.test", 
                     label = "p.signif", exact = FALSE, label.y = c(5, 7, 6), 
                     symnum.args = symnum.args) +  # Use approximate p-value
  stat_compare_means(label.y = 7)

print(p6)

####### Post glucose salivary pH #######
kruskal.test(richness$Observed ~ sample_data(ps)$post.gluocose_di)
kruskal.test(richness$Shannon ~ sample_data(ps)$post.gluocose_di)
comps <- make_pairs(sample_data(ps)$post.gluocose_di)

p7 <- plot_richness(ps, x = "post.gluocose_di", measures = c("Observed", "Shannon"), 
                    color = "post.gluocose_di") +
  geom_boxplot(aes(fill = post.gluocose_di), color = "black") + 
  scale_fill_manual(values = c("High pH" = "#0072B2", "Low pH" = "#D55E00")) + 
  labs(x = "Saliva pH post-glucose challenge") +
  theme_minimal() + 
  theme(strip.background = element_blank(),
        legend.position = "none") +
  geom_jitter(aes(color = post.gluocose_di), position = position_jitter(width = 0.2)) +
  stat_compare_means(comparisons = comps, method = "wilcox.test", 
                     label = "p.signif", exact = FALSE, label.y = c(5, 7, 6), 
                     symnum.args = symnum.args) +  # Use approximate p-value
  stat_compare_means(label.y = 7)

print(p7)


####### Professionally applied Fluoride #######

kruskal.test(richness$Observed ~ sample_data(ps)$fluoride_prof_status)
kruskal.test(richness$Shannon ~ sample_data(ps)$fluoride_prof_status)
comps <- make_pairs(sample_data(ps)$fluoride_prof_status)

p8 <- plot_richness(ps, x = "fluoride_prof_status", measures = c("Observed", "Shannon"), 
                    color = "fluoride_prof_status") +
  geom_boxplot(aes(fill = fluoride_prof_status), color = "black") + 
  scale_fill_manual(values = c("Fluoride applied" = "#0072B2", "No fluoride" = "#D55E00")) + 
  labs(x = "Professionally applied fluoride") +
  theme_minimal() + 
  theme(strip.background = element_blank(),
        legend.position = "none") +
  geom_jitter(aes(color = fluoride_prof_status), position = position_jitter(width = 0.2)) +
  stat_compare_means(comparisons = comps, method = "wilcox.test", 
                     label = "p.signif", exact = FALSE, label.y = c(5, 7, 6), 
                     symnum.args = symnum.args) +  # Use approximate p-value
  stat_compare_means(label.y = 7)

print(p8)

library(ggplot2)
library(cowplot)


combined_alpha_s_plot <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 3, labels = c("A", "B", "C", "D", "E","F", "G"),
                           label_size = 20)
combined_alpha_s_plot <- ggdraw() + draw_plot(combined_alpha_s_plot) + 
  theme(plot.background = element_rect(fill = "white", color = "white"))
print(combined_alpha_s_plot)

ggsave("combined_alpha_s_plot.png", plot = combined_alpha_s_plot, width = 20, height = 13, dpi = 300)

# Save the plot as SVG
ggsave(filename = "alpha_diversity_sig.svg", plot = p_alpha, width = 20, height = 10)


########
##############################################
#############LINEAR REGRESSION MODELS#########
#Alpha diversity of the continous variable, 
#1. Energy consumption
#2. Water consumption
#3. Carbohydrate consumption
#4. Sugar consumption
#5. Caries score
#6. Salivary pH
#7. Salivary flow rate
#8. Professionally applied fluoride
##################################
set.seed(123)
otu_mf <- lm(Observed ~ mf + age + gender.factor.x, data = alpha_diversity)
summary(otu_mf)

sha_mf <- lm(Shannon ~ mf + age + gender.factor.x, data = alpha_diversity)
summary(sha_mf)

otu_energy <- lm(Observed ~ Energy.inc.fibre + age + gender.factor.x, data = alpha_diversity)
summary(otu_energy)
sha_energy <- lm(Shannon ~ Energy.inc.fibre + age + gender.factor.x, data = alpha_diversity)
summary(sha_energy)

otu_water <- lm(Observed ~ Water + age + gender.factor.x, data = alpha_diversity)
summary(otu_water)
sha_water <- lm(Shannon ~ Water + age + gender.factor.x, data = alpha_diversity)
summary(sha_water)


otu_carb <- lm(Observed~ Carbohydrate + age + gender.factor.x, data = alpha_diversity)
summary(otu_carb)
sha_carb <- lm(Shannon ~ Carbohydrate + age + gender.factor.x, data = alpha_diversity)
summary(sha_carb)

otu_sugar <- lm(Observed~ Sugars + age + gender.factor.x, data = alpha_diversity)
summary(otu_sugar)
sha_sugar <- lm(Shannon ~ Sugars + age + gender.factor.x, data = alpha_diversity)
summary(sha_sugar)

otu_flow <- lm(Observed ~ saliva_flow_rate + age + gender.factor.x, data = alpha_diversity)
summary(otu_flow)
sha_flow <- lm(Shannon ~ saliva_flow_rate + age + gender.factor.x, data = alpha_diversity)
summary(sha_flow)

otu_pg2 <- lm(Observed ~ post.glucose_2 + age + gender.factor.x, data = alpha_diversity)
summary(otu_pg2)
sha_pg2 <- lm(Shannon ~ post.glucose_2 + age + gender.factor.x, data = alpha_diversity)
summary(sha_pg2)

otu_pg2 <- lm(Observed ~ fluoride_prof_status+ age + gender.factor.x, data = alpha_diversity)
summary(otu_pg2)
sha_pg2 <- lm(Shannon ~ fluoride_prof_status + age + gender.factor.x, data = alpha_diversity)
summary(sha_pg2)

otu_all <- lm(Observed ~ mf + Energy.inc.fibre +  Carbohydrate + Sugars + 
                 saliva_flow_rate + post.glucose_2 + fluoride_prof_status +
                age + gender.factor.x, data = alpha_diversity)
summary(otu_all)

shannon_all <- lm(Shannon ~ mf + Energy.inc.fibre +  Carbohydrate + Sugars + 
                saliva_flow_rate + post.glucose_2 +
                age + gender.factor.x, data = alpha_diversity)
summary(shannon_all)

##########################################
ggplot(alpha_diversity, aes(x = mf, y = Shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(title = "Linear Regression", x = "mf", y = "Shannon")


##########################################
#############Beta diversity#########
###########################################
ps1<- physeq
ps1.rel <- microbiome::transform(ps1, "compositional")
otu <- abundances(ps1.rel)
meta <- meta(ps1.rel)

# Run PERMANOVA on aitchison
adonis2(t(otu) ~ Energy.inc.fibre, data = meta, permutations = 999, method = "euclidean")
adonis2(t(otu) ~ Water, data = meta, permutations = 999, method = "euclidean")
adonis2(t(otu) ~ Carbohydrate, data = meta, permutations = 999, method = "euclidean")
adonis2(t(otu) ~ Sugars, data = meta, permutations = 999, method = "euclidean")
adonis2(t(otu) ~ saliva_flow_rate, data = meta, permutations = 999, method = "euclidean")
adonis2(t(otu) ~ post.glucose_2, data = meta, permutations = 999, method = "euclidean")
adonis2(t(otu) ~ caries_experience, data = meta, permutations = 999, method = "euclidean")
adonis2(t(otu) ~ fluoride_prof_status, data = meta, permutations = 999, method = "euclidean")

adonis_eu_all<- adonis2(t(otu) ~ Energy.inc.fibre+  Water+ Carbohydrate+ Sugars+
                            saliva_flow_rate + post.glucose_2 + caries_experience+ fluoride_prof_status, 
                        data = meta, permutations = 999, method = "euclidean")


adonis2(t(otu) ~ Energy.inc.fibre, data = meta, permutations = 999, method = "bray")
adonis2(t(otu) ~ Water, data = meta, permutations = 999, method = "bray")
adonis2(t(otu) ~ Carbohydrate, data = meta, permutations = 999, method = "bray")
adonis2(t(otu) ~ Sugars, data = meta, permutations = 999, method = "bray")
adonis2(t(otu) ~ saliva_flow_rate, data = meta, permutations = 999, method = "bray")
adonis2(t(otu) ~ post.glucose_2, data = meta, permutations = 999, method = "bray")
adonis2(t(otu) ~ caries_experience, data = meta, permutations = 999, method = "bray")
adonis2(t(otu) ~ fluoride_prof_status, data = meta, permutations = 999, method = "bray")                        


adonis_bray_all<- adonis2(t(otu) ~  Energy.inc.fibre+  Water+ Carbohydrate_di + Sugars_di+
                            saliva_flow_rate + post.glucose_2+ caries_experience + fluoride_prof_status, 
                          data = meta, permutations = 999, method = "bray")



library(microbiome)
library(ggplot2)
library(dplyr)
library(vegan)
library(microViz)
##Fixing the tables
ps1<- physeq
tax_table(ps1) %>% head(3)
phyloseq_validate(ps1)
#tax_fix_interactive(ps1)

ps1 %>%
  tax_fix(
    min_length = 4,
    unknowns = c(""),
    sep = " ", anon_unique = TRUE,
    suffix_rank = "classified"
  )

ps3 <- tax_fix(ps1)  
ps3 %>% tax_fix(unknowns = c("uncultured"))
ps3 <- ps3 %>% tax_fix(unknowns = c("uncultured"))


ps_a <- ps3 %>%
  tax_filter(min_prevalence = 2.5 / 100, verbose = FALSE) %>%
  tax_agg(rank = "Genus") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "Energy_di") +
  theme_classic(12) +
  coord_fixed(0.7) +
  stat_ellipse(aes(color = Energy_di)) +
  scale_color_brewer(palette = "Set1")
ps_a

ps_b <- ps3 %>%
  tax_filter(min_prevalence = 2.5 / 100, verbose = FALSE) %>%
  tax_agg(rank = "Genus") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "Carbohydrate_di") +
  theme_classic(12) +
  coord_fixed(0.7) +
  stat_ellipse(aes(color = Carbohydrate_di)) +
  scale_color_brewer(palette = "Set1")
ps_b

ps_c <- ps3 %>%
  tax_filter(min_prevalence = 2.5 / 100, verbose = FALSE) %>%
  tax_agg(rank = "Genus") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "Sugars_di") +
  theme_classic(12) +
  coord_fixed(0.7) +
  stat_ellipse(aes(color = Sugars_di)) +
  scale_color_brewer(palette = "Set1")
ps_c

ps_d <- ps3 %>%
  tax_filter(min_prevalence = 2.5 / 100, verbose = FALSE) %>%
  tax_agg(rank = "Genus") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "post.gluocose_di") +
  theme_classic(12) +
  coord_fixed(0.7) +
  stat_ellipse(aes(color = post.gluocose_di)) +
  scale_color_brewer(palette = "Set1")
ps_d


ps_e <- ps3 %>%
  tax_filter(min_prevalence = 2.5 / 100, verbose = FALSE) %>%
  tax_agg(rank = "Genus") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "fluoride_prof_status") +
  theme_classic(12) +
  coord_fixed(0.7) +
  stat_ellipse(aes(color = fluoride_prof_status)) +
  scale_color_brewer(palette = "Set1")
ps_e

ps3 %>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Carbohydrate_di", size = 2) +  stat_ellipse(aes(color = Carbohydrate_di)) +
  scale_colour_brewer(palette = "Set1", aesthetics = c("fill", "colour")) +
  theme_bw() +
  ggside::geom_ysideboxplot(aes(fill = Carbohydrate_di, x = Carbohydrate_di), orientation = "x") +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::theme_ggside_void()



######Figure 2#########
Figure_2 <- plot_grid(p4, p5, perm_eu, perm_bray, ps_b, ps_c, ncol = 2,
                      labels = c("A", "B", "C", "D", "E", "F"), 
                      label_size = 20)
print(Figure_2)
Figure_2 <- ggdraw() + draw_plot(Figure_2) + 
  theme(plot.background = element_rect(fill = "white", color = "white"))
ggsave("Figure_2.svg", plot = Figure_2, device = "svg", width = 16, height = 12, units = "in")
ggsave("Figure_2.png", plot = Figure_2, width = 16, height = 12, dpi = 600)

################################################
##########Correlation analysis plots##############
###########################################
correlation_plot <- ps3 %>%
  tax_agg("Genus") %>%
  tax_sort(by = prev, at = "Genus") %>%
  cor_heatmap(
    seriation_method = "Identity",
    seriation_method_col = "OLO_ward",
    taxa = tax_top(ps3, 20, by = max, rank = "Genus"),
    vars = c( "mf", "Energy.inc.fibre", "Water", "Carbohydrate", "Sugars",
              "saliva_flow_rate",
              "post.glucose_2", "fluoride"),
    tax_anno = taxAnnotation(
      Prev. = anno_tax_prev(ylim = 0:1),
      CLR = anno_tax_box(trans = "clr", zero_replace = "halfmin")
    )
  )

# calculate distances
aitchison_dists <- ps3%>%
  tax_transform("identity", rank = "Genus") %>%
  dist_calc("aitchison")

# the more permutations you request, the longer it takes
# but also the more stable and precise your p-values become
aitchison_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = "Energy.inc.fibre"
  )
aitchison_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = "Water"
  )
aitchison_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = "Carbohydrate"
  )
aitchison_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = "Sugars"
  )
aitchison_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = "saliva_flow_rate"
  )
aitchison_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = "post.glucose_2"
  )

aitchison_dists %>%dist_permanova(
  seed = 1234, # for set.seed to ensure reproducibility of random process
  n_processes = 1, n_perms = 999, # you should use at least 999!
  variables = "fluoride"
)

perm2 <- aitchison_dists %>%
  dist_permanova(variables = c("mf", "Energy.inc.fibre", "Water", "Carbohydrate", "Sugars",
                                 "saliva_flow_rate", "post.glucose_2", "fluoride", "age"), seed = 321)


par_ord <- perm2 %>%
  ord_calc(
    constraints = c("Energy.inc.fibre", "Water", "Carbohydrate", 
                    "Sugars", "saliva_flow_rate", "fluoride",
                    "post.glucose_2", "mf"),
    condition = "age"
  ) %>%
  ord_plot(
    colour = "country_born.factor", 
    alpha = 0.35,
    auto_caption = 7,
    constraint_vec_length = 4,
    constraint_vec_style = vec_constraint(1.5, colour = "grey15"),
    constraint_lab_style = constraint_lab_style(
      max_angle = 90, 
      size = 3, 
      aspect_ratio = 0.8, 
      colour = "black"
    )
  ) +
  stat_ellipse(aes(colour = country_born.factor), linewidth = 0.2) +  # Ensure Sugars_di is a valid column in your data
  scale_color_brewer(palette = "Set1", guide = guide_legend(position = "right")) +
  coord_fixed(ratio = 0.8, clip = "off", xlim = c(-4, 4)) +
  theme(legend.position = "bottom", legend.background = element_rect(fill = "white"))





#############################################
##### Differential Abundance Analysis ########
#################################################
library(Maaslin2)
df_input_data = read.table(file = "table.txt", header = TRUE, sep = "\t",
                           row.names = 1,
                           stringsAsFactors = FALSE)
df_input_data[1:5, 1:5]
df_input_metadata = read.table(file= "metadata_phy.tsv", header = TRUE, sep = "\t",
                               row.names = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]

fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "neg-csffs",
                    fixed_effects = c("Energy.inc.fibre", "Water", "Carbohydrate", 
                                      "Sugars", "saliva_flow_rate", "fluoride",
                                      "post.glucose_2", "mf")
)
#2.
fit_data = Maaslin2(input_data = df_input_data, 
                    input_metadata = df_input_metadata,
                    analysis_method = "NEGBIN",
                    min_prevalence = 0.5,
                    min_abundance = 0.1,
                    normalization = "CSS",
                    transform = "NONE",
                    output = "neg-cat",
                    fixed_effects = c("Energy_di", "Water_di", "Carbohydrate_di", 
                                      "Sugars_di", "Saliva_flow_di", "fluoride_prof_status",
                                      "post.gluocose_di", "caries_experience")
)
#######################################################################################
#################Creating a plot from the results of diffrential abundance analyis#####
######################################################################################
library(ggplot2)
library(dplyr)

data <- read.csv("mas_da.csv") # I have saved the file, and has feature, coff

# Convert metadata to a factor for consistent ordering
data$metadata <- factor(data$metadata, levels = c("Diet", "Physiological", "Clinical"))

# Arrange data by metadata category and then by coefficient
data <- data %>%
  arrange(metadata, coef)

# Create a combined label for the y-axis (feature + value)
data$combined_label <- paste0(data$feature, " (", data$value, ")")

# Create the forest plot
da_plot <-ggplot(data, aes(x = coef, y = reorder(combined_label, interaction(metadata, coef)), color = metadata)) +
  geom_point(size = 3) + # Plot points
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgreen") + # Reference line at 0
  scale_color_manual(values = c("Diet" = "blue", "Physiological" = "orange", "Clinical" = "red")) + # Customize colors
  labs(
    x = "Coefficient",
    y = "Taxa (Metadata)",
    color = "Metadata Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10), # Adjust y-axis text size
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "right"
  )
###################################
#########Figure 1#################
###################################
# Ensure the ComplexHeatmap package is loaded
library(ComplexHeatmap)
library(grid)
library(cowplot)
# Convert the ComplexHeatmap object to a grob
corr_plot_grob <- grid.grabExpr(draw(correlation_plot))
# Use the grob in the plot_grid
col2 <- plot_grid(corr_plot_grob, da_plot, ncol = 1, labels = c("B", "C"), label_size = 20)
Figure_1 <- plot_grid(par_ord, col2, ncol = 2, labels = c("A", ""), label_size = 20)
print(Figure_1)
Figure_1 <- ggdraw() + draw_plot(Figure_1) + 
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("Figure_1.svg", plot = Figure_1, device = "svg", width = 15, height = 8, units = "in")
ggsave("Figure_1.png", plot = Figure_1, width = 15, height = 8, dpi = 600)
################################

################################### 
######### MEDIATION ANALYSIS#######
###################################
#data <-read_tsv("metadata_phy.tsv")
data <- data %>%
  mutate(
    caries_experience_binary = ifelse(caries_experience == "Caries experience", 1, 0),
    country_born_binary = ifelse(country_born.factor == "Australia", 1, 0), 
    education_binary = ifelse(education_status == "Tertiary education", 1, 0),
    dental_visit_binary = ifelse(last_visit_dent_status == "Less than 12 months", 1, 0),
    sugar_binary = ifelse(Sugars_di == "High sugar", 1, 0),
    saliva_binary = ifelse(saliva_flow_status == "Normal saliva flow rate", 1, 0)
  )

################################################
############# PROPENSITY SCORE AND MATCHING#####
################################################
library(MatchIt)
library(cobalt)

m.out <- matchit(sugar_binary ~ I(age^2) + 
                   dental_visit_binary + fluoride + 
                   saliva_binary+ country_born_binary +  education_binary,
                 data=data, distance="logit", estimand="ATT", method='nearest', 
                 replace=TRUE, ratio=16)
summary(m.out, standardize=TRUE)

love.plot(bal.tab(m.out, m.threshold=0.25), stat = "mean.diffs", grid=TRUE, stars="raw", abs = F)


m.out$match.matrix

bal.tab(m.out, m.threshold = 0.3, un = TRUE)
#bal.tab(m.out, v.threshold = 2)

mdata <- match.data(m.out)
head(mdata)


##Plots
plot(m.out)
plot(m.out, type = 'jitter')
plot(m.out, type = 'hist')
love.plot(bal.tab(m.out, m.threshold=0.25), stat = "mean.diffs", grid=TRUE, stars="raw", abs = F)

write.table(mdata, "mdata.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#####################################################################################
####PHYLOSEQ OBJECT######################################
#########################################################
library(phyloseq)
library(qiime2R)
physeq<-qza_to_phyloseq(
  features="table.qza",
  tree="rooted-tree.qza",
  "taxonomy.qza", 
  metadata= "mdata.tsv"
)
physeq
ps <- physeq

# Assuming you have a phyloseq object `ps`
# Step 1: Create a relative abundance table
# Extract the OTU/ASV table and transform it to relative abundance
otu_table_rel <- transform_sample_counts(ps, function(x) x / sum(x))
# Convert to a data frame
abundance_table <- as.data.frame(t(otu_table(otu_table_rel)))

# Step 2: Ensure row names (taxa names) are column names
# Optionally remove metadata if present in phyloseq sample_data
rownames(abundance_table) <- NULL

# Step 3: Define filtering functions
abundance_filter <- function(df, min_abundance = 0.01) {
  # Ensure the data is numeric
  df <- df %>% select_if(is.numeric)
  # Retain taxa with mean relative abundance above the specified threshold
  df_filtered <- df[, colMeans(df) > min_abundance]
  return(df_filtered)
}

prevalence_filter <- function(df, min_prevalence = 0.05) {
  # Ensure the data is numeric
  df <- df %>% select_if(is.numeric)
  # Retain taxa present in at least a specified percentage of samples
  prevalence_threshold <- min_prevalence * nrow(df)
  df_filtered <- df[, colSums(df > 0) >= prevalence_threshold]
  return(df_filtered)
}

# Step 4: Apply the filters
# Adjust these thresholds as needed
min_abundance_threshold <- 0.01 
min_prevalence_threshold <- 0.05  
#####################
# Filter by abundance
abundance_filtered_table <- abundance_filter(abundance_table, min_abundance = min_abundance_threshold)

# Filter by prevalence
final_filtered_table <- prevalence_filter(abundance_filtered_table, min_prevalence = min_prevalence_threshold)

# Step 5: Display the resulting filtered table
print(final_filtered_table)
#######################################################
treatment <- mdata$sugar_binary
outcome <- mdata$post.glucose_2

library(SparseMCMM)
set.seed(1234)
res=SparseMCMM(treatment, final_filtered_table, outcome, n.split=100, num.per=NULL)
View(res)
res

reg <- lm(outcome ~ treatment)

res_summary <- as.data.frame(res$`Esitmated Causal Effects`)
res_component_me <- as.data.frame(res$`Compontent-wise ME`)


# Save the main results to CSV files
write.csv(res_summary, "SparseMCMM_Causal_Effects.csv", row.names = FALSE)
write.csv(res_component_me, "SparseMCMM_Component_ME.csv", row.names = FALSE)

# If you need to save the entire `res` object as an R data file for later use
save(res, file = "SparseMCMM_results.RData")



# Load the dataset
component_ME <- read.csv("component_ME.csv")  # Replace with the path to your file

# Add a unique identifier for each row to avoid merging duplicates
component_ME$unique_id <- paste0(component_ME$feature, "_", seq_len(nrow(component_ME)))

# Create the forest plot
forest_plot <- ggplot(component_ME, aes(x = mean, y = reorder(unique_id, mean))) +
  geom_point(size = 3, color = "black") + # Points for mean values
  geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2, color = "black") + # Error bars for CI
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") + # Reference line at 0
  theme_minimal() +
  labs(
    x = "Effect Size (Mean with 95% CI)",
    y = "Taxa",
  ) +
  scale_y_discrete(labels = component_ME$feature) + # Display only taxa names on the y-axis
  theme(
    axis.text.y = element_text(size = 10), # Adjust y-axis text size
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12)
  )
print(forest_plot)

######STACKED BAR PLOT######
# Load the necessary library
library(ggplot2)

data <- data.frame(
  Effect = c("Direct Effect", "Indirect Effect"),
  Value = c(-0.026, -0.407)
)

stacked_bar_plot <- ggplot(data, aes(x = "Effects", y = Value, fill = Effect)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Direct Effect" = "#1E88E5", "Indirect Effect" = "#FFC107")) + # Custom colors
  labs(
    x = "Effects Sugar consumption on saliva pH",
    y = "Effect Size",
    title = "Stacked Vertical Bar Plot of Direct and Indirect Effects",
    fill = "Effect Type" # Legend title
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), # Remove x-axis labels for a cleaner look
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right" # Position of the legend
  )
# Print the plot
print(stacked_bar_plot)


##############FIGURE 3#########
# Load necessary libraries
library(cowplot)
library(ggplot2)
library(grid)
library(png)

# Load the PNG image
img <- readPNG("mediation_analysis.png")  # Replace with your PNG file path
img_grob <- rasterGrob(img, interpolate = TRUE)  # Convert to grob for use in cowplot

# Arrange the plots into a grid
column_1 <- plot_grid(
  ggdraw() + draw_grob(img_grob, scale = 0.8),  # Insert PNG at the top
  stacked_bar_plot,                             # Stacked bar plot below PNG
  ncol = 1, labels = c("A", "B"), label_size = 20,
  rel_heights = c(2, 1)  # Adjust the size: PNG (3 parts), Stacked Bar (1 part)
)

column_2 <- forest_plot  # Forest plot in the second column

# Combine the two columns
Figure_3 <- plot_grid(
  column_1, column_2,
  ncol = 2,
  labels = c("", "C"),  # Add labels
  label_size = 20
)
print(Figure_3)

Figure_3 <- ggdraw() + draw_plot(Figure_3) + 
  theme(plot.background = element_rect(fill = "white", color = "white"))
ggsave("Figure_3.svg", plot = Figure_3, device = "svg", width = 12, height = 10, units = "in")
ggsave("Figure_3.png", plot = Figure_3, width = 12, height = 10, dpi = 600)
