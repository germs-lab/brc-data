#Inter-BRC Workshop Sept 9, 2024
#sequences went through DADA2 pipeline to get them into an OTU table

#helpful R basics:  https://riffomonas.org/minimalR/
  
library(vegan)
library(phyloseq)
library(tidyverse)
library(ggplot2)

setwd("/Users/asuratt/Desktop/InterBRC_Workshop/brc-data-main")
ps <- readRDS('phyloseq-object.rds')
ps
otu_table(ps)[1:5, 1:5]  #Note that this code prints the first 5 rows and 5 columns of the matrix
tax_table(ps)[1:5]
sample_data(ps)[1:15] #this section just helps preview the data

#now we want to discover and pull out a potential "core" set of taxa
taxa_sums(ps)
sort(taxa_sums(ps), decreasing = TRUE)
names(sort(taxa_sums(ps), decreasing = TRUE)[1:10])
core_taxa_names <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:10]) #this pulled out the taxa names for the top 10 most abundant taxa across all samples

ps_core = prune_taxa(core_taxa_names, ps) #makes new taxa table that only has the core taxa we care about, using names list we just made above.
taxa_sums(ps_core) #sum of all observations of each taxa
sample_sums(ps_core) #sum of all taxa observe in each sample

#pull out samples where the core taxa are present in abundance of at least 100 ASVs
ps_core_alt <- filter_taxa(ps, function(x) mean(x) > 100, TRUE)

prevalance = apply(X = otu_table(ps_core),
                   MARGIN = 2,
                   FUN = function(x){sum(x > 0)})
tot_samples <- dim(sample_data(ps_core))[1] #this is the total number of samples
core_data = data.frame(Prevalence = prevalance, Prop_Prev = prevalance/tot_samples,
                       TotalAbundance = taxa_sums(ps_core),
                       tax_table(ps_core))

#You can explore distributions
library(ggplot2)
ggplot(data = core_data, aes(x = Prop_Prev)) +
  geom_histogram(binwidth = .1, fill = "blue", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Prevalence",
       x = "Prevalence",
       y = "Frequency")

ggplot(data = core_data, aes(x = TotalAbundance, y = Prop_Prev)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Scatter Plot of Prevalence vs Total Abundance",
       x = "Prevalence",
       y = "Total Abundance")

#sorting taxa by prevalence and/or abundance - to further subset our top 10 taxa if we want
threshold_prev = 0.05
threshold_abund = 1000
prev_cutoff <- core_data[core_data$Prop_Prev > threshold_prev,] # Filter first by only the prevalence threshold
prev_cutoff
abund_cutoff <- prev_cutoff[prev_cutoff$TotalAbundance > threshold_abund,] #further filter by abundance threshold
taxa_to_filter <- rownames(abund_cutoff)
length(taxa_to_filter)
mycore = prune_taxa(taxa_to_filter, ps) #filters original ASV table to save only the "core" taxa that meet a certain threshold of both prevalence AND abundance
#before we had 10 core taxa, but after the thresholds we cut it down to 6 taxa located within 25 samples

#Defining the core a different way, maybe with more taxa instead of just the top 10:
other_taxa = names(sort(taxa_sums(ps), decreasing = TRUE)[1:100])
other_ps = prune_taxa(other_taxa, ps)
prev = apply(X = otu_table(other_ps), #pull out taxa that have prevalence of more than 0
             MARGIN = 2,
             FUN = function(x){sum(x > 0)})
dim(sample_data(other_ps))[1] #gets number of taxa that met above prevalence threshold (of >0) for the taxa of interest
#next make a new dataframe that contains prevalence, proportion prevalence (prev/number of samples), total abundance - and adds it to the "other_ps" phyloseq object 
data_table = data.frame(Prevalence = prev, Prop_Prev = prev/dim(sample_data(other_ps))[1],
                        TotalAbundance = taxa_sums(other_ps),
                        tax_table(other_ps)) #taxa table subsetted "ps" by the top 100 taxa
threshold_prev = 0.05 #now we're once again subsetting the top 100 taxa by our thresholds of choice
threshold_abund = 1000
prev_cutoff <- core_data[core_data$Prop_Prev > threshold_prev,] # Pulls out taxa that meet the prevalence threshold
prev_cutoff
abund_cutoff <- prev_cutoff[prev_cutoff$TotalAbundance > threshold_abund,] #same as above, but for abundance instead
core_taxa <- rownames(abund_cutoff) #list of names of the taxa that meet the abundance cutoff
data_table$core_id <- ifelse(rownames(data_table) %in% core_taxa, "core", "other") #new column for taxa with "core" or "other" for graphing in next step, so we can color code

ggplot(data = data_table, aes(x = TotalAbundance, y = Prop_Prev, color=core_id)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Scatter Plot of Prevalence vs Total Abundance",
       x = "Total Abundance",
       y = "Prevalence")

#try making a taxa barplot of the core taxa
mycore
mycore_norm = transform_sample_counts(mycore, function(x) x / sum(x)) #transforming absolute counts to relative abundance
plot_bar(mycore_norm, fill="Genus")
f <- psmelt(mycore_norm) #converting phyloseq object to a dataframe that's easier to view/graph
ggplot(f, aes(x=Sample, y=Abundance, fill = Class))+geom_bar(stat="identity", position="stack")+facet_grid(~crop, scales="free_x")

