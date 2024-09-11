
# Load Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggsci)
library(ggpubr)

# Load phyloseq object
ps <- readRDS('phyloseq-object.RDS')
ps

## TUTOTIAL #############################
otu_table(ps)[1:5, 1:5]  #Note that this code prints the first 5 rows and 5 columns of the matrix
tax_table(ps)[1:5]
sample_data(ps)[1:15]

# This is new, this changes the original phyloseq object to relative abundance, we should do this prior to the core selection
ps2  = transform_sample_counts(ps, function(x) x / sum(x) )
ps <- ps2 #I did this to cheat, it makes all my future code usable

# Explore data and subsetting core
taxa_sums(ps)
sort(taxa_sums(ps), decreasing = TRUE)
names(sort(taxa_sums(ps), decreasing = TRUE)[1:100])
core_taxa_names <- names(sort(taxa_sums(ps), decreasing = TRUE)[1:100])

ps_core <- prune_taxa(core_taxa_names, ps) # This core is defined as the top ten must abundant ASVs present in the dataset
taxa_sums(ps_core) #sum of all observations of each taxa
sample_sums(ps_core) #sum of all taxa observe in each sample

#ps_core_alt <- filter_taxa(ps, function(x) mean(x) > 100, TRUE) # This core is defined as ASVs with a greater mean abundance of 100

#You can add prevelance and abundance data to a dataframe
prevalance <- apply(X = otu_table(ps_core),
                   MARGIN = 2,
                   FUN = function(x){sum(x > 0)})
tot_samples <- dim(sample_data(ps_core))[1] #this is the total number of samples
core_data <- data.frame(Prevalence = prevalance, Prop_Prev = prevalance/tot_samples,
                       TotalAbundance = taxa_sums(ps_core),
                       tax_table(ps_core))

#You can explore distributions
library(ggplot2)
ps <- ps_core # sets the ps variable to represent core
ps_subset <- subset_samples(ps, row.names(sample_data(ps)) != "Control-1-12-25-root-16S_F_filt.fastq")
ps_subset <- subset_samples(ps_subset, row.names(sample_data(ps_subset)) != "Control-1-29-25-root-16S_F_filt.fastq")
psm <- psmelt(ps_subset)

f <- ddply(psm, .(crop, Class), summarise, MEAN=mean(Abundance))


ggplot(f, aes(x=crop, y=MEAN, fill=Class))+geom_bar(position="stack", stat="identity")+theme_bw()+xlab("Feedstock") + ylab("Average Relative abundance")

