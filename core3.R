setwd("~/Box Sync/sro-analysis-core/")
library(phyloseq)
library(microbiome)
otu_csv <- read.csv("abundance-table-final.csv", header=TRUE)
rownames(otu_csv) <- otu_csv[,1]
otu_csv[,1] <- NULL
otu <- otu_table(as.matrix(otu_csv), taxa_are_rows = FALSE)

taxa_csv <- read.csv("tax-final.csv", header=TRUE)
rownames(taxa_csv) <- taxa_csv[,1]
taxa_csv[,1] <- NULL

tax <- tax_table(as.matrix(taxa_csv))
meta <- read.csv("meta.txt", header=TRUE, sep="\t")
rownames(meta) <- meta$X
meta$X <- NULL
meta_phy <- sample_data(meta)

phy <- phyloseq(tax, otu, meta_phy)
taxa_names(phy) <- paste0("Seq", seq(ntaxa(phy)))

pseq.rel <- microbiome::transform(phy, "compositional")

otu_tab <- otu_table(pseq.rel)

hist(taxa_sums(pseq.rel))
sample_sums(pseq.rel)

head(prevalence(pseq.rel, detection = 1/100, sort = TRUE))
head(prevalence(pseq.rel, detection = 1/100, sort = TRUE, count = TRUE))

core.taxa.standard <- core_members(pseq.rel, detection = 0, prevalence = 50/100)
core.taxa.standard



pseq.core <- core(pseq.rel, detection = 0.0, prevalence = .05) 

pseq.core

core.taxa <- taxa(pseq.core)
                             
core.abundance <- sample_sums(core(pseq.rel, detection = 2/25, prevalence = .05))

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))


plot_core(pseq.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")

data(dietswap)
p <- plot_core(transform(dietswap, "compositional"),
  prevalences=seq(0.1, 1, .1), detections=seq(0.01, 1, length = 10))det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
plot_core(pseq.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")
print(p)

library(RColorBrewer)
library(reshape)

prevalences <- seq(.05, 1, .05)

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

#Added pseq.rel, I thin... must be checked if it was in the the rednred version,; where it is initialized
#pseq.rel<- microbiome::transform(pseq, 'compositional')
#min-prevalence gets the 100th highest prevalence
p <- plot_core(pseq.rel,
               plot.type = "heatmap", 
               colours = gray,
               prevalences = prevalences, 
               detections = detections, 
               min.prevalence = prevalence(pseq.rel, sort = TRUE)[100]) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  
  #Adjusts axis text size and legend bar height
  theme(axis.text.y= element_text(size=2, face="italic"),
        axis.text.x.bottom=element_text(size=2),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)

# Core with absolute counts and horizontal view:
# and minimum population prevalence (given as percentage)
detections <- seq(from = 50, to = round(max(abundances(pseq))/10, -1), by = 100)

p <- plot_core(pseq.rel, plot.type = "heatmap",
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .2, horizontal = TRUE) +
  theme(axis.text.x= element_text(size=1, face="italic", hjust=1),
        axis.text.y= element_text(size=1),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)

ps_bray <- ordinate(pseq.rel, "NMDS", "bray")
plot_ordination(pseq.rel, ps_bray, type="samples", color="crop") + geom_point(size = 3) 
plot_bar(pseq.core, fill="Phylum")

det <- c(0, 1, 2, 5, 20, 25)/100
prevalences <- seq(.01, .03, length=6)

plot_core(pseq.rel, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()

plot_core(pseq.rel, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()

