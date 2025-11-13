# load libraries
library(reshape2)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(matrixStats)
library(vegan)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")

species.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
species.oscillo  = species.metadata[grepl("g__CAG-170", species.metadata$Lineage), ]
genome.metadata  = read.delim("metadata/genomes-nr_metadata.tsv")

# parse taxonomy
ranks = c("Domain","Phylum","Class","Order","Family","Genus","Species")
species.oscillo  = separate(species.oscillo, col = Lineage, into = ranks, sep = ";")
rownames(species.oscillo) = species.oscillo$Genome

genome.oscillo = genome.metadata[genome.metadata$Species_rep %in% species.oscillo$Species_rep, c("Genome","Species_rep", "Completeness")]
genome.oscillo$Genus = species.oscillo[genome.oscillo$Species_rep, "Genus"]
genome.oscillo = genome.oscillo[genome.oscillo$Genus != "g__", ]
rownames(genome.oscillo) = genome.oscillo$Genome

# load kegg presence/absence data
keggor.df = read.delim("oscillo/kegg_orthologs.tsv", sep = "", header = FALSE)
keggor.df$V1 = gsub("_\\d+$", "", keggor.df$V1)
keggor.matrix = data.frame(acast(keggor.df, V2 ~ V1, length))

# filter dataset
selected.genomes = intersect(colnames(keggor.matrix), rownames(genome.oscillo))
genome.oscillo   = genome.oscillo[selected.genomes, , drop = FALSE]
keggor.matrix    = keggor.matrix[, selected.genomes, drop = FALSE]
keggor.matrix = keggor.matrix[rowSums(keggor.matrix > 0) > ncol(keggor.matrix) * 0.01, , drop = FALSE]

# calculate KO distances
ko.dist = as.matrix(vegdist(keggor.matrix, method="jaccard"))

# extract B12 KOs
b12.kos = scan("oscillo/B12_KOs.txt", what="")
b12.dist = ko.dist[b12.kos,setdiff(colnames(ko.dist), b12.kos)]

# histogram
ko.freq = data.frame(rowSums(keggor.matrix > 0))
b12.means = data.frame(colMeans(b12.dist))
b12.fi = merge(b12.means, ko.freq, by="row.names")
colnames(b12.fi) = c("KO", "Mean_dist", "Frequency")
b12.fi$Prop = b12.fi$Frequency/ncol(keggor.matrix)
b12.fi$Classification = ifelse(b12.fi$Prop > 0.9, "Core", "Accessory")
b12.hist = ggplot(b12.fi, aes(x=Mean_dist, fill=Classification)) +
  geom_histogram(colour="grey40", alpha=0.7, bins=50, linewidth=0.1) +
  theme_classic() +
  ylab("Number of KOs") +
  xlab("Average Jaccard distance to B12 biosynthesis KOs") +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))
ggsave(filename = "figures/cag170_b12_corr.pdf", width=9, height=3)
