# load libraries
library(ggplot2)
library(ggridges)
library(ggrastr)
library(ggpubr)
library(data.table)
library(scales)
library(tidyverse)
library(RColorBrewer)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/")
species.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
metadata = read.delim("metadata/metagenomes_03-2024_samples.tsv")
rownames(metadata) = metadata$Sample

# taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.df = separate(data = species.metadata, col = Lineage, sep = ";", into = ranks)
tax.df$Taxon = NA
for (n in 1:nrow(tax.df)){
  for (r in ranks){
    if(!grepl("__$", tax.df[n,r])){
      tax.df[n,"Taxon"] = tax.df[n,r] }}}

# load counts and convert to relab
counts.filt = as.data.frame(fread("bwa/bwa_counts-filtered_samples.csv"))
rownames(counts.filt) = counts.filt$Genome
counts.filt = counts.filt[,-1]
source("../scripts/alex/metadata_disease-numbers.R")
counts.filt = counts.filt[species.metadata$Genome,rownames(metadata.disease)]
counts.norm = counts.filt/species.metadata$Length
counts.relab = t(t(counts.norm)/colSums(counts.norm)*100)

# filter for cag170
cag170.species = tax.df[which(tax.df$Genus == "g__CAG-170"),"Species_rep"]
cag170.df = counts.relab[cag170.species,]
cag170.prev = data.frame(rowSums(cag170.df > 0)/nrow(metadata.disease)*100)
cag170.prev$Species = rownames(cag170.prev)
colnames(cag170.prev)[1] = "Prevalence"

cag170.abund = reshape2::melt(cag170.df)
cag170.abund = cag170.abund[which(cag170.abund$value != 0),]
colnames(cag170.abund) = c("Species", "Sample", "Abundance")

# plot prevalence and abundance
order.genomes = cag170.prev[order(cag170.prev$Prevalence, decreasing=FALSE),"Species"]
prev.plot = ggplot(cag170.prev, aes(x=Species, y=Prevalence)) +
  geom_bar(stat="identity", fill="palegreen4", alpha=0.8) +
  coord_flip() +
  theme_classic() +
  ylab("Prevalence (%)") +
  scale_x_discrete(limits=order.genomes) +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_blank())

abund.plot = ggplot(cag170.abund, aes(x=Species, y=Abundance)) +
  geom_point_rast(alpha=0.05, size=0.3, colour="grey", position = position_jitter(width=0.2)) + 
  geom_boxplot(outlier.shape = NA, fill="palegreen4", alpha=0.6) +
  coord_flip() +
  theme_classic() +
  scale_y_log10(labels = label_number()) +
  scale_x_discrete(limits=order.genomes) +
  ylab("Relative abundance (%)") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_blank())

ggarrange(prev.plot, abund.plot, widths=c(1,0.75))
