# load libraries
library(ggplot2)
library(ggridges)
library(ggpubr)
library(data.table)
library(tidyverse)
library(RColorBrewer)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/")
species.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
uncultured = species.metadata[which(species.metadata$Status == "Uncultured"),"Species_rep"]
metadata = read.delim("metadata/metagenomes_03-2024_samples.tsv")
rownames(metadata) = metadata$Sample

counts.filt = as.data.frame(fread("bwa/bwa_counts-filtered_samples.csv"))
rownames(counts.filt) = counts.filt$Genome
counts.filt = counts.filt[,-1]
source("../scripts/alex/metadata_disease-numbers.R")
counts.filt = counts.filt[,rownames(metadata.disease)]

# taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.df = unique(species.metadata[,c("Species_rep", "Lineage", "Status")])
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
tax.df$Taxon = NA
for (n in 1:nrow(tax.df)){
  for (r in ranks){
    if(!grepl("__$", tax.df[n,r])){
      tax.df[n,"Taxon"] = tax.df[n,r] }}}

# check most prevalent uncult
uncult.df = counts.filt[uncultured,]
uncult.prev = data.frame(sort(rowSums(uncult.df > 0), decreasing=TRUE))
colnames(uncult.prev) = "Prevalence"
uncult.prev$Prevalence = uncult.prev$Prevalence/ncol(uncult.df)*100
uncult.prev$Genus = tax.df[match(rownames(uncult.prev), tax.df$Species_rep),"Genus"]
uncult.prev$Family = tax.df[match(rownames(uncult.prev), tax.df$Species_rep),"Family"]
uncult.min1perc = length(which(uncult.prev$Prevalence > nrow(metadata.disease)*0.01))

# get top genera
nonsingle.genus = names(which(table(uncult.prev$Genus) > 1))
gen.prev = aggregate(Prevalence ~ Genus, data=uncult.prev, FUN=median)
gen.prev = gen.prev[order(gen.prev$Prevalence, decreasing=TRUE)[1:50],"Genus"]
uncult.prev = uncult.prev[which(uncult.prev$Genus %in% gen.prev & uncult.prev$Genus %in% nonsingle.genus),]
uncult.prev$Genus = gsub("g__", "", uncult.prev$Genus)
uncult.prev$Genus = gsub("_", " ", uncult.prev$Genus)
uncult.prev$Family = gsub("f__", "", uncult.prev$Family)
uncult.prev$Family = gsub("_", " ", uncult.prev$Family)
uncult.prev$N=table(uncult.prev$Genus)[uncult.prev$Genus]
uncult.prev$Label = paste(uncult.prev$Genus, " (n = ", uncult.prev$N, ")", sep="")

# plot top genera
gen.plot = ggplot(uncult.prev, aes(x=reorder(Label, -Prevalence), y=Prevalence, fill=Family)) +
  geom_boxplot(alpha=0.8, width=0.5) +
  geom_point(alpha=0.4, size=1, position = position_jitterdodge(jitter.width = 0.25)) +
  theme_classic() +
  xlab("") +
  ylab("Prevalence (%)") +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(size=11, angle=75, vjust=1, hjust=1)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title = element_text(size=14)) +
  scale_fill_manual(values=rev(brewer.pal(11, "Set3")))

# combine plots
source("../scripts/alex/mapping_rates.R")
ggarrange(gen.plot, uncult.plot, nrow=2, labels=c("a", "b"), font.label=list(size=16))
ggsave(file="../../figures/top-prev_uncult.pdf", width=10, height=10)
