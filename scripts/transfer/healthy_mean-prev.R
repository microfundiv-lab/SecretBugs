# load libraries
library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(robustbase)

# input and output
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")
metadata = read.delim("metadata/metagenomes_03-2024_samples.tsv", check.names = TRUE)
rownames(metadata) = metadata$Sample
abund.data = fread("bwa/bwa_counts-filtered_batch-corr_samples.csv")
abund.data.df = data.frame(abund.data[,-1], check.names = FALSE)
rownames(abund.data.df) = abund.data$Genome
abund.data.df = abund.data.df[,metadata$Sample]
abund.data.df = abund.data.df[,which(colSums(abund.data.df) > 0)]
metadata = metadata[colnames(abund.data.df),]

# select only healthy individuals
metadata = metadata[which(metadata$Disease.name == "Healthy" & metadata$Age.group %in% c("Elderly", "Adult", "Teenager")),]
metadata$Continent = gsub(" ", "_", metadata$Continent)
abund.data.df = abund.data.df[,rownames(metadata)]

# parse taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.df = read.delim("metadata/species_uhgg_v1.2.tsv", stringsAsFactors = FALSE)
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
rownames(tax.df) = tax.df$Species_rep
tax.df$Taxon = NA
for (n in 1:nrow(tax.df)){
  for (r in ranks){
    if(!grepl("__$", tax.df[n,r])){
      tax.df[n,"Taxon"] = tax.df[n,r] }}}

# calculate species prevalence per continent
continents = sort(unique(metadata$Continent))
prev.df = data.frame(matrix(NA, nrow=nrow(abund.data.df), ncol=length(continents)))
rownames(prev.df) = rownames(abund.data.df)
colnames(prev.df) = continents

for (c in continents) {
  samples = metadata[which(metadata$Continent == c),"Sample"]
  prev.df[,c] = rowSums(abund.data.df[,samples] > 0)/length(samples)*100
}
prev.df$mean_prevalence = rowMeans(prev.df)

# add species metadata
prev.fi = merge(prev.df, tax.df, by="row.names")
prev.fi = prev.fi[order(prev.fi$mean_prevalence, decreasing=TRUE),]

# save table
write.table(prev.fi, file = "healthy_analysis/prevalence_mean.tsv", sep="\t", row.names=FALSE, quote=FALSE)