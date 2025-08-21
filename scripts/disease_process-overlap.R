# load libraries
library(tidyverse)
library(data.table)
library(scales)
library(reshape2)

# input and output
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/meta-analysis")
genome.paths = read.delim("../metadata/genomes_uhgg-v1.2.tsv")[,c("Genome", "FTP_download")]
rownames(genome.paths) = genome.paths$Genome
species.metadata = read.delim("../metadata/species_uhgg_v1.2.tsv")
rownames(species.metadata) = species.metadata$Genome
species_list = scan("species_list.txt", what="")

# load results
input.files =  list.files(pattern = ".tsv")
disease.list = lapply(input.files, function(i) {
  disease = strsplit(i, split="\\.")[[1]][1]
  df = read.delim(i)
  colnames(df)[which(colnames(df) == paste0("Disease_name",disease,".Est"))] = "ALDEx2:Est"
  colnames(df)[which(colnames(df) == paste0("Disease_name",disease,".pval"))] = "ALDEx2:pval"
  return(df)
})

# combine all results into a dataframe
all.results = as.data.frame(rbindlist(disease.list))
all.results$mean = rowMeans(all.results[,c("coef", "ALDEx2:Est")])

# reshape df
reshape.df = as.data.frame(acast(feature ~ value, value.var = "mean", data=all.results))
reshape.df[is.na(reshape.df)] = 0
scores = as.vector(apply(reshape.df, 1, function(x) {sum(x)/length(which(x == 0))}))
reshape.df$DS = scores

# parse taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.df = unique(species.metadata[,c("Species_rep", "Lineage", "Status")])
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
tax.df$Taxon = NA
for (n in 1:nrow(tax.df)){
  for (r in ranks){
    if(!grepl("__$", tax.df[n,r])){
      tax.df[n,"Taxon"] = tax.df[n,r] }}}

# save candidate lists
candidate.list = genome.paths[rownames(reshape.df),]
write.table(candidate.list$Genome, file="candidate_species.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)
write.table(candidate.list$FTP_download, file="candidate_species_paths.txt", row.names=FALSE, quote=FALSE, col.names=FALSE)

# add missing species
missing.species = species_list[which(!species_list %in% rownames(reshape.df))]
reshape.df[missing.species,] = 0
ds.final = merge(reshape.df, tax.df, by="row.names")
colnames(ds.final)[1] = "Genome"
write.table(ds.final, file="species_ds.csv", row.names = F, quote=FALSE, sep=",")
