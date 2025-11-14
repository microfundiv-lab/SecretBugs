# load libraries
library(reshape2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(lme4)
library(vegan)
library(ggpubr)
library(ggrastr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")
top1.pos = scan("fastspar/top-pos-1_cag-170.txt", what="")
genome.metadata  = read.delim("metadata/genomes-nr_metadata.tsv")

# parse taxonomy
cag170.species = unique(genome.metadata[grep("g__CAG-170", genome.metadata$Lineage),"Species_rep"])
genome.selected = genome.metadata[genome.metadata$Species_rep %in% c(cag170.species, top1.pos), c("Genome","Species_rep", "Completeness")]
rownames(genome.selected) = genome.selected$Genome
genome.selected$Classification = NA
genome.selected$Classification = ifelse(genome.selected$Species_rep %in% cag170.species, "CAG-170", genome.selected$Classification)
genome.selected$Classification = ifelse(genome.selected$Species_rep %in% top1.pos, "Co-colonizers", genome.selected$Classification)

# load kegg presence/absence data
keggor.df = read.delim("genofan/cag170_ecology/kegg_orthologs.tsv", sep = "", header = FALSE)
keggor.df$V1 = gsub("_\\d+$", "", keggor.df$V1)
keggor.matrix = data.frame(acast(keggor.df, V2 ~ V1, length))

# filter dataset
selected.genomes = intersect(colnames(keggor.matrix), rownames(genome.selected))
genome.oscillo   = genome.selected[selected.genomes, , drop = FALSE]
keggor.matrix    = keggor.matrix[, selected.genomes, drop = FALSE]
metH = data.frame(t(keggor.matrix["K00548",]))
metH.fi = merge(metH, genome.oscillo, by="row.names")

# Chi-Sq test
fisher.df = table(metH.fi[,c("K00548", "Classification")])
rownames(fisher.df) = c("Absent", "Present")
res = fisher.test(chisq.df)
