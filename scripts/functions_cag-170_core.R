# load libraries
library(reshape2)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(matrixStats)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")
species.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
species.cag170  = species.metadata[grepl("g__CAG-170", species.metadata$Lineage), ]
genome.metadata  = read.delim("metadata/genomes-nr_metadata.tsv")
rownames(species.cag170) = species.cag170$Genome
genome.cag170 = genome.metadata[genome.metadata$Species_rep %in% species.cag170$Species_rep, c("Genome","Species_rep", "Completeness")]
genome.cag170 = genome.cag170[which(genome.cag170$Completeness > 90),]
rownames(genome.cag170) = genome.cag170$Genome

# load kegg presence/absence data
keggor.df = read.delim("oscillo/genofan/kegg_modules.tsv", header = TRUE)
keggor.matrix = data.frame(acast(keggor.df, module_accession ~ genome, value.var="completeness"))

# load gutsmash data
gutsmash.df = read.delim("genofan/cag170_ecology/gutsmash_results.tsv", header = FALSE)
gutsmash.matrix = t(data.frame(acast(gutsmash.df, V3 ~ V1)))
gutsmash.matrix = gutsmash.matrix[intersect(rownames(genome.cag170), rownames(gutsmash.matrix)),]
gutsmash.matrix = gutsmash.matrix[,names(which(colSums(gutsmash.matrix) > 0))]

# load antismash data
antismash.df = read.delim("genofan/cag170_ecology/antismash_results.tsv", header = FALSE)
antismash.matrix = t(data.frame(acast(antismash.df, V2 ~ V1)))
antismash.matrix = antismash.matrix[intersect(rownames(genome.cag170), rownames(antismash.matrix)),]
antismash.matrix = antismash.matrix[,names(which(colSums(antismash.matrix) > 0))]

# combine smash results
comb.smash = merge(gutsmash.matrix, antismash.matrix, by="row.names")
rownames(comb.smash) = comb.smash$Row.names
comb.smash = comb.smash[,-1]
comb.smash = comb.smash[,names(which(colSums(comb.smash > 0)/nrow(genome.cag170)*100 > 1))]

# define metab colours
metab.df = data.frame(colnames(comb.smash))
metab.df$Type = ifelse(metab.df[,1] %in% colnames(gutsmash.matrix), "Primary metabolism", "Secondary metabolism")
rownames(metab.df) = metab.df[,1]
metab.df = metab.df[,-1, drop=FALSE]

metab.colours = c("steelblue", "#D7C7F9")
names(metab.colours) = c("Primary metabolism", "Secondary metabolism")
metab.class = list(Type=metab.colours)

# rename columns
renamed.cols = colnames(comb.smash)
renamed.cols = gsub("_", " ", renamed.cols)
renamed.cols = gsub("betalactone", "Beta-lactone", renamed.cols)
renamed.cols = gsub("Others HGD unassigned", "2-hydroxyglutaryl-CoA dehydratase (other)", renamed.cols)
renamed.cols = gsub("cyclic-lactone-autoinducer", "Cyclic lactone autoinducer", renamed.cols)
renamed.cols = gsub("OD", "Oxidative decarboxylation", renamed.cols)
renamed.cols = gsub("histidine2", "Histidine to ", renamed.cols)
renamed.cols = gsub("Arginine2", "Arginine to", renamed.cols)
renamed.cols = gsub("porA", "R-pyruvate to R-acetate (porA)", renamed.cols)
renamed.cols = gsub("unknown", "(unknown)", renamed.cols)

# plot
comb.smash[comb.smash > 0] = 1
pheatmap(t(comb.smash), labels_row = renamed.cols, show_colnames = FALSE,
         annotation_row = metab.df, annotation_colors = metab.class, color=c("#fdfde4", "#f3756b"),
         filename = "figures/cag170_core_heatmap.pdf", width=9, height=5)

# load cazy data
cazy.df = read.delim("genofan/cag170_ecology/cazy_results.tsv", header = FALSE, sep=" ")
cazy.df = cazy.df[which(grepl("GUT_GENOME", cazy.df$V1)),]
cazy.df$V1 = sub("_[^_]+$", "", cazy.df$V1)
cazy.matrix = t(data.frame(acast(cazy.df, V2 ~ V1)))
cazy.matrix = cazy.matrix[intersect(rownames(genome.cag170), rownames(cazy.matrix)),]
cazy.sums = data.frame(colSums(cazy.matrix > 0)/nrow(genome.cag170)*100)
colnames(cazy.sums) = "CAG170_prop"
cazy.classes = read.csv("genofan/cag170_ecology/cazy_classes.csv", row.names=1)
cazy.fi = merge(cazy.sums, cazy.classes, by="row.names")
cazy.fi = cazy.fi[order(cazy.fi$CAG170_prop, decreasing=TRUE),]

# plot
cazy.plot = ggplot(cazy.fi, aes(x=reorder(Row.names, CAG170_prop), y=CAG170_prop, fill=factor(Broad_Functional_Group))) +
  geom_bar(stat="identity", alpha=0.75) +
  theme_classic() +
  coord_flip() +
  xlab("CAZyme") +
  ylab("Prevalence in CAG-170 genomes (%)") +
  scale_fill_manual(values=brewer.pal(10, "Set3"), name="Functional class") +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size=12)) +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))
ggsave(filename = "figures/cag170_core_cazy.pdf", height=9, width=6)
