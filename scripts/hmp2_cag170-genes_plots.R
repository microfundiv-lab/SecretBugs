# load libraries
library(ggplot2)
library(RColorBrewer)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/timeseries/")
glm = read.delim("sign_genes_glmer.tsv")
eggnog = read.delim("sign_genes_functions.tsv", header=FALSE)[,c(1,21,22)]
freq = read.delim("sign_genes_tax.tsv", header=FALSE)
colnames(freq) = c("Gene", "N_proteins", "N_genomes", "N_genomes_hq", "N_species", "Species_rep", "Taxonomy")

# COG dictionary
COG_descriptions = c(
  B = "Chromatin Structure and dynamics",
  J = "Translation",
  L = "Replication, recombination & repair",
  K = "Transcription",
  O = "Molecular chaperones and related functions",
  M = "Cell wall structure & outer membrane",
  N = "Secretion, motility and chemotaxis",
  T = "Signal transduction",
  P = "Inorganic ion transport and metabolism",
  U = "Intracellular trafficking and secretion",
  C = "Energy production and conversion",
  G = "Carbohydrate metabolism and transport",
  E = "Amino acid metabolism and transport",
  F = "Nucleotide metabolism and transport",
  H = "Coenzyme metabolism",
  I = "Lipid metabolism",
  D = "Cell division and chromosome partitioning",
  R = "General functional prediction only",
  S = "No functional prediction",
  Q = "Secondary Structure",
  V = "Defense mechanisms",
  W = "Extracellular structures",
  Z = "Cytoskeleton"
)

# clean eggnNOG COG annotations
rownames(eggnog) = eggnog$V1
missing = glm$Gene[which(!glm$Gene %in% eggnog$V1)]
eggnog[missing,"V21"] = "S"
eggnog[missing,"V1"] = missing
eggnog$V21 = ifelse(eggnog$V21 == "", "S", eggnog$V21)
eggnog$V21 = ifelse(nchar(eggnog$V21) > 1, "R", eggnog$V21)
eggnog$V22 = ifelse(eggnog$V22 == "", NA, eggnog$V22)
eggnog$COG = COG_descriptions[eggnog$V21]
colnames(eggnog) = c("Gene", "COG", "eggNOG_annotation", "Description")

# combine data
merged = merge(glm, eggnog, by = "Gene", all = TRUE)
merged = merge(merged, freq, by = "Gene", all = TRUE)

# keep only top 10 COGs
keep.cogs = names(sort(table(eggnog$Description), decreasing=TRUE)[1:10])
merged.filt = merged
merged.filt$Description = ifelse(merged.filt$Description %in% keep.cogs, merged.filt$Description, "Other")
merged.filt$COG = ifelse(merged.filt$Description == "Other", "Other", merged.filt$COG)
cog.colors = brewer.pal(length(keep.cogs),"Set3")
cog.colors = c(cog.colors, "grey50")
names(cog.colors) = c(keep.cogs, "Other")

# plot gene freq vs effect size
gene.plot = ggplot(merged.filt, aes(x=log(N_genomes), y=Estimate, colour=Description)) +
  geom_point(size=2) +
  theme_classic() +
  ylab("log(Odds Ratio)") +
  xlab("log(Number of genomes)") +
  scale_colour_manual(values=cog.colors, name = "COG category") +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14)) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size=10))

# combine and save plots
source("../../scripts/alex/hmp2_cag170-dysb.R")
ggarrange(dysb.comb, gene.plot, nrow=2, labels=c("", "c"), font.label = list(size=18))
ggsave(file="../figures/cag170_dysbiosis.pdf", height=9, width=10)
write.table(merged, file="sign_genes_merged.tsv", sep="\t", quote=FALSE, row.names=FALSE)