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

# Optional parallel (auto-fallback to lapply if not installed)
future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
glmer_apply = function(X, FUN) future.apply::future_lapply(X, FUN, future.seed = TRUE)
options(contrasts = c("contr.treatment", "contr.poly"))

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")
top1.pos = scan("fastspar/top-pos-1_cag-170.txt", what="")
genome.metadata  = read.delim("metadata/genomes-nr_metadata.tsv")
cobrapy.res = read.delim("genofan/cag170_ecology/cobrapy_M3.tsv", header=FALSE)
cobrapy.uptake = cobrapy.res[which(cobrapy.res$V4 == "Uptake"),]
cobrapy.secretion = cobrapy.res[which(cobrapy.res$V4 == "Secretion"),]

# parse taxonomy
cag170.species = unique(genome.metadata[grep("g__CAG-170", genome.metadata$Lineage),"Species_rep"])
genome.selected = genome.metadata[genome.metadata$Species_rep %in% c(cag170.species, top1.pos), c("Genome","Species_rep", "Completeness")]
rownames(genome.selected) = genome.selected$Genome
selected = intersect(genome.selected$Genome, unique(cobrapy.res$V1))
genome.selected = genome.selected[selected,]
genome.selected$Classification = NA
genome.selected$Classification = ifelse(genome.selected$Species_rep %in% cag170.species, "CAG-170", genome.selected$Classification)
genome.selected$Classification = ifelse(genome.selected$Species_rep %in% top1.pos, "Co-colonizers", genome.selected$Classification)

# compare metabolites
compare_metab = function(x, type) {
  cobrapy.matrix = data.frame(acast(x, V2 ~ V1, length))
  cobrapy.matrix = t(cobrapy.matrix[,selected])
  
  # calculate shannon diversity
  div.df = data.frame(specnumber(cobrapy.matrix))
  div.df$Classification = genome.selected[rownames(div.df),"Classification"]
  colnames(div.df) = c("Richness", "Classification")
  
  # plot boxplot
  box.plot = ggplot(div.df, aes(x=Classification, y=Richness, fill=Classification)) +
    geom_point_rast(colour="grey50", size=0.5, position = position_jitterdodge(jitter.width = 0.5)) +
    geom_boxplot(alpha=0.7, outlier.shape=NA) +
    ylab(paste0("Number of metabolites (", type, ")")) +
    scale_fill_manual(values = c("CAG-170" = 'palegreen3', "Co-colonizers" = "steelblue")) +
    guides(fill="none") +
    scale_x_discrete(limits=c("CAG-170", "Co-colonizers")) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size=14)) +
    theme(axis.text.x = element_text(size=12)) +
    theme(axis.text.y = element_text(size=12)) +
    geom_signif(
      comparisons = list(c("CAG-170", "Co-colonizers")),
      map_signif_level = FALSE
    )
}

metab.uptake = compare_metab(cobrapy.uptake, "uptake")
metab.secretion = compare_metab(cobrapy.secretion, "secretion")

# combine plots
ggarrange(metab.secretion, metab.uptake, nrow=2, font.label = list(size = 18))
ggsave(filename = "figures/cag170_metabolites_vs_coccur.pdf", width=3, height=7)
