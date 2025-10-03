#load libraries
library(ggplot2)

# load input
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/")
metadata.genomes = read.delim("metadata/genomes_uhgg-v1.2.tsv")
metadata.genomes = as.data.frame(table(metadata.genomes$Genome_type))
metadata.genomes$Prop = metadata.genomes$Freq/sum(metadata.genomes$Freq)*100
metadata.genomes$Data = "Genomes"
metadata.species = read.delim("metadata/species_uhgg_v1.2.tsv")
metadata.species = as.data.frame(table(metadata.species$Status))
metadata.species$Prop = metadata.species$Freq/sum(metadata.species$Freq)*100
metadata.species$Data = "Species"
metadata.df = rbind(metadata.genomes, metadata.species)

# plot donuts
donut.plot <- ggplot(metadata.df, aes(x = Data, y = Prop, fill = Var1)) +
  geom_col(width = 0.65, colour = "grey") +
  scale_x_discrete(limits = c(" ", " ", "Species","Genomes")) +
  scale_fill_manual(values = c("#b4c8ff", "#d9ffc1", "#b4c8ff", "#d9ffc1")) +
  coord_polar(theta = "y", start = 0, direction = 1) +  # still polar
  ylim(0, sum(metadata.df$Prop) * (240 / 360)) +         # restrict to 270°
  theme_void() +
  theme(legend.position = "none")