# load necessary libraries
library(tidyverse)
library(data.table)
library(ggpubr)
library(RColorBrewer)

# load input file that has results of overlap Maaslin/Aldex
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/meta-analysis")
all_overlap_species = as.data.frame(fread("species_ds.csv"))

# add disease class
all_overlap_species$Classification = ifelse(all_overlap_species$DS > 0, "Disease", ifelse(all_overlap_species$DS < 0, "Health", "NS"))

# add genus stats
genus.status = as.data.frame.matrix(table(all_overlap_species[,c("Genus", "Status")]))
genus.ds = aggregate(DS ~ Genus, data=all_overlap_species, FUN=sum)
rownames(genus.ds) = genus.ds$Genus
genus.ds = genus.ds[,-1, drop=FALSE]
genus.class = as.data.frame.matrix(table(all_overlap_species[,c("Genus", "Classification")]))
genus.fi = cbind(genus.status, genus.ds, genus.class)
genus.fi$Total = rowSums(genus.fi[,c("Cultured", "Uncultured")])
genus.fi$Significant = genus.fi$Total-genus.fi$NS
genus.fi$Prop_uncult = genus.fi$Uncultured/genus.fi$Total
genus.fi$Norm_ES = genus.fi$DS/genus.fi$Total
genus.fi$Prop_health = genus.fi$Health/genus.fi$Total
genus.fi$Prop_disease = genus.fi$Disease/genus.fi$Total

# calculate uncultured health/disease scores
genus.fi$UHS = genus.fi$Prop_uncult*genus.fi$Norm_ES*log10(genus.fi$Significant)
genus.fi = genus.fi[order(genus.fi$UHS, decreasing=FALSE),]

# add species number class
genus.fi$N_class = ifelse(genus.fi$Significant > 10, ">10", ifelse(genus.fi$Significant > 5, "6-10", ifelse(genus.fi$Significant > 1, "2-5", "1")))
genus.fi$N_class = factor(genus.fi$N_class, levels=c("1", "2-5", "6-10", ">10"))
genus.fi = cbind(rownames(genus.fi), genus.fi)
colnames(genus.fi)[1] = "Genus"

# add family data
tax.data = read.delim("../metadata/species_uhgg_v1.2.tsv")
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.data = separate(tax.data, col="Lineage", into=ranks, sep=";")
genus.fi$Family = tax.data[match(genus.fi$Genus, tax.data$Genus),"Family"]
genus.fi$Genus = gsub("g__", "", genus.fi$Genus)
genus.fi$Family = gsub("f__", "", genus.fi$Family)

# plot genus-level stats
genus.health = genus.fi[which(genus.fi$UHS < 0),]
scatter = ggplot(genus.health, aes(x=Prop_uncult, y=Norm_ES, size=Significant, colour=Significant)) +
  geom_point() +
  theme_classic() +
  ylab("Normalized Effect Size") +
  xlab("% Uncultured species") +
  scale_colour_gradient(low="lightblue1", high="steelblue4") +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.position = "inside") +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12))

genus.health$UHS = abs(genus.health$UHS)
lolli = ggplot(genus.health[1:25,], aes(x=reorder(Genus, UHS), y=UHS, fill=Family)) +
  geom_linerange(aes(ymin=0, ymax=UHS), colour="grey", linewidth=1) +
  geom_point(size=4, shape=21, colour="darkgrey") +
  theme_classic() +
  coord_flip() +
  ylab("Uncultured health score") +
  scale_fill_manual(values=brewer.pal(length(unique(genus.fi[1:25,"Family"])), "Set3")) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.position = "inside") +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12))

genus.plot = ggarrange(scatter, lolli, align="h", widths=c(0.8,1), labels=c("A", "B"), font.label=list(size=18))

# plot cag-170 species
source("../../scripts/alex/cag170_prev-abund.R")
cag170.df = all_overlap_species[which(all_overlap_species$Genus == "g__CAG-170"),c("Genome", "CD", "ME", "Obesity", "UC", "DS", "Taxon")]
cag170.df$Label = paste(cag170.df$Genome, " (", gsub("s__CAG-170 ","", cag170.df$Taxon), ")", sep="")
rownames(cag170.df) = cag170.df$Genome
cag170.melt = reshape2::melt(cag170.df)
cag170.melt$value = round(cag170.melt$value, digits=2)
order.labels = cag170.df[order.genomes, "Label"]
  
cag170.heat = ggplot(cag170.melt, aes(y = Genome, x = variable, fill = value)) +
  geom_tile(color = "lightgrey", linewidth = 0.4) +
  geom_text(aes(label = ifelse(value == 0, "N.S.", value), fontface = ifelse(variable == "DS", "bold", "plain")), colour = "black") +
  scale_fill_gradient(low = "palegreen4", high = "white", name = "Effect Size") +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA), guide = "none") +  # black only for last col
  theme_minimal() +
  scale_y_discrete(limits = order.genomes, labels = order.labels) +
  scale_x_discrete(limits = c("CD", "ME", "Obesity", "UC", "DS"), 
                   labels = c("Crohn's disease", "Myalgic\nencephalomyelitis", "Obesity", 
                              "Ulcerative\ncolitis", "Normalized\nEffect Size")) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.position = "bottom")

cag170.plot = ggarrange(cag170.heat, prev.plot, abund.plot, ncol=3, widths=c(1,0.3,0.4), align="h")

# combine final plots
ggarrange(genus.plot, cag170.plot, ncol=1, labels=c("", "C"), heights=c(0.8,1), font.label=list(size=18))
ggsave(file="figures/cag170_stats.pdf", height=12, width=12)

# save table
write.csv(genus.fi, file="meta-analysis/genus_stats.csv", row.names = F, quote=FALSE)