# load necessary libraries
library(tidyverse)
library(data.table)
library(ggpubr)
library(RColorBrewer)

# load input file that has results of overlap Maaslin/Aldex
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/aldex2")
all_overlap_species = as.data.frame(fread("combined_across_diseases_meta.csv"))
tax.data = read.delim("../metadata/species_uhgg_v1.2.tsv")
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.data = separate(tax.data, col="Lineage", into=ranks, sep=";")

# add metadata
all_overlap_species$Classification = ifelse(all_overlap_species$FDR < 0.05, ifelse(all_overlap_species$meta_z > 0, "Disease", "Health"), "NS")
all_overlap_species$Genus = tax.data[match(all_overlap_species$species, tax.data$Genome),"Genus"]
all_overlap_species$Status = tax.data[match(all_overlap_species$species, tax.data$Genome),"Status"]

# add genus stats
genus.status = as.data.frame.matrix(table(all_overlap_species[,c("Genus", "Status")]))
genus.ds = aggregate(meta_effect ~ Genus, data=all_overlap_species, FUN=sum)
rownames(genus.ds) = genus.ds$Genus
genus.ds = genus.ds[,-1, drop=FALSE]
genus.class = as.data.frame.matrix(table(all_overlap_species[,c("Genus", "Classification")]))
genus.fi = cbind(genus.status, genus.ds, genus.class)
genus.fi$Total = rowSums(genus.fi[,c("Cultured", "Uncultured")])
genus.fi$Significant = genus.fi$Total-genus.fi$NS
genus.fi$Prop_uncult = genus.fi$Uncultured/genus.fi$Total
genus.fi$Norm_ES = genus.fi$meta_effect/genus.fi$Total
genus.fi$Prop_health = genus.fi$Health/genus.fi$Total
genus.fi$Prop_disease = genus.fi$Disease/genus.fi$Total

# calculate uncultured health/disease scores
genus.fi$UHS = genus.fi$Prop_uncult*genus.fi$Norm_ES*log10(genus.fi$Significant)
genus.fi = genus.fi[order(genus.fi$UHS, decreasing=FALSE),]
genus.fi = genus.fi[which(is.finite(genus.fi$UHS)),]

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

# plot genus stats
genus.health = genus.fi[which(genus.fi$UHS < 0),]
scatter = ggplot(genus.health, aes(x=Prop_uncult, y=Norm_ES, colour=N_class)) +
  geom_point(size=1.5) +
  theme_classic() +
  ylab("Normalized Effect Size") +
  xlab("% Uncultured species") +
  scale_colour_manual(values=c("steelblue", "darkgreen", "tomato", "purple"), name="Number of biomarkers") +
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

# arrange plots
ggarrange(scatter, lolli, align="h", widths=c(0.8,1), labels=c("a", "b"), font.label=list(size=18))
ggsave(file="../figures/cag170_validation.pdf", height=6, width=12)

# save table
write.csv(genus.fi, file="genus_stats.csv", row.names = F, quote=FALSE)