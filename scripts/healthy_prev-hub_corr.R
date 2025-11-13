# load libraries
library(data.table)
library(matrixStats)
library(ComplexUpset)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(RColorBrewer)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")

# parse taxonomy
species.data = read.delim("metadata/species_uhgg_v1.2.tsv")
rownames(species.data) = species.data$Species_rep
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
species.data = separate(data=species.data, col=Lineage, into=ranks, sep=";")

# process disease data
disease.df = read.csv("meta-analysis/species_ds.csv")
disease.df[is.na(disease.df)] = 0
disease.df = disease.df[disease.df$DS != 0,]
disease.df$Classification = ifelse(disease.df$DS < 0, "Health", "Disease")
disease.df = disease.df[,c("Species_rep", "Classification")]
colnames(disease.df)[2] = "Variable"

# process hub/core data
central = read.delim("healthy_analysis/central_mean.tsv")
central = central[,c("Species_rep", "mean_centrality", "Status")]
prevalence = read.delim("healthy_analysis/prevalence_mean.tsv")
prevalence = prevalence[,c("Species_rep", "mean_prevalence")]
hub.core = merge(central, prevalence, by="Species_rep")
rownames(hub.core) = hub.core$Species_rep
hub_thresh = quantile(hub.core$mean_centrality, 0.99, na.rm=TRUE)
core_thresh = quantile(hub.core$mean_prevalence, 0.99, na.rm=TRUE)

# centrality vs prevalence
hub.min30 = hub.core[which(hub.core$mean_prevalence > 30),]
corr = cor.test(hub.core$mean_prevalence, hub.core$mean_centrality, method = "spearman", exact=TRUE)
corr_min30 = cor.test(hub.min30$mean_prevalence, hub.min30$mean_centrality, method = "spearman", exact=TRUE)

scatter.corr = ggplot(hub.core, aes(x=mean_prevalence, y=mean_centrality, colour=Status)) +
  geom_point_rast(size=1, alpha=0.5) +
  geom_hline(yintercept = hub_thresh, linetype = "dashed") +
  geom_vline(xintercept = core_thresh, linetype = "dashed") +
  theme_classic() +
  scale_colour_manual(values=c("steelblue", "palegreen4")) +
  ylab("Centrality") +
  xlab("Prevalence (%)") +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size=10))

# define central and core
hub.core = hub.core[which(hub.core$mean_prevalence != 0 | !is.na(hub.core$mean_centrality)),]
hub.core$Hub_class = ifelse(hub.core$mean_centrality > hub_thresh, "Network hub", NA)
hub.core$Core_class = ifelse(hub.core$mean_prevalence > core_thresh, "Core", NA)
hub.df = hub.core[,c("Species_rep", "Hub_class")]
hub.df = hub.df[which(!is.na(hub.df$Hub_class)),]
colnames(hub.df)[2] = "Variable"
core.df = hub.core[,c("Species_rep", "Core_class")]
core.df = core.df[which(!is.na(core.df$Core_class)),]
colnames(core.df)[2] = "Variable"

# merge data
all.long = rbind(disease.df, hub.df, core.df)
all.matrix = acast(Species_rep ~ Variable, data=all.long, fun.aggregate=length)
all.matrix.meta = merge(all.matrix, species.data, by="row.names")
rownames(all.matrix.meta) = all.matrix.meta$Species_rep
upset.df = all.matrix.meta[,c("Core", "Disease", "Health", "Network hub", "Status", "Genus", "Family", "Species")]
upset.df$Genus = gsub("g__", "", upset.df$Genus)

# function to generate upset
generate_upset = function(input_df) {
  
  # filter taxa
  taxon = "Genus"
  upset.df = input_df
  upset.df = upset.df[rowSums(upset.df[,1:4]) > 1,]
  keep.taxa = names(sort(table(upset.df[,taxon]), decreasing=TRUE)[1:12])
  upset.df[,taxon] = ifelse(upset.df[,taxon] %in% keep.taxa, upset.df[,taxon], "Other")
  upset.df[,taxon] = factor(upset.df[,taxon], levels=c(sort(keep.taxa), "Other"))
  
  # upset plot
  upset.plot = upset(upset.df, colnames(upset.df)[1:4], n_intersections=10, name="",
                     base_annotations=list('Intersection size'=intersection_size(
                       mapping=aes(fill=Status), width=0.5) 
                       + scale_fill_manual(values=c('Uncultured'='palegreen4', 'Cultured'='steelblue'), name="Status")),
                     annotations = list(
                       taxon=(
                         ggplot(mapping=aes(fill=Genus))
                         + geom_bar(stat='count', position='fill', width=0.5)
                         + scale_y_continuous(labels=scales::percent_format())
                         + scale_fill_manual(values=c(brewer.pal(12, "Set3"), "darkgrey"), name=taxon)
                         + ylab('% of species'))),
                     set_sizes = upset_set_size() + ylab('Number of species'),
                     width_ratio=0.25,
                     themes=upset_default_themes(text=element_text(size=12)))
  return(upset.plot)
}

upset.plot = generate_upset(upset.df)
corr.comb = ggarrange(scatter.corr, upset.plot, widths=c(1,0.8), labels=c("C", "D"), font.label=list(size=18))
