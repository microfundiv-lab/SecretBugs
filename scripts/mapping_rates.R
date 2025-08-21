# load libraries
library(ggplot2)
library(ggridges)
library(ggpubr)
library(ggrastr)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggpointdensity)
library(viridis)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/")
species.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
uncultured = species.metadata[which(species.metadata$Status == "Uncultured"),"Species_rep"]
metadata = read.delim("metadata/metagenomes_03-2024_samples.tsv")
rownames(metadata) = metadata$Sample

# all counts
counts.all = as.data.frame(fread("bwa/bwa_counts-total_samples.csv"))
rownames(counts.all) = counts.all$Genome
counts.all = counts.all[,-1]
counts.all = counts.all[,rownames(metadata)]

# filtered counts
counts.filt = as.data.frame(fread("bwa/bwa_counts-filtered_samples.csv"))
rownames(counts.filt) = counts.filt$Genome
counts.filt = counts.filt[,-1]
counts.filt = counts.filt[,rownames(metadata)]

# calculate mapping rate
sample.rates = data.frame(matrix(ncol=8, nrow=nrow(metadata)))
rownames(sample.rates) = metadata$Sample
colnames(sample.rates) = c("Reads_count", "Reads_raw", "Reads_filt", "Mapping_raw", "Mapping_filt", "Mapping_uncultured", "Species_total", "Species_uncultured")
sample.rates$Reads_count = metadata$Read.count
sample.rates$Reads_raw = colSums(counts.all)
sample.rates$Reads_filt = colSums(counts.filt)
sample.rates$Mapping_raw = sample.rates$Reads_raw/metadata$Read.count*100
sample.rates$Mapping_filt = sample.rates$Reads_filt/metadata$Read.count*100
sample.rates$Mapping_uncultured = colSums(counts.filt[uncultured,])/metadata$Read.count*100
sample.rates$Species_total = colSums(counts.filt > 0)
sample.rates$Species_uncultured = colSums(counts.filt[uncultured,] > 0)

# exclude healthy only
source("../scripts/alex/metadata_disease-numbers.R")
sample.rates = sample.rates[rownames(metadata.disease),]
sample.rates$Health.state = metadata.disease$Health.state

# plot mapping rates
sample.rates.melt = reshape2::melt(sample.rates)
abs.reads = sample.rates.melt[which(sample.rates.melt$variable %in% c("Reads_count", "Reads_raw", "Reads_filt")),]
reads.boxplot = ggplot(abs.reads, aes(x=variable, y=log10(value), fill=variable)) +
  geom_point_rast(colour = "darkgrey", size=0.2, position = position_jitter(width = 0.2)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  theme_classic() +
  ylab(expression("Number of reads (log" [10] * ")")) +
  scale_x_discrete(labels=c("Total", "Mapped\n(raw)", "Mapped\n(filtered)")) +
  scale_fill_manual(values=c("#d9f0d3", "#a6dba0", "#5aae61")) +
  theme(
    legend.position="none",
    axis.text = element_text(size=12),
    axis.title.y = element_text(size=14),
    axis.title.x = element_blank())

perc.mapp = sample.rates.melt[which(sample.rates.melt$variable %in% c("Mapping_raw", "Mapping_filt")),]
state.boxplot = ggplot(perc.mapp, aes(x=variable, y=value, fill=Health.state)) +
  geom_point_rast(colour = "darkgrey", size=0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width=0.7)) +
  geom_boxplot(alpha=0.7, outlier.shape=NA) +
  theme_classic() +
  ylab("% Reads mapped") +
  coord_cartesian(ylim = c(0,100)) +
  scale_fill_manual(values=c("tomato", "steelblue"), name="Health state") +
  scale_x_discrete(labels=c("Raw", "Filtered")) +
  theme(
    legend.position="top",
    legend.title = element_text(size=12),
    legend.text = element_text(size=10),
    axis.text = element_text(size=12),
    axis.title.y = element_text(size=14),
    axis.title.x = element_blank())

ggarrange(reads.boxplot, state.boxplot, labels=c("a", "b"), font.label=list(size=16), align="h")
ggsave("figures/mapping_rates.pdf", height=5, width=8)

# list of comparisons to analyse
disease.data = list()
diseases = c("Adenoma", "AS", "BoneDisease", "CD", "CRC","ME",
             "MS", "Obesity", "Parkinson", "RA", "T1D", "T2D", "UC")

for (d in diseases) {
  # check in which studies that disease is and if it also has healthy samples
  disease.studies = unique(metadata[which(metadata$Disease.name == d),"Study"])
  healthy.studies = unique(metadata[which(metadata$Disease.name == "Healthy"),"Study"])
  which_studies = intersect(disease.studies, healthy.studies)
  
  # filter metadata for those studies
  metadata.filt = metadata[which(metadata$Study %in% which_studies),]
  
  # filter so that only the specific disease + healthy are included
  metadata.filt = metadata.filt %>% filter(!is.na(Disease.name)) %>%
    filter(Disease.name == ifelse((Health.state=="Diseased" & Disease.name != d), NA, Disease.name))
  
  rates.filt = sample.rates[metadata.filt$Sample,]
  rates.filt$Disease = d
  rates.filt$Health.state = metadata.filt$Health.state
  disease.data[[d]] = rates.filt
}

# combine data frames
all.data = do.call(what = rbind, args = disease.data)

# rename diseases
all.data$Disease = gsub("CRC", "Colorectal cancer", all.data$Disease)
all.data$Disease = gsub("CD", "Crohn's disease", all.data$Disease)
all.data$Disease = gsub("UC", "Ulcerative colitis", all.data$Disease)
all.data$Disease = gsub("T1D", "Type 1 diabetes", all.data$Disease)
all.data$Disease = gsub("T2D", "Type 2 diabetes", all.data$Disease)
all.data$Disease = gsub("MS", "Multiple sclerosis", all.data$Disease)
all.data$Disease = gsub("ME", "Myalgic encephalomyelitis", all.data$Disease)
all.data$Disease = gsub("T2D", "Type 2 diabetes", all.data$Disease)
all.data$Disease = gsub("AS", "Ankylosing spondylitis", all.data$Disease)
all.data$Disease = gsub("RA", "Rheumatoid arthritis", all.data$Disease)
all.data$Disease = gsub("Parkinson", "Parkinson's disease", all.data$Disease)
all.data$Disease = gsub("BoneDisease", "Bone disease", all.data$Disease)
all.data$Disease = factor(all.data$Disease, levels=sort(unique(all.data$Disease), decreasing=TRUE))
all.melt = reshape2::melt(all.data)

# plot uncultured abundance and prevalence
mean.abund = mean(all.data$Mapping_uncultured/all.data$Mapping_filt*100)
uncult.abund = ggplot(all.data, aes(x=(Mapping_uncultured/Mapping_filt*100), y=Disease, fill=Health.state)) +
  geom_density_ridges(alpha=0.7, colour="darkgrey") +
  geom_vline(xintercept = mean.abund, linetype="dashed", colour="black") +
  theme_classic() +
  ylab("") +
  xlab("% Uncultured abundance") +
  scale_fill_manual(values=c("tomato", "steelblue"), name="Health state") +
  theme(legend.position="top", legend.text = element_text(size=10), legend.title = element_text(size=10)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14))

mean.species = mean(all.data$Species_uncultured/all.data$Species_total*100)
uncult.species = ggplot(all.data, aes(x=(Species_uncultured/Species_total*100), y=Disease, fill=Health.state)) +
  geom_density_ridges(alpha=0.7, colour="darkgrey") +
  geom_vline(xintercept = mean.species, linetype="dashed", colour="black") +
  theme_classic() +
  ylab("") +
  xlab("% Uncultured species") +
  scale_fill_manual(values=c("tomato", "steelblue"), name="Health state") +
  theme(legend.position="top", legend.text = element_text(size=10), legend.title = element_text(size=10)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14))

uncult.plot = ggarrange(uncult.species, uncult.abund, ncol=2, common.legend = TRUE, widths=c(1.5,1)) +
  theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm"))

# perform Wilcoxon test for each disease
wilcox.df = all.data %>%
  group_by(Disease) %>%
  summarise(
    wilcox_p = tryCatch({
      wilcox.test(Species_uncultured ~ Health.state)$p.value
    }, error = function(e) NA)
  )
wilcox.df$FDR = p.adjust(wilcox.df$wilcox_p)
wilcox.df = wilcox.df[which(wilcox.df$FDR < 0.05),]
