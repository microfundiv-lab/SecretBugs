# load libraries
library(data.table)
library(CoDaSeq)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(cvequality)

# input and output
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/timeseries/")
abund.data = fread("../bwa/bwa_counts-filtered_samples.csv")
abund.data.df = data.frame(abund.data[,-1], check.names = FALSE)
output.name = "hmp2_dysbiosis"
rownames(abund.data.df) = abund.data$Genome

# filter metadata
metadata = read.delim("../metadata/hmp2_ibd/metagenomes_dysb-score.tsv")
metadata = unique(metadata[,-1])
rownames(metadata) = metadata$Sample

metadata$Disease.name = ifelse(metadata$Health.state == "Healthy", "Healthy", metadata$Disease.name)
abund.data.df = abund.data.df[,metadata$Sample]

# load taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.df = read.delim("../metadata/species_uhgg_v1.2.tsv")
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
rownames(tax.df) = tax.df$Species_rep
cag.species = tax.df[which(tax.df$Genus == "g__CAG-170"),"Species_rep"]

# calculate centred log-ratio
abund.data.genus = abund.data.df[tax.df$Species_rep,]
abund.data.genus$Genus = tax.df[rownames(abund.data.genus), "Genus"]
abund.data.genus = abund.data.genus[which(abund.data.genus$Genus != "g__"),]
abund.data.agg = aggregate(. ~ Genus, data=abund.data.genus, FUN=sum)
rownames(abund.data.agg) = abund.data.agg$Genus
abund.data.agg = abund.data.agg[,-1]
abund.data.agg.clr = codaSeq.clr(abund.data.agg + 0.5, samples.by.row = FALSE)
cag170.abund = t(abund.data.agg.clr["g__CAG-170",,drop=FALSE])
cag170.prev = abund.data.agg["g__CAG-170",]
cag170.prev[cag170.prev > 0] = 1
cag170.prev = t(cag170.prev)

# add to metadata
metadata$CAG170_abund = cag170.abund[metadata$Sample,]
metadata$CAG170_prev = cag170.prev[metadata$Sample,]

# aggregate prevalence/abundance by time and health status
keep_weekn = names(which(table(metadata$week_num) > 10))

cag170.prev.agg = aggregate(CAG170_prev ~ week_num + Health.state, data=metadata, FUN = function(x) sum(x) / length(x))
cag170.abund.agg = aggregate(CAG170_abund ~ week_num + Health.state, data=metadata, FUN = mean)

cag170.prev.filt = cag170.prev.agg[which(cag170.prev.agg$week_num %in% keep_weekn),]
cag170.abund.filt = cag170.abund.agg[which(cag170.abund.agg$week_num %in% keep_weekn),]

# plot
abund.plot = ggplot(cag170.abund.filt, aes(x=week_num, y=CAG170_abund, colour=Health.state)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm") +
  xlab("Week number") +
  ylab("CAG-170 mean abundance (CLR)") +
  scale_colour_manual(values=c("tomato", "steelblue"), name="Health status") +
  theme_classic() +
  theme(axis.title = element_text(size=14)) +
  theme(axis.text = element_text(size=12))

prev.plot = ggplot(cag170.prev.filt, aes(x=week_num, y=CAG170_prev*100, colour=Health.state)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm") +
  xlab("Week number") +
  ylab("CAG-170 prevalence (%)") +
  scale_colour_manual(values=c("tomato", "steelblue"), name="Health status") +
  theme_classic() +
  theme(axis.title = element_text(size=14)) +
  theme(axis.text = element_text(size=12))

# combine plots
ggarrange(prev.plot, abund.plot, ncol=1, labels=c("A", "B"), font.label = list(size=16), common.legend = TRUE)
ggsave(file = "../figures/cag170_longitudinal.pdf", height=6, width=8)

# calculate corr
prev.health = cag170.prev.filt[which(cag170.prev.filt$Health.state == "Healthy"),]
prev.health.corr = cor.test(prev.health$week_num, prev.health$CAG170_prev)
prev.disease = cag170.prev.filt[which(cag170.prev.filt$Health.state == "Diseased"),]
prev.disease.corr = cor.test(prev.disease$week_num, prev.disease$CAG170_prev)

abund.health = cag170.abund.filt[which(cag170.abund.filt$Health.state == "Healthy"),]
abund.health.corr = cor.test(abund.health$week_num, abund.health$CAG170_abund)
abund.disease = cag170.abund.filt[which(cag170.abund.filt$Health.state == "Diseased"),]
abund.disease.corr = cor.test(abund.disease$week_num, abund.disease$CAG170_abund)
