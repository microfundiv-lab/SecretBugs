# load libraries
library(data.table)
library(CoDaSeq)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)

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
genus_stats = function(genus) {
  ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax.df = read.delim("../metadata/species_uhgg_v1.2.tsv")
  tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
  rownames(tax.df) = tax.df$Species_rep
  genus.species = tax.df[which(tax.df$Genus == genus),"Species_rep"]
  
  # calculate centred log-ratio
  abund.data.genus = abund.data.df[tax.df$Species_rep,]
  abund.data.genus$Genus = tax.df[rownames(abund.data.genus), "Genus"]
  abund.data.genus = abund.data.genus[which(abund.data.genus$Genus != "g__"),]
  abund.data.agg = aggregate(. ~ Genus, data=abund.data.genus, FUN=sum)
  rownames(abund.data.agg) = abund.data.agg$Genus
  abund.data.agg = abund.data.agg[,-1]
  abund.data.agg.clr = codaSeq.clr(abund.data.agg + 0.5, samples.by.row = FALSE)
  genus.abund = t(abund.data.agg.clr[genus,,drop=FALSE])
  genus.prev = abund.data.agg[genus,]
  genus.prev[genus.prev > 0] = 1
  genus.prev = t(genus.prev)
  return(list(
    genus_abund = as.data.frame(genus.abund),
    genus_prev  = as.data.frame(genus.prev)
  ))
}

cag170.stats = genus_stats("g__CAG-170")

# add to metadata
metadata$CAG170_abund = cag170.stats$genus_abund[metadata$Sample,]
metadata$CAG170_prev = cag170.stats$genus_prev[metadata$Sample,]

# average duplicate subject + weeks
metadata = metadata %>%
  group_by(Subject, week_num) %>%
  summarise(CAG170_abund = mean(CAG170_abund, na.rm = TRUE),
            across(everything(), ~ first(.)), .groups = "drop")

# aggregate prevalence/abundance by time and health status
keep_weekn = names(which(table(metadata$week_num) >= 10))

cag170.prev = aggregate(CAG170_prev ~ Subject, data=metadata, FUN=sum)
keep_subjects_min1 = cag170.prev[which(cag170.prev$CAG170_prev > 0),"Subject"]
keep_subjects_min11 = cag170.prev[which(cag170.prev$CAG170_prev > 10),"Subject"]

cag170.prev.agg = aggregate(CAG170_prev ~ week_num + Health.state, data=metadata, FUN = function(x) sum(x) / length(x))
cag170.abund.agg = aggregate(CAG170_abund ~ week_num + Health.state, data=metadata, FUN = mean)

cag170.prev.filt = cag170.prev.agg[which(cag170.prev.agg$week_num %in% keep_weekn),]
cag170.abund.filt = cag170.abund.agg[which(cag170.abund.agg$week_num %in% keep_weekn),]

# calculate corrs
prev.health = cag170.prev.filt[which(cag170.prev.filt$Health.state == "Healthy"),]
prev.health.corr = cor.test(prev.health$week_num, prev.health$CAG170_prev)
prev.disease = cag170.prev.filt[which(cag170.prev.filt$Health.state == "Diseased"),]
prev.disease.corr = cor.test(prev.disease$week_num, prev.disease$CAG170_prev)

abund.health = cag170.abund.filt[which(cag170.abund.filt$Health.state == "Healthy"),]
abund.health.corr = cor.test(abund.health$week_num, abund.health$CAG170_abund)
abund.disease = cag170.abund.filt[which(cag170.abund.filt$Health.state == "Diseased"),]
abund.disease.corr = cor.test(abund.disease$week_num, abund.disease$CAG170_abund)

# per subject stability
subject.btwn = metadata %>%
  arrange(Subject, week_num) %>% 
  group_by(Subject) %>%
  mutate(abundance_diff = CAG170_abund - lag(CAG170_abund)) %>% ungroup()
subject.btwn$Health.state = metadata[match(subject.btwn$Subject, metadata$Subject),"Health.state"]

# filter for min CAG-170 prev
subject.btwn.min1 = subject.btwn[which(subject.btwn$Subject %in% keep_subjects_min1),]
subject.btwn.min11 = subject.btwn[which(subject.btwn$Subject %in% keep_subjects_min11),]

# build plots
med1 = median(subject.btwn.min1$abundance_diff, na.rm = TRUE)
btwn.plot.min1 = ggplot(subject.btwn.min1, aes(x = abundance_diff)) +
  geom_histogram(fill = "darkgreen", alpha = 0.6) +
  theme_classic() +
  xlab("CAG-170 abundance difference\nbetween consecutive timepoints") +
  ylab("Frequency") +
  geom_vline(xintercept = med1, linetype = "dashed") +
  annotate("text", x = 1, y = Inf, label = paste0("Median = ", signif(med1, digits=2)), vjust = 1.5, hjust = 0, size = 4.5) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

med11 = median(subject.btwn.min11$abundance_diff, na.rm = TRUE)
btwn.plot.min11 = ggplot(subject.btwn.min11, aes(x = abundance_diff)) +
  geom_histogram(fill = "darkgreen", alpha = 0.6) +
  theme_classic() +
  xlab("CAG-170 abundance difference\nbetween consecutive timepoints") +
  ylab("Frequency") +
  geom_vline(xintercept = med11, linetype = "dashed") +
  annotate("text", x = 1, y = Inf, label = paste0("Median = ", signif(med11, digits=2)), vjust = 1.5, hjust = 0, size = 4.5) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
subject.plot = ggarrange(btwn.plot.min1, btwn.plot.min11, ncol=2, align="h", labels=c("A", "B"), font.label=list(size=18))

# combine plots and save
source("../../scripts/alex/hmp2_cag170-genes_plots.R")
ggarrange(subject.plot, gene.plot, ncol=1, labels=c("", "C"), heights=c(1,1), font.label=list(size=18))
ggsave(filename = "../figures/cag170_stab-genes.pdf", height=9, width=10)
