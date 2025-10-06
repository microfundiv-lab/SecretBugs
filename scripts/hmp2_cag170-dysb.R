# load libraries
library(data.table)
library(CoDaSeq)
library(tidyr)
library(grid)
library(ggplot2)
library(nlme)
library(RColorBrewer)
library(ggsignif)
library(ggpubr)

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
cag.present = colnames(abund.data.agg)[which(abund.data.agg["g__CAG-170",] > 0)]
abund.data.agg.clr = codaSeq.clr(abund.data.agg + 0.5, samples.by.row = FALSE)
cag.relab = data.frame(abund.data.agg.clr["g__CAG-170",])

# prepare correlation dataset
corr.select = cag.relab
corr.select$Dysb = metadata[rownames(corr.select),"Dysb_score"]
corr.select$Subject = metadata[rownames(corr.select),"Subject"]
corr.fi = reshape2::melt(corr.select, id.var=c("Dysb", "Subject"))
corr.fi$Class = ifelse(corr.fi$Dysb > quantile(corr.fi$Dysb, 0.5), "High", "Low")

# correlate dysbiosis score for all subjects
cor_test = lme(value ~ Dysb, random = ~1|Subject, data=corr.fi)
cor_test_pvalue = anova(cor_test)$`p-value`[2]
cor_test_pvalue = ifelse(cor_test_pvalue == 0, 2.2e-16)
grob = grobTree(textGrob(paste("Beta =", signif(cor_test$coefficients$fixed[2], digits=3), "\nP =", signif(cor_test_pvalue, digits=2)), x=0.7,  y=0.8, hjust=0, gp=gpar(col="black", fontsize=12)))
corr.plot = corr.fi[which(corr.fi$value > 0),]
dysb.scatter = ggplot(corr.plot, aes(x=Dysb, y=value)) +
  geom_point(size=0.8, alpha=0.8) +
  geom_smooth(aes(colour = NULL), method="lm", fullrange=TRUE) +
  theme_classic() +
  ylab("CAG-170 abundance (CLR)") +
  xlab("Dysbiosis score") +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12)) +
  annotation_custom(grob)

# compare dysbiosis class for all subjects
corr.fi$Class = ifelse(corr.fi$Dysb > quantile(corr.fi$Dysb, 0.5), "High", "Low")
cor_test = lme(value ~ Class, random = ~1|Subject, data=corr.fi)
cor_test_pvalue = anova(cor_test)$`p-value`[2]
dysb.boxplot = ggplot(corr.plot, aes(x=Class, y=value, fill=Class)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(alpha=0.2, size=0.5, position = position_jitterdodge(jitter.width = 0.5)) +
  geom_signif(comparisons = list(c("High", "Low")), annotations=signif(cor_test_pvalue, digits=3)) +
  theme_classic() +
  ylab("CAG-170 abundance (CLR)") +
  xlab("Dysbiosis classification") +
  guides(fill="none") +
  scale_fill_manual(values=c("tomato", "steelblue"), name="Dysbiosis class") +
  scale_x_discrete(limits=c("Low", "High")) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12))

dysb.comb = ggarrange(dysb.scatter, dysb.boxplot, widths=c(2,1), labels=c("C", "D"), font.label = list(size=18), align="h")