# load libraries
library(tidyverse)
library(data.table)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)

# set working directory
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/ml_analysis/")

# load performance data
pool.files =  list.files(path = "pool-study/", pattern = "*_both_results.csv", recursive = TRUE, full.names=TRUE)
cross.files =  list.files(path = "cross-study/", pattern = "*_both_results.csv", recursive = TRUE, full.names=TRUE)
input.files = c(pool.files, cross.files)
pf.list = lapply(input.files, function(i) {
  analysis = strsplit(i, "/")[[1]][1]
  filename = basename(i)
  disease = strsplit(filename, "_")[[1]][1]
  meta = strsplit(filename, "_")[[1]][2]
  df = read.csv(i); df$Analysis = analysis; df$Disease = disease; df$Dataset = meta
  return(df)
})

# combine all results into a dataframe
pf.combined = as.data.frame(rbindlist(pf.list))
pf.combined = pf.combined[which(pf.combined$Analysis %in% c("cross-study", "pool-study")),]

# select best parameters
selected.method = "glmnet"
pf.combined = pf.combined[which(pf.combined$method == selected.method),]

# rename diseases
pf.combined$Disease = gsub("CRC", "Colorectal cancer", pf.combined$Disease)
pf.combined$Disease = gsub("CD", "Crohn's disease", pf.combined$Disease)
pf.combined$Disease = gsub("UC", "Ulcerative colitis", pf.combined$Disease)
pf.combined$Disease = gsub("T1D", "Type 1 diabetes", pf.combined$Disease)
pf.combined$Disease = gsub("T2D", "Type 2 diabetes", pf.combined$Disease)
pf.combined$Disease = gsub("MS", "Multiple sclerosis", pf.combined$Disease)
pf.combined$Disease = gsub("ME", "Myalgic encephalomyelitis", pf.combined$Disease)
pf.combined$Disease = gsub("T2D", "Type 2 diabetes", pf.combined$Disease)
pf.combined$Disease = gsub("AS", "Ankylosing spondylitis", pf.combined$Disease)
pf.combined$Disease = gsub("RA", "Rheumatoid arthritis", pf.combined$Disease)
pf.combined$Disease = gsub("Parkinson", "Parkinson's disease", pf.combined$Disease)
pf.combined$Disease = gsub("BoneDisease", "Bone disease", pf.combined$Disease)

# aggregate by disease and analysis
pf.combined.agg = aggregate(AUC ~ Disease + Analysis, data=pf.combined, FUN=mean)
pf.combined.agg.long = as.data.frame(acast(Disease ~ Analysis, value=AUC, data=pf.combined.agg))
pf.combined.agg.long$Disease = rownames(pf.combined.agg.long)
pf.combined.agg.long = pf.combined.agg.long[which(!is.na(pf.combined.agg.long$`cross-study`)),]
pf.combined.agg.long$Ratio = pf.combined.agg.long$`cross-study`/pf.combined.agg.long$`pool-study`
pf.combined.agg.long[order(pf.combined.agg.long$`cross-study`, decreasing=TRUE),]

# plot scatter
scatter.auc = ggplot(pf.combined.agg.long, aes(x=`pool-study`, y = `cross-study`)) +
  geom_point(alpha=1, size=2) +
  geom_text_repel(aes(label = Disease), colour="black") +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="grey") +
  ylab("AUROC (cross-study)") +
  xlab("AUROC (pooled studies)") +
  xlim(0.45,1) +
  ylim(0.45,1) +
  theme_classic() +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))
