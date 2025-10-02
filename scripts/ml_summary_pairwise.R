# load libraries
library(tidyverse)
library(data.table)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)

# set working directory
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/ml_analysis/all_diseases_pairwise/")

# load performance data
input.files =  list.files(pattern = "results.csv", recursive = FALSE, full.names=TRUE)
pf.list = lapply(input.files, function(i) {
  filename = strsplit(i, "/")[[1]][2]
  disease1 = strsplit(filename, "_")[[1]][1]
  disease2 = strsplit(filename, "_")[[1]][3]
  df = read.csv(i); df$Training = disease1; df$Testing = disease2
  return(df)
})

# combine all results into a dataframe
pf.combined = as.data.frame(rbindlist(pf.list))
pf.agg = aggregate(AUC ~ Training + Testing, FUN=median, data=pf.combined)
pf.agg$AUC = round(pf.agg$AUC, digits=2)
pf.agg$Class = ifelse(pf.agg$AUC >= 0.7, ">=0.7", ifelse(pf.agg$AUC > 0.6, ">=0.6", ifelse(pf.agg$AUC >0.5, ">=0.5", "<0.5")))

# rename diseases
rename_disease = function(df, column) {
  df[,column] = gsub("CRC", "Colorectal cancer", df[,column])
  df[,column] = gsub("CD", "Crohn's disease", df[,column])
  df[,column] = gsub("UC", "Ulcerative colitis", df[,column])
  df[,column] = gsub("T1D", "Type 1 diabetes", df[,column])
  df[,column] = gsub("T2D", "Type 2 diabetes", df[,column])
  df[,column] = gsub("MS", "Multiple sclerosis", df[,column])
  df[,column] = gsub("ME", "Myalgic encephalomyelitis", df[,column])
  df[,column] = gsub("T2D", "Type 2 diabetes", df[,column])
  df[,column] = gsub("AS", "Ankylosing spondylitis", df[,column])
  df[,column] = gsub("RA", "Rheumatoid arthritis", df[,column])
  df[,column] = gsub("Parkinson", "Parkinson's disease", df[,column])
  df[,column] = gsub("BoneDisease", "Bone disease", df[,column])
  df2 = df
  return(df2)
}

pf.agg.ren1 = rename_disease(pf.agg, "Training")
pf.agg.fi = rename_disease(pf.agg.ren1, "Testing")

# plot heatmap
ml.heat = ggplot(pf.agg.fi, aes(x=Training, y=Testing, fill=Class, label=AUC)) +
  geom_tile() +
  geom_text() +
  theme_classic() +
  scale_fill_manual(values=c("grey90", "lightblue", "steelblue", "pink2"), name="AUROC") +
  theme(axis.text.x = element_text(size=12, angle=45, vjust=1, hjust=1)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14))