# load libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(rstatix)
library(ggrepel)
library(ggpubr)

# set working directory
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/ml_analysis/all_diseases_pairwise")

# load performance data
gen_ml_df = function(x) {
  input.files =  list.files(path = paste0(getwd(), "/", x), pattern = "results.csv", recursive = FALSE, full.names=TRUE)
  pf.list = lapply(input.files, function(i) {
    filename = basename(i)
    disease1 = strsplit(filename, "_")[[1]][1]
    disease2 = strsplit(filename, "_")[[1]][3]
    df = read.csv(i); df$Training = disease1; df$Testing = disease2
    return(df)
  })
}

# combine all results into a dataframe
pf.cult.combined = as.data.frame(rbindlist(gen_ml_df("cult_only/results")))[,c("AUC", "Training", "Testing")]
pf.cult.combined$Dataset = "Cultured only"
pf.both.combined = as.data.frame(rbindlist(gen_ml_df(".")))[,c("AUC", "Training", "Testing")]
pf.both.combined$Dataset = "Cultured + Uncultured"
pf.rbind = rbind(pf.cult.combined, pf.both.combined)

# perform stats
wilcox_results <- pf.rbind %>%
  group_by(Training, Testing) %>%
  wilcox_test(AUC ~ Dataset) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")        
ns.results = wilcox_results[which(wilcox_results$p.adj >= 0.05),c("Training", "Testing")]

# aggregate for plotting
pf.cult.agg = aggregate(AUC ~ Training + Testing + Dataset, FUN=median, data=pf.cult.combined)
pf.both.agg = aggregate(AUC ~ Training + Testing + Dataset, FUN=median, data=pf.both.combined)
pf.merge = merge(pf.cult.agg, pf.both.agg, by=c("Training", "Testing"))
pf.merge$Diff = round(pf.merge$AUC.x-pf.merge$AUC.y, digits=3)
pf.merge$Diff_class = ifelse(pf.merge$Diff > 0, "Higher without uncultured", "Higher with uncultured")
pf.merge$Diff_class = ifelse(paste(pf.merge$Testing, pf.merge$Training) %in% paste(ns.results$Testing, ns.results$Training) | abs(pf.merge$Diff) < 0.001, "Non-significant", pf.merge$Diff_class)

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

pf.agg.ren1 = rename_disease(pf.merge, "Training")
pf.agg.fi = rename_disease(pf.agg.ren1, "Testing")

# quantify benefit of uncultured
pf.agg.fi = pf.agg.fi[which(pf.agg.fi$AUC.y > 0.6),]
prop.class = as.data.frame(table(pf.agg.fi$Diff_class))

# plot heatmap
diff.heat = ggplot(pf.agg.fi, aes(x=Training, y=Testing, fill=Diff_class, label=abs(Diff))) +
  geom_tile(alpha=0.5) +
  geom_text() +
  theme_classic() +
  scale_fill_manual(values=c("darkgreen", "steelblue", "darkgrey"), name="AUROC difference") +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(size=12, angle=45, vjust=1, hjust=1)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_text(size=14))

# plot barplot
bar.class = ggplot(prop.class, aes(x=Var1, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", alpha=0.5) +
  theme_classic() +
  ylab("Number of cross-disease comparisons\n(AUROC >0.6)") +
  scale_fill_manual(values=c("darkgreen", "steelblue", "darkgrey"), name="AUROC difference") +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size=12, angle=45, vjust=1, hjust=1)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=14))

# arrange plots
ggarrange(diff.heat, bar.class, ncol=2, widths=c(2,0.5), align="h", labels=c("A", "B"), font.label = list(size=18))
