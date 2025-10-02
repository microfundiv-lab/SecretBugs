# load libraries
library(data.table)
library(stringr)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggpubr)

# set working directory
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/ml_analysis/pool-study/")

# load performance data
input.files =  list.files(pattern = "results.csv", recursive = TRUE, full.names=TRUE)
pf.list = lapply(input.files, function(i) {
  filename = strsplit(i, "/")[[1]][2]
  disease = strsplit(filename, "_")[[1]][1]
  meta = strsplit(filename, "_")[[1]][2]
  df = read.csv(i); df$Disease = disease; df$Dataset = meta
  return(df)
})
pf.combined = as.data.frame(rbindlist(pf.list))

# rename categories
pf.combined$method = gsub("glmnet", "Ridge Regression", pf.combined$method)
pf.combined$method = gsub("rf", "Random Forest", pf.combined$method)
pf.combined$method = gsub("xgbTree", "Gradient Boosting", pf.combined$method)
pf.combined$Dataset = str_to_title(pf.combined$Dataset)

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

# select best
pf.agg = aggregate(AUC ~ Dataset + method, data=pf.combined, FUN=median)
pf.agg = pf.agg[order(pf.agg$AUC, decreasing=TRUE),]
order.method = as.vector(unique(pf.agg[,"method"]))
pf.combined$method = factor(pf.combined$method, levels=rev(order.method))

# get delta values
pf.delta = pf.combined[which(pf.combined$Dataset %in% c("Cultured", "Both") & pf.combined$method == order.method[1]),]
pf.delta.agg = aggregate(AUC ~ Dataset + Disease, data=pf.delta, FUN=median)
pf.delta.agg$Delta = -ave(as.numeric(pf.delta.agg$AUC), factor(pf.delta.agg$Disease), FUN=function(x) c(NA,diff(x)))
pf.delta.agg = pf.delta.agg[!is.na(pf.delta.agg$Delta),]

# bar plot delta uncultured
pf.delta.agg$Diff_class = ifelse(pf.delta.agg$Delta > 0, "Higher with uncultured", "Higher without uncultured")
both.delta = ggplot(pf.delta.agg, aes(x=Disease, y = Delta, fill=Diff_class)) +
  geom_bar(stat="identity", alpha=0.5) +
  ylab(expression(Delta~"AUROC (both - cultured)")) +
  coord_flip() +
  scale_x_discrete(limits=disease.order) +
  scale_fill_manual(values=c("darkgreen", "steelblue"), name="AUROC difference") +
  scale_y_continuous(breaks=c(-0.05,0,0.05,0.1), limits=c(-0.05,0.075)) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_blank())

# perform stats
gtype.stats = pf.delta %>%
  group_by(Disease) %>%
  summarise(
    test = list(wilcox.test(AUC ~ Dataset)),
    median_cultured = median(AUC[Dataset == "Cultured"]),
    median_both = median(AUC[Dataset == "Both"]),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(test, `[[`, "p.value"),
    W_stat  = sapply(test, function(x) x$statistic),
    direction = case_when(
      median_cultured > median_both ~ "Higher without uncultured",
      median_cultured < median_both ~ "Higher with uncultured",
      TRUE ~ "No difference"
    )
  ) %>%
  select(Disease, median_cultured, median_both, direction, W_stat, p_value)
gtype.stats$FDR = p.adjust(gtype.stats$p_value, method="BH")
gtype.stats$Final_result = ifelse(gtype.stats$FDR < 0.05, "Significant", "Non-significant")