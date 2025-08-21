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
pf.delta = pf.combined[which(pf.combined$Dataset %in% c("Cultured", "Uncultured") & pf.combined$method == order.method[1]),]
pf.delta.agg = aggregate(AUC ~ Dataset + Disease, data=pf.delta, FUN=median)
pf.delta.agg$Delta = ave(as.numeric(pf.delta.agg$AUC), factor(pf.delta.agg$Disease), FUN=function(x) c(NA,diff(x)))
pf.delta.agg = pf.delta.agg[!is.na(pf.delta.agg$Delta),]

# get stats per disease
gtype.stats = pf.delta %>%
  group_by(Disease) %>%
  summarise(
    test = list(wilcox.test(AUC ~ Dataset)),
    median_cultured = median(AUC[Dataset == "Cultured"]),
    median_uncultured = median(AUC[Dataset == "Uncultured"]),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(test, `[[`, "p.value"),
    W_stat  = sapply(test, function(x) x$statistic),
    direction = case_when(
      median_cultured > median_uncultured ~ "Higher in cultured",
      median_cultured < median_uncultured ~ "Higher in uncultured",
      TRUE ~ "No difference"
    )
  ) %>%
  select(Disease, median_cultured, median_uncultured, direction, W_stat, p_value)
gtype.stats$FDR = p.adjust(gtype.stats$p_value, method="BH")
gtype.stats$Final_result = ifelse(gtype.stats$FDR < 0.05, "Significant", "Non-significant")

# boxplot both per disease
pf.both.gb = pf.combined[which(pf.combined$Dataset == "Uncultured" & pf.combined$method == order.method[1]),]
pf.both.gb.agg = aggregate(AUC ~ Disease, data=pf.both.gb, FUN=median)
disease.order = pf.both.gb.agg[order(pf.both.gb.agg$AUC, decreasing=FALSE),"Disease"]
box.disease = ggplot(pf.both.gb, aes(x=Disease, y = AUC)) +
  geom_point(alpha=0.5, size=0.8, position = position_jitter(width=0.1)) +
  geom_boxplot(alpha=0.5, outlier.colour=NA, width=0.7, fill="steelblue") +
  coord_flip() +
  ylab("AUROC") +
  scale_x_discrete(limits=disease.order) +
  theme_classic() +
  theme(legend.text = element_text(size=10), legend.title = element_text(size=12)) +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

# bar plot delta uncultured
bar.delta = ggplot(pf.delta.agg, aes(x=Disease, y = Delta)) +
  geom_bar(stat="identity", alpha=0.5, fill="steelblue") +
  ylab(expression(Delta~"AUROC (uncultured - cultured)")) +
  coord_flip() +
  scale_x_discrete(limits=disease.order) +
  scale_y_continuous(breaks=c(-0.1,-0.05,0,0.05,0.1), limits=c(-0.12,0.05)) +
  theme_classic() +
  theme(axis.title.x = element_text(size=14)) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_blank())

# combine plots
source("../../../scripts/alex/ml_summary_cross-vs-pool.R")
ml.combined = ggarrange(box.disease, bar.delta, scatter.auc, widths=c(1,0.5,0.8), ncol=3, align="h", labels=c("a", "", "b"), font.label=list(size=18)) +
  theme(plot.margin = margin(0.2,1,0.2,0.2, "cm"))
source("../../scripts/alex/ml_summary_pairwise.R")
ggarrange(ml.combined, ml.heat, nrow=2, labels=c("", "c"), font.label=list(size=18), heights=c(0.8,1))
ggsave(file="../../figures/ml_summary.pdf", width=14, height=12)
