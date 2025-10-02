# load libraries
library(data.table)
library(stringr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(tidyverse)

# set working directory
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/ml_analysis/pool-study/")

# load performance data
input.files =  list.files(pattern = "results.csv", recursive = FALSE, full.names=TRUE)
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

# select best
pf.agg = aggregate(AUC ~ Dataset + method, data=pf.combined, FUN=median)
pf.agg = pf.agg[order(pf.agg$AUC, decreasing=TRUE),]
order.method = as.vector(unique(pf.agg[,"method"]))
pf.combined$method = factor(pf.combined$method, levels=rev(order.method))
pf.combined$method = as.vector(pf.combined$method)

# plot boxplot by method
comparisons = combn(unique(pf.combined$method), 2, simplify = FALSE)
ml.box = ggplot(pf.combined, aes(x=method, y=AUC, fill=method)) +
  geom_boxplot(alpha=0.5, outlier.colour=NA, width=0.5) +
  geom_point(alpha=0.2, size=0.3, position = position_jitter(width=0.1)) + 
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format"
  ) +
  ylab("AUROC") +
  scale_x_discrete(limits=c("Ridge Regression", "Gradient Boosting", "Random Forest"),
                   labels=c("Ridge\nRegression", "Gradient\nBoosting", "Random\nForest")) +
  scale_fill_manual(values=c("#0072B2", "#009E73", "#E69F00"), name="Method") +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.text.y = element_text(size=12))

# kruskal-wallis test for method
kruskal.model = kruskal.test(AUC ~ method, data=pf.combined)

# plot summary statistics split by disease in uncultured
pf.agg.median = aggregate(AUC ~ Dataset + method + Disease, data=pf.combined, FUN=median)
pf.agg.uncult = pf.agg.median[which(pf.agg.median$Dataset == "Uncultured"),]
source("../../../scripts/alex/ml_summary_pairwise.R")
pf.agg.uncult = rename_disease(pf.agg.uncult, "Disease")
pf.agg.uncult = pf.agg.uncult %>%
  group_by(Disease) %>%
  mutate(is_top = AUC == max(AUC, na.rm = TRUE)) %>%
  ungroup()

pf.agg.uncult$AUC = round(pf.agg.uncult$AUC, digits=2)
ml.heat = ggplot(pf.agg.uncult, aes(y=Disease, x=method, fill=AUC)) +
  geom_tile(color = NA, linewidth = 0.4) +
  geom_text(
    aes(label = AUC, fontface = ifelse(is_top, "bold", "plain")),
    colour = "black"
  ) +
  scale_fill_gradient(low = "white", high = "#c2a5cf", name="AUROC") +  # gradient for AUC
  scale_x_discrete(limits=c("Ridge Regression", "Gradient Boosting", "Random Forest"),
                   labels=c("Ridge\nRegression", "Gradient\nBoosting", "Random\nForest")) +
  scale_y_discrete(limits = rev(unique(pf.agg.uncult$Disease))) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12, angle=45, vjust=1, hjust=1),
        axis.text.y = element_text(size=12),
        axis.title = element_blank())

# combine plots
ggarrange(ml.box, ml.heat, widths=c(0.5,1), labels=c("A", "B"), font.label = list(size=18))
ggsave(file="../../figures/ml_summary_best-method.pdf", width=10, height=6)