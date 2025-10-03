# load libraries
library(ggplot2)
library(ggpubr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/healthy_analysis/")
prev.data = read.delim("prevalence_mean.tsv")
prev.data = prev.data[which(prev.data$mean_prevalence > 70),]
key.data = read.delim("central_mean.tsv")
key.data = key.data[which(key.data$mean_centrality > quantile(key.data$mean_centrality, 0.99, na.rm=TRUE)),]

# define function
heat_bar = function(df, metric, label) {
  
  # clean names
  df$Taxon = gsub("s__","",df$Taxon)
  df$Taxon = gsub("g__","",df$Taxon)
  df$Taxon = gsub("_"," ",df$Taxon)
  
  # order taxon
  order.taxon = order(df[,metric], decreasing=FALSE)
  df$Taxon = factor(df$Taxon, levels=df$Taxon[order.taxon])
  
  # parse table
  cont_cols = c("Africa", "Asia", "Europe", "North_America", "Oceania", "South_America")
  select.cols = colnames(df)[which(colnames(df) %in% cont_cols)]
  heat = reshape2::melt(df[, c("Taxon", select.cols)], id.vars = "Taxon", variable.name = "Continent", value.name = "Value")
  heat$Continent = gsub("_", " ", heat$Continent)
  
  # heatmap
  p_heat = ggplot(heat, aes(x = Continent, y = Taxon, fill = Value)) +
    geom_tile() +
    scale_fill_gradient(low="lavenderblush", high="darkorchid4", name = label) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 12, face = "italic"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      legend.position = "bottom")
  
  # barplot
  p_bar = ggplot(df, aes(x = Taxon, y = !!sym(metric), fill=Status)) +
    geom_bar(stat="identity", width = 0.8, alpha=0.8) +
    coord_flip(ylim=c(min(df[,metric])*0.99,max(df[,metric]))) +
    scale_fill_manual(values=c("#b4c8ff", "#a5e27f")) +
    labs(y = label, x = NULL) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size=14))
  
  # combine and return
  comb.plot = ggarrange(p_heat, p_bar, align="h", widths=c(1,0.5))
  return(comb.plot)
}

# run functions and combine
prev.plots = heat_bar(prev.data, "mean_prevalence", "Prevalence (%)")
key.plots = heat_bar(key.data, "mean_centrality", "Centrality")
heat.comb = ggarrange(prev.plots, key.plots, nrow=1, widths=c(1,0.9), labels=c("A", "B"), font.label = list(size=18))

# combine with correlation
source("../../scripts/alex/healthy_prev-hub_corr.R")
ggarrange(heat.comb, corr.comb, ncol=1, heights=c(1,0.8))
ggsave(file="figures/healthy_analysis.pdf", width=14, height=12)
