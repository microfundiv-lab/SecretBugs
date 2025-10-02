# load libraries
library(ggplot2)
library(ggpubr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/healthy_analysis/")
key.data1 = read.delim("centrality_mean_corr01.tsv")
key.data2 = read.delim("centrality_mean_corr02.tsv")

# filter for keystone
key.data1 = key.data1[which(key.data1$mean_centrality > quantile(key.data1$mean_centrality, 0.99, na.rm=TRUE)),]
key.data2 = key.data2[which(key.data2$mean_centrality > quantile(key.data2$mean_centrality, 0.99, na.rm=TRUE)),]

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
    scale_fill_gradient2(low="lavenderblush", mid="#C39BD3", high="darkorchid4", midpoint = ((max(heat$Value)-min(heat$Value))/2)+min(heat$Value), name = label) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.y = element_text(size = 12),
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
key.plot1 = heat_bar(key.data1, "mean_centrality", "Centrality")
key.plot2 = heat_bar(key.data2, "mean_centrality", "Centrality")
ggarrange(key.plot1, key.plot2, nrow=1, widths=c(1,1), labels=c("A", "B"), font.label = list(size=18))
ggsave(file="../figures/centrality_strict.pdf", width=13, height=7)
