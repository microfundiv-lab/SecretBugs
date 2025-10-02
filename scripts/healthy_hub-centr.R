# load libraries
library(reshape2)
library(tidyr)
library(igraph)
library(data.table)
library(scales)

# set working dir
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/fastspar/exp06_corr/")

# parse taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.df = read.delim("../../metadata/species_uhgg_v1.2.tsv", stringsAsFactors = FALSE)
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
rownames(tax.df) = tax.df$Species_rep
tax.df$Taxon = NA
for (n in 1:nrow(tax.df)){
  for (r in ranks){
    if(!grepl("__$", tax.df[n,r])){
      tax.df[n,"Taxon"] = tax.df[n,r] }}}

get_central = function(continent){
  
  # read fastspar data
  cat(paste0("Running function for ", continent,"\n"))
  d.corr_full = read.table(paste0(continent, "_corr.tsv"), sep='\t', header=TRUE, comment.char='', row.names=1, check.names=FALSE, stringsAsFactors = FALSE)
  d.pval_full = read.table(paste0(continent, "_pvalues.tsv"), sep='\t', header=TRUE, comment.char='', row.names=1, check.names=FALSE, stringsAsFactors = FALSE)
  
  # mask upper triangle with NAs then melt, excluding upper triangle values
  d.corr_full[upper.tri(d.corr_full, diag=TRUE)] = NA
  d.pval_full[upper.tri(d.pval_full, diag=TRUE)] = NA
  
  d.corr = reshape2::melt(as.matrix(d.corr_full), na.rm=TRUE, varnames=c('otu_1', 'otu_2'), value.name='correlation')
  d.pval = reshape2::melt(as.matrix(d.pval_full), na.rm=TRUE, varnames=c('otu_1', 'otu_2'), value.name='pvalue')
  
  # merge correlations and pvalues
  d = merge(d.corr, d.pval)
  d$FDR = p.adjust(d$pvalue, method="BH")
  d = d[which(d$FDR < 0.05),]
  if (nrow(d) > 0) { 
    d$otu_1 = as.character(d$otu_1)
    d$otu_2 = as.character(d$otu_2)
    d = d[,c("otu_1", "otu_2", "correlation")]
    
    # construct graph
    d.dist = graph.data.frame(d, directed=FALSE)
    E(d.dist)$weight = 1-abs(E(d.dist)$correlation)
    central.df = data.frame(matrix(0, ncol=3, nrow=length(V(d.dist))))
    colnames(central.df) = c("betweenness", "closeness", "degree")
    rownames(central.df) = names(V(d.dist))
    central.df$species = names(V(d.dist))
    central.df$betweenness = -betweenness(d.dist)
    central.df$rscl_btwn = rescale(central.df$betweenness)
    central.df$closeness = closeness(d.dist)
    central.df$rscl_clsn = rescale(central.df$closeness)
    central.df$degree = degree(d.dist)
    central.df$rscl_degree = rescale(central.df$degree)
    central.df$mean_centrality = rowMeans(central.df[,c("rscl_btwn", "rscl_clsn", "rscl_degree")])
    central.df$continent = continent
    return(central.df)
  } else { return(NULL)
  }
}

# run function per continent
continents = c("Africa", "Asia", "Europe", "North_America", "South_America")
results.list = lapply(continents, function(x) { get_central(x) })
results.df = rbindlist(results.list)
results.cast = as.data.frame(acast(results.df, species ~ continent, value.var="mean_centrality"))
results.cast[is.na(results.cast)] = 0
results.cast$mean_centrality = as.vector(rowMeans(results.cast))
results.agg.tax = merge(results.cast, tax.df, by="row.names")
results.agg.tax = results.agg.tax[order(results.agg.tax$mean_centrality, decreasing=TRUE),]

# top 1%
results.top1perc = results.agg.tax[which(results.agg.tax$mean_centrality > quantile(results.agg.tax$mean_centrality, 0.99)),]
top.genera = sort(table(results.top1perc$Genus), decreasing=TRUE)[1:5]

# check correlation between centrality measures
btwn_clsn = cor.test(results.df$rscl_btwn, results.df$rscl_clsn, method="spearman", exact=FALSE)
btwn_degree = cor.test(results.df$rscl_btwn, results.df$rscl_degree, method="spearman", exact=FALSE)
clsn_degree = cor.test(results.df$rscl_clsn, results.df$rscl_degree, method="spearman", exact=FALSE)

# save tables
write.table(results.agg.tax, file = "../healthy_analysis/central_mean_corr02_exp06-corr.tsv", sep="\t", row.names=FALSE, quote=FALSE)