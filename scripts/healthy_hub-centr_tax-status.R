# load libraries
library(reshape2)
library(tidyr)
library(igraph)
library(data.table)
library(scales)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

# set working dir
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/fastspar/")

# parse taxonomy
ranks = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax.df = read.delim("../metadata/species_uhgg_v1.2.tsv", stringsAsFactors = FALSE)
tax.df = separate(data = tax.df, col = Lineage, sep = ";", into = ranks)
rownames(tax.df) = tax.df$Species_rep
tax.df$Taxon = NA
for (n in 1:nrow(tax.df)){
  for (r in ranks){
    if(!grepl("__$", tax.df[n,r])){
      tax.df[n,"Taxon"] = tax.df[n,r] }}}

# get health/disease association
hd.info = read.csv("../meta-analysis/species_ds.csv")
health = hd.info[which(hd.info$DS < 0), "Species_rep"]
disease = hd.info[which(hd.info$DS > 0), "Species_rep"]

# read fastspar data
neighbours_df = function(continent) {
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
  }
  
  # filter for CAG-170
  cag170 = tax.df[which(tax.df$Genus == "g__CAG-170"),"Genome"]
  d = d[d$otu_1 %in% cag170 | d$otu_2 %in% cag170,]
  d$neighbour = ifelse(d$otu_1 %in% cag170, d$otu_2, d$otu_1)
  d = d[!d$neighbour %in% cag170,]
  d$Continent = continent
  return(d)
}

cag170.africa = neighbours_df("Africa")
cag170.asia = neighbours_df("Asia")
cag170.northa = neighbours_df("North_America")
cag170.europe = neighbours_df("Europe")
cag170.southa = neighbours_df("South_America")

# combine and plot
comb.df = rbind(cag170.africa, cag170.asia, cag170.northa, cag170.europe, cag170.southa)
comb.agg = aggregate(correlation ~ neighbour, data=comb.df, FUN=median)
comb.agg$Status = tax.df[match(comb.agg$neighbour, tax.df$Species_rep),"Status"]
comb.agg$Order = tax.df[match(comb.agg$neighbour, tax.df$Species_rep),"Order"]
comb.agg$Family = tax.df[match(comb.agg$neighbour, tax.df$Species_rep),"Family"]
comb.agg$Association = "Non-significant"
comb.agg$Association = ifelse(comb.agg$neighbour %in% health, "Health", comb.agg$Association)
comb.agg$Association = ifelse(comb.agg$neighbour %in% disease, "Disease", comb.agg$Association)

# define colors for top orders
top.orders = gsub("o__", "", names(sort(table(comb.agg$Order), decreasing=TRUE)[1:12]))
order.colours = c(brewer.pal(length(top.orders), "Set3"), "darkgrey")
names(order.colours) = c(top.orders, "Other")

comb.agg$Order_class = gsub("o__", "", comb.agg$Order)
comb.agg$Order_class = ifelse(comb.agg$Order_class %in% top.orders, comb.agg$Order_class, "Other")
assoc.colours = c("steelblue", "tomato", "darkgrey")
names(assoc.colours) = c("Health", "Disease", "Non-significant")

# create df of props
pos_thresh = quantile(comb.agg$correlation, 0.90)
neg_thresh = quantile(comb.agg$correlation, 0.1)

order.pos.df = as.data.frame(table(comb.agg[which(comb.agg$correlation > pos_thresh),"Order_class"])/length(comb.agg[which(comb.agg$correlation > pos_thresh),"Order_class"])*100)
order.pos.df$Direction = "Positive"
order.neg.df = as.data.frame(table(comb.agg[which(comb.agg$correlation < neg_thresh),"Order_class"])/length(comb.agg[which(comb.agg$correlation < neg_thresh),"Order_class"])*100)
order.neg.df$Direction = "Negative"
order.df = rbind(order.pos.df, order.neg.df)
order.df$Var1 = factor(order.df$Var1, levels = unique(order.df[order(order.df$Freq, decreasing=TRUE),"Var1"]))

status.pos.df = as.data.frame(table(comb.agg[which(comb.agg$correlation > pos_thresh),"Status"])/length(comb.agg[which(comb.agg$correlation > pos_thresh),"Status"])*100)
status.pos.df$Direction = "Positive"
status.neg.df = as.data.frame(table(comb.agg[which(comb.agg$correlation < neg_thresh),"Status"])/length(comb.agg[which(comb.agg$correlation < neg_thresh),"Status"])*100)
status.neg.df$Direction = "Negative"
status.df = rbind(status.pos.df, status.neg.df)

assoc.pos.df = as.data.frame(table(comb.agg[which(comb.agg$correlation > pos_thresh),"Association"])/length(comb.agg[which(comb.agg$correlation > pos_thresh),"Association"])*100)
assoc.pos.df$Direction = "Positive"
assoc.neg.df = as.data.frame(table(comb.agg[which(comb.agg$correlation < neg_thresh),"Association"])/length(comb.agg[which(comb.agg$correlation < neg_thresh),"Association"])*100)
assoc.neg.df$Direction = "Negative"
assoc.df = rbind(assoc.pos.df, assoc.neg.df)

# tax plot
order.hist = ggplot(comb.agg, aes(x=correlation, fill=Order_class)) +
  geom_histogram(alpha=0.9, colour="darkgrey", linewidth=0.2) +
  theme_classic() +
  xlab("Correlation") +
  ylab("Number of species") +
  scale_fill_manual(values=order.colours, name="Order") +
  scale_x_continuous(breaks=c(-0.2,-0.1,0,0.1,0.20,0.30,0.40), limits=c(-0.25,0.45)) +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))

order.bar = ggplot(order.df, aes(x=Direction, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", alpha=0.9) +
  theme_classic() +
  xlab("Direction") +
  ylab("% Species within the top 10%") +
  scale_fill_manual(values=order.colours, name="Order") +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))

order.plot = ggarrange(order.hist, order.bar, widths=c(1,0.3), common.legend=TRUE, align="h", labels=c("A", ""), font.label = list(size=18))

# status plot
status.hist = ggplot(comb.agg, aes(x=correlation, fill=Status)) +
  geom_histogram(alpha=0.6, colour="darkgrey", linewidth=0.2) +
  theme_classic() +
  xlab("Correlation") +
  ylab("Number of species") +
  scale_fill_manual(values=c("steelblue", "darkgreen"), name="Status") +
  scale_x_continuous(breaks=c(-0.2,-0.1,0,0.1,0.20,0.30,0.40), limits=c(-0.25,0.45)) +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))

status.bar = ggplot(status.df, aes(x=Direction, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", alpha=0.6) +
  theme_classic() +
  xlab("Direction") +
  ylab("% Species within the top 10%") +
  scale_fill_manual(values=c("steelblue", "darkgreen"), name="Status") +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))

status.plot = ggarrange(status.hist, status.bar, widths=c(1,0.3), common.legend=TRUE, align="h", labels=c("B", ""), font.label = list(size=18))

# health plot
assoc.hist = ggplot(comb.agg, aes(x=correlation, fill=Association)) +
  geom_histogram(alpha=0.6, colour="darkgrey", linewidth=0.2) +
  theme_classic() +
  xlab("Correlation") +
  ylab("Number of species") +
  scale_fill_manual(values=assoc.colours, name="Association") +
  scale_x_continuous(breaks=c(-0.2,-0.1,0,0.1,0.20,0.30,0.40), limits=c(-0.25,0.45)) +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))

assoc.bar = ggplot(assoc.df, aes(x=Direction, y=Freq, fill=Var1)) +
  geom_bar(stat="identity", alpha=0.6) +
  theme_classic() +
  xlab("Direction") +
  ylab("% Species within the top 10%") +
  scale_fill_manual(values=assoc.colours, name="Association") +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14))

assoc.plot = ggarrange(assoc.hist, assoc.bar, widths=c(1,0.3), common.legend=TRUE, align="h", labels=c("C", ""), font.label = list(size=18))

# final plot and save
ggarrange(order.plot, status.plot, assoc.plot, ncol=1)
ggsave(filename = "../figures/cag170_correlations.pdf", height=12, width=10)
