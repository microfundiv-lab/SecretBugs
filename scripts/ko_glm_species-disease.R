# load libraries
library(reshape2)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# load data 
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")
disease.summ = read.csv("meta-analysis/species_ds.csv")
rownames(disease.summ) = disease.summ$Species_rep
disease.summ = disease.summ[,-1]

# filter out species with mixed signals
health.species = names(which(rowSums(disease.summ[,2:9] < 0) > 0))
disease.species = names(which(rowSums(disease.summ[,2:9] > 0) > 0))
mixed.signal = intersect(health.species, disease.species)
disease.summ = disease.summ[which(disease.summ$DS != 0 & !disease.summ$Species_rep %in% mixed.signal),]
disease.summ$Classification = ifelse(disease.summ$DS > 0, "Disease", "Health")

# parse KEGG cateories
ko.classes.ori = read.delim("genofan/KO_Orthology_ko00001.txt", header=FALSE)
ko.classes = separate(ko.classes.ori, V4, "V5", sep=" ")
ko.desc = cbind(ko.classes.ori, ko.classes)
ko.desc = unique(ko.desc[,c("V5", "V4")])
rownames(ko.desc) = ko.desc$V5

# load kegg data
keggor.df = read.delim("genofan/kegg_orthologs.tsv", sep="", header = F)
keggor.df$V1 = gsub("_\\d+$", "", keggor.df$V1)
keggor.matrix = data.frame(acast(keggor.df, V2 ~ V1, length))
keggor.matrix = keggor.matrix[,rownames(disease.summ)]
keggor.matrix = t(keggor.matrix[which(rowSums(keggor.matrix > 0) > ncol(keggor.matrix)*0.01),])
keggor.matrix[keggor.matrix > 0] = 1
kegg.features = colnames(keggor.matrix)

# prepare data
final.df = merge(keggor.matrix, disease.summ, by="row.names")
rownames(final.df) = final.df$Row.names
final.df = final.df[,-1]
final.df = final.df[order(final.df$Classification),]

# add genome type
species.metadata = unique(read.delim("metadata/genomes_uhgg-v1.2.tsv")[,c("Genome", "Genome_type")])
rownames(species.metadata) = species.metadata$Genome
final.df$Gtype = species.metadata[rownames(final.df), "Genome_type"]

# run glm
glm.summary = lapply(kegg.features, function(x) {
  cat(paste("Running glm for", x, "\n", sep=" "))
  glm.formula = formula(paste(paste0(x, " ~ Classification + Order + Gtype")))
  glm.out = summary(glm(glm.formula, data = final.df, family = "binomial", maxit = 1000))
  return(glm.out$coefficients["ClassificationHealth",-3])
})

# process output
raw.output = data.frame(t(data.frame(glm.summary)))
rownames(raw.output) = kegg.features
raw.output$FDR = p.adjust(raw.output[,3])
sign.output = raw.output[which(raw.output$FDR < 0.05),]
sign.output = sign.output[order(sign.output$FDR),]
sign.output$Estimate = -sign.output$Estimate
sign.fi = cbind(rownames(sign.output), sign.output)
colnames(sign.fi) = c("KEGG_Ortholog", "Estimate", "Std.Error", "Pvalue", "FDR")
sign.fi$Description = ko.desc[sign.fi$KEGG_Ortholog,"V4"]
write.table(sign.fi, file="genofan/ko-glm_results.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# select top ko and parse description
ko.pos = rownames(sign.fi)[which(sign.fi$Estimate > 0)]
ko.pos.top = sign.fi[order(sign.fi$Estimate, decreasing=TRUE)[1:20],]
ko.neg = rownames(sign.fi)[which(sign.fi$Estimate < 0)]
ko.neg.top = sign.fi[order(sign.fi$Estimate, decreasing=FALSE)[1:20],]
ko.top = rbind(ko.pos.top, ko.neg.top)
ko.top = ko.top[order(ko.top$Estimate, decreasing=TRUE),]
ko.top$Description = gsub(".*;","",ko.top$Description)
ko.top$Description = gsub("\\[EC.*","",ko.top$Description)
ko.top$Description = gsub("\\(.*","",ko.top$Description)
ko.top$Description = str_to_title(ko.top$Description)
ko.top$Description = gsub("16s","16S",ko.top$Description)
ko.top$Description = gsub("30s","30S",ko.top$Description)
ko.top$Description = gsub("Trna","tRNA",ko.top$Description)
ko.top$Description = gsub("Rna","RNA",ko.top$Description)
ko.top$Description = gsub("Mfs","MFS",ko.top$Description)
ko.top$Description = gsub("Udp","UDP",ko.top$Description)
ko.top$Description = gsub("Rrna","rRNA",ko.top$Description)
ko.top$Description = gsub("Topoisomerase Iv","Topoisomerase IV",ko.top$Description)
ko.top$Description = gsub("Dnai","DnaI",ko.top$Description)
ko.top$Description = gsub("Nss","NSS",ko.top$Description)
ko.top$Description = gsub("Hth","HTH",ko.top$Description)
ko.top$Description = gsub("Involved In","involved in",ko.top$Description)
ko.top$Description = gsub("Gtpase","GTPase",ko.top$Description)
ko.top$Description = gsub("Htpg","HtpG",ko.top$Description)
ko.top$Description = gsub("Dna ","DNA ",ko.top$Description)
ko.top$Description = sub("^\\s+", "", ko.top$Description)

# heatmap colors
neg = rownames(final.df)[which(final.df$Classification == "Health")]
pos = rownames(final.df)[which(final.df$Classification == "Disease")]
annot = final.df[c(neg, pos), c("Classification", "Status", "Order"), drop=FALSE]
class.colors = c("steelblue", "tomato")
names(class.colors) = c("Health", "Disease")
status.colors = c("steelblue1", "palegreen2")
names(status.colors) = c("Cultured", "Uncultured")
top.taxa = names(sort(table(annot$Order), decreasing = TRUE)[1:5])
annot$Order = ifelse(annot$Order %in% top.taxa, annot$Order, "Other")
top.taxa.colors = c(brewer.pal(length(top.taxa), "Set3"), "lightgrey")
names(top.taxa.colors) = c(top.taxa, "Other")
annot.colors = list(Classification=class.colors, Status=status.colors, Order=top.taxa.colors)
kegg.annot = data.frame(ifelse(ko.top$Estimate > 0, "Disease", "Health"))
rownames(kegg.annot) = rownames(ko.top)
colnames(kegg.annot) = "Classification"

# heatmap plot
sign.df = t(final.df[c(neg, pos), ko.top$KEGG_Ortholog])
pheatmap(sign.df, show_colnames = FALSE, show_rownames = TRUE, cluster_rows = TRUE,
         annotation_col = annot, annotation_colors = annot.colors,
         annotation_names_row = FALSE, annotation_names_col=TRUE,
         color=c("lightgrey", "darkgreen"), labels_row = ko.top$Description,
         filename="figures/ko-heatmap_top20.pdf", width = 15, height=7)
