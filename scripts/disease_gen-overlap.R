# load libraries
library(stringr)
library(tidyr)
library(data.table)
library(matrixStats)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/")
all.files = list.files("aldex2/", pattern=".tsv")
variables = str_match(all.files, "aldex2_\\s*(.*?)\\s*.tsv")
              
for (variable in variables[,2]) {
  aldex2.file = list.files(path = "aldex2/", pattern = variable, full.names=TRUE)
  maaslin2.file = list.files(path = "maaslin2/", pattern = variable, full.names=TRUE)
  
  # filter data
  maaslin2.all = read.delim(maaslin2.file, stringsAsFactors = FALSE)
  maaslin2.all = maaslin2.all[which(maaslin2.all$metadata == "Disease_name"),]
  maaslin2.all$qval = p.adjust(maaslin2.all$pval, method="fdr")
  maaslin2.sign = maaslin2.all[which(maaslin2.all$qval < 0.05),]
  
  aldex2.all = read.delim(aldex2.file, stringsAsFactors = FALSE)
  aldex2.sign = aldex2.all[which(aldex2.all$FDR < 0.05),]
  overlap = merge(maaslin2.sign, aldex2.sign, by="feature", all=FALSE)
  overlap = overlap[which(sign(overlap$coef) == sign(overlap[,colnames(aldex2.sign)[2]])),]
  
  # save table
  write.table(overlap, file=paste0("meta-analysis/", variable, ".tsv"), sep="\t", row.names=FALSE, quote=FALSE)
}
