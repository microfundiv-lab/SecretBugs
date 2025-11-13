# Load libraries
library(lme4)
library(dplyr)
library(reshape2)
library(parallel)
library(ggplot2)

# Set working directory and load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/timeseries")
df.all = read.delim("cag170_genes.tsv", header=FALSE)
df.all$Gene = 1
df = acast(df.all, V1 ~ V2, fun.aggregate = sum, value.var = "Gene")
metadata = read.delim("../metadata/hmp2_ibd/metagenomes_dysb-score.tsv")
rownames(metadata) = metadata$Run

# Filter subjects
subjects = names(which(table(metadata$Subject) > 1))
metadata = metadata[which(metadata$Subject %in% subjects),]

# Subset to common samples
common.samples = intersect(rownames(df), rownames(metadata))
df = as.data.frame(df[common.samples, ])
df[df > 0] = 1
rownames(df) = metadata[rownames(df), "Sample"]

# filter based on prevalence
samples = names(which(rowSums(df) > 0))
df = df[samples,]
genes = names(which(colSums(df) > nrow(df)*0.01))

# fit glmm
fit_glmm = function(gene) {
  df$gene_status = df[[gene]]
  
  pres_dysb = mean(df[which(df[, gene] > 0),"Score"])
  abs_dysb = mean(df[which(df[, gene] == 0),"Score"])
  model = glmer(gene_status ~ Score + (1 | Subject), data = df, family = binomial)
  coef_summary = summary(model)$coefficients
  
  estimate = coef_summary["Score", "Estimate"]
  std_error = coef_summary["Score", "Std. Error"]
  pval = coef_summary["Score", "Pr(>|z|)"]
  
  return(c(pres_dysb, abs_dysb, estimate, std_error, pval))
}

# Set number of cores to use (adjust if needed)
n_cores = detectCores() - 1

# Run models in parallel
rownames(metadata) = metadata$Sample
df$Score = metadata[rownames(df), "Dysb_score"]
df$Subject = metadata[rownames(df), "Subject"]
glm.summary = mclapply(genes, fit_glmm, mc.cores = n_cores)

# Process output
raw.output = data.frame(t(data.frame(glm.summary)))
colnames(raw.output) = c("Dysb_pres", "Dysb_abs", "Estimate", "Std.error", "P-value")
rownames(raw.output) = genes

# Adjust p-values and add log-transformed p-values
raw.output$FDR = p.adjust(raw.output$`P-value`, method = "fdr")

# Filter significant results
sign.output = raw.output[which(raw.output$FDR < 0.05 & raw.output$`P-value` != 0), ]
sign.output = sign.output[order(sign.output$FDR), ]
sign.output$Sanity_check = ifelse(sign.output$Dysb_pres < sign.output$Dysb_abs, "Pass", "Fail")
sign.output = cbind(Gene = rownames(sign.output), sign.output)

# Save results
write.table(sign.output, file = "sign_genes_glmer.tsv", quote = FALSE, row.names = FALSE, sep = "\t")