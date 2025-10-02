# load libraries
library(ggplot2)
library(ggridges)
library(ggrastr)
library(ggpubr)
library(data.table)
library(scales)
library(tidyverse)
library(RColorBrewer)
library(vegan)
library(lme4)
library(lmerTest)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/")
species.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
cultured = species.metadata[which(species.metadata$Status == "Cultured"),"Species_rep"]
uncultured = species.metadata[which(species.metadata$Status == "Uncultured"),"Species_rep"]
metadata = read.delim("metadata/metagenomes_03-2024_samples.tsv")
rownames(metadata) = metadata$Sample

# load counts and convert to relab
counts.filt = as.data.frame(fread("bwa/bwa_counts-filtered_samples.csv"))
rownames(counts.filt) = counts.filt$Genome
counts.filt = counts.filt[,-1]
source("../scripts/alex/metadata_disease-numbers.R")
counts.filt = counts.filt[species.metadata$Genome,rownames(metadata.disease)]
counts.norm = counts.filt/species.metadata$Length
counts.relab = t(t(counts.norm)/colSums(counts.norm)*100)

# calculate alpha div
div.cult = as.data.frame(diversity(t(counts.relab[cultured,]), index="shannon"))
cult.rscl = t(t(counts.relab[cultured,])/colSums(counts.relab[cultured,]))*100
div.cult.rscl = as.data.frame(diversity(t(cult.rscl), index="shannon"))
div.uncult = as.data.frame(diversity(t(counts.relab[uncultured,]), index="shannon"))
uncult.rscl = t(t(counts.relab[uncultured,])/colSums(counts.relab[uncultured,]))*100
div.uncult.rscl = as.data.frame(diversity(t(uncult.rscl), index="shannon"))

# list of diseases to analyse
diseases = c("Adenoma", "AS", "BoneDisease", "CD", "CRC","ME",
             "MS", "Obesity", "Parkinson", "RA", "T1D", "T2D", "UC")

# add metadata to diversity dfs
add_meta = function(x) {
  disease.data = list()
  for (d in diseases) {
    # check in which studies that disease is and if it also has healthy samples
    disease.studies = unique(metadata[which(metadata$Disease.name == d),"Study"])
    healthy.studies = unique(metadata[which(metadata$Disease.name == "Healthy"),"Study"])
    which_studies = intersect(disease.studies, healthy.studies)
    
    # filter metadata for those studies
    metadata.filt = metadata[which(metadata$Study %in% which_studies),]
    
    # filter so that only the specific disease + healthy are included
    metadata.filt = metadata.filt %>% filter(!is.na(Disease.name)) %>%
      filter(Disease.name == ifelse((Health.state=="Diseased" & Disease.name != d), NA, Disease.name))
    
    rates.filt = x[metadata.filt$Sample,, drop=FALSE]
    colnames(rates.filt)[1] = "Shannon"
    rates.filt$Disease = d
    rates.filt$Health.state = metadata.filt$Health.state
    rates.filt$Study = metadata.filt$Study
    rates.filt$Subject = metadata.filt$Subject
    rates.filt$Read.count = metadata.filt$Read.count
    disease.data[[d]] = rates.filt
  }
  
  # combine data frames
  all.data = do.call(what = rbind, args = disease.data)
  return(all.data)
}

div.cult.fi = add_meta(div.cult)
div.cult.rscl.fi = add_meta(div.cult.rscl)
div.uncult.fi = add_meta(div.uncult)
div.uncult.rscl.fi = add_meta(div.uncult.rscl)

# perform stats tests
corr.stats = summary(lm(div.cult[,1] ~ div.uncult[,1]))
corr.stats.rscl = summary(lm(div.cult.rscl[,1] ~ div.uncult.rscl[,1]))

alpha_sign = function(div.combined) {
  disease.test = lapply(unique(div.combined$Disease), function(x) {
    sub_data = div.combined[div.combined$Disease == x, ]
    studies = length(unique(sub_data$Study))
    subjects = max(table(sub_data$Subject))
    
    if (subjects > 1) {
      # Use mixed-effects model with Subject as a random effect
      if (studies > 1) {
        model = lmer(Shannon ~ Health.state + Study + log10(Read.count) + (1 | Subject), data = sub_data)
      } else {
        model = lmer(Shannon ~ Health.state + log10(Read.count) + (1 | Subject), data = sub_data)
      }
      p_val = summary(model)$coefficients["Health.stateHealthy", "Pr(>|t|)"]
    } else {
      # Use linear model without random effect
      if (studies > 1) {
        model = lm(Shannon ~ Health.state + Study + log10(Read.count), data = sub_data)
      } else {
        model = lm(Shannon ~ Health.state + log10(Read.count), data = sub_data)
      }
      p_val = summary(model)$coefficients["Health.stateHealthy", "Pr(>|t|)"]
    }
    
    return(p_val)
  })
  
  # Format results
  disease.test = data.frame(t(data.frame(disease.test)))
  rownames(disease.test) = unique(div.combined$Disease)
  disease.test$FDR = p.adjust(disease.test[,1], method="fdr")
  disease.test$Sign = ifelse(disease.test$FDR < 0.05, "Significant", "NS")
  return(disease.test)
}

sign_cult_relab = alpha_sign(div.cult.fi)
sign_cult_relab_rscl = alpha_sign(div.cult.rscl.fi)
sign_uncult_relab = alpha_sign(div.uncult.fi)
sign_uncult_relab_rscl = alpha_sign(div.uncult.rscl.fi)