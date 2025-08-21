# load libraries
library(data.table)
library(ggplot2)
library(vegan)
library(lme4)
library(lmerTest)
library(ggpubr)
library(ggrastr)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/alpha_div")
both =  list.files(pattern = "*_both.csv", recursive = TRUE, full.names=TRUE)
cultured = list.files(pattern = "*_cultured.csv", recursive = TRUE, full.names=TRUE)
uncultured = list.files(pattern = "*_uncultured.csv", recursive = TRUE, full.names=TRUE)

# function to generate alpha div dataframe
gen_alphadiv = function(x) {
  div.list = lapply(x, function(i) {
    analysis = strsplit(i, "/")[[1]][2]
    disease = strsplit(analysis, "_")[[1]][1]
    df = read.csv(i)
    df$Shannon = diversity(df[,!colnames(df) %in% c("Sample", "Variable", "Study")], index = "shannon")
    df$InvSimp = diversity(df[,!colnames(df) %in% c("Sample", "Variable", "Study")], index = "invsimpson")
    df$Disease = disease
    df$Variable = gsub(disease, "Diseased", df$Variable)
    return(df[,c("Sample", "Variable", "Shannon", "InvSimp", "Disease")])
  })
  
  # merge data
  div.combined = as.data.frame(rbindlist(div.list))
  
  # add read counts
  metadata = read.delim("../metadata/metagenomes_03-2024_samples.tsv")
  div.combined$Read.count = metadata[match(div.combined$Sample, metadata$Sample),"Read.count"]
  div.combined$Study = metadata[match(div.combined$Sample, metadata$Sample),"Study"]
  div.combined$Subject = metadata[match(div.combined$Sample, metadata$Sample),"Subject"]
  
  # rename diseases
  div.combined$Disease = gsub("CRC", "Colorectal cancer", div.combined$Disease)
  div.combined$Disease = gsub("CD", "Crohn's disease", div.combined$Disease)
  div.combined$Disease = gsub("UC", "Ulcerative colitis", div.combined$Disease)
  div.combined$Disease = gsub("T1D", "Type 1 diabetes", div.combined$Disease)
  div.combined$Disease = gsub("T2D", "Type 2 diabetes", div.combined$Disease)
  div.combined$Disease = gsub("MS", "Multiple sclerosis", div.combined$Disease)
  div.combined$Disease = gsub("ME", "Myalgic encephalomyelitis", div.combined$Disease)
  div.combined$Disease = gsub("T2D", "Type 2 diabetes", div.combined$Disease)
  div.combined$Disease = gsub("AS", "Ankylosing spondylitis", div.combined$Disease)
  div.combined$Disease = gsub("RA", "Rheumatoid arthritis", div.combined$Disease)
  div.combined$Disease = gsub("Parkinson", "Parkinson's disease", div.combined$Disease)
  div.combined$Disease = gsub("BoneDisease", "Bone disease", div.combined$Disease)
  div.combined$Variable = factor(div.combined$Variable, levels=c("Healthy", "Diseased"))
  return(div.combined)
}

# generate dataframes
div_both = gen_alphadiv(both)
div_cultured = gen_alphadiv(cultured)
div_cultured$Type = "Cultured"
div_uncultured = gen_alphadiv(uncultured)
div_uncultured$Type = "Uncultured"
div_long = div_cultured
div_long$Shannon.uncult = div_uncultured$Shannon
div_long$Ratio = log10(div_long$Shannon.uncult/div_long$Shannon)
div_long_unique = unique(div_long[,c("Sample", "Variable", "Shannon", "Shannon.uncult")])
div_melt = rbind(div_cultured, div_uncultured)

# add disease counts
source("../../scripts/alex/metadata_disease-numbers.R")
disease.counts = rowSums(disease.df[,c("Healthy", "Diseased")])
div_uncultured$Counts = disease.counts[div_uncultured$Disease]
div_uncultured$Label = paste(div_uncultured$Disease, " (n = ", div_uncultured$Counts, ")", sep="")

# scatterplot uncult vs cult
div.scatter = ggplot(div_long_unique, aes(x=Shannon, y = Shannon.uncult, colour=Variable)) +
  geom_point_rast(alpha=0.5) +
  geom_smooth(method = "lm", colour="grey20", linetype="dashed") +
  theme_classic() +
  scale_colour_manual(values=c("steelblue", "tomato"), name="Health state") +
  ylab("Shannon diversity (uncultured)") +
  xlab("Shannon diversity (cultured)") +
  theme(axis.title = element_text(size=18)) +
  theme(axis.text = element_text(size=16)) +
  theme(legend.text = element_text(size=16)) +
  theme(legend.title= element_text(size=16))

# uncult diversity per disease
div.box = ggplot(div_uncultured, aes(x=reorder(Label, -Shannon), y = Shannon, fill=Variable)) +
  geom_boxplot(alpha=0.6, outlier.colour=NA, width=0.7) +
  geom_point_rast(alpha=0.05, size=0.3, position = position_jitterdodge(jitter.width=0.2)) + 
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values=c("steelblue", "tomato"), name="Health state") +
  ylab("Shannon diversity\n(uncultured)") +
  theme(axis.title.y = element_blank()) +
  theme(axis.title.x = element_text(size=18)) +
  theme(axis.text = element_text(size=16)) +
  theme(legend.text = element_text(size=16)) +
  theme(legend.title= element_text(size=16))

# combine plots
ggarrange(div.box, div.scatter, ncol=2, common.legend=TRUE, widths=c(1,0.8), align="h", labels=c("a", "b"), font.label=list(size=22), legend = "bottom")
ggsave(file="figures/alpha_diversity.pdf", height=8, width=13)

# perform stats tests
corr.stats = summary(lm(div_long_unique$Shannon.uncult ~ div_long_unique$Shannon))

alpha_sign = function(div.combined) {
  disease.test = lapply(unique(div.combined$Disease), function(x) {
    sub_data = div.combined[div.combined$Disease == x, ]
    studies = length(unique(sub_data$Study))
    subjects = max(table(sub_data$Subject))
    sub_data$Shannon_rank = rank(sub_data$Shannon)
    
    if (subjects > 1) {
      # Use mixed-effects model with Subject as a random effect
      if (studies > 1) {
        model = lmer(Shannon_rank ~ Variable + Study + log10(Read.count) + (1 | Subject), data = sub_data)
      } else {
        model = lmer(Shannon_rank ~ Variable + log10(Read.count) + (1 | Subject), data = sub_data)
      }
      p_val = summary(model)$coefficients["VariableDiseased", "Pr(>|t|)"]
    } else {
      # Use linear model without random effect
      if (studies > 1) {
        model = lm(Shannon_rank ~ Variable + Study + log10(Read.count), data = sub_data)
      } else {
        model = lm(Shannon_rank ~ Variable + log10(Read.count), data = sub_data)
      }
      p_val = summary(model)$coefficients["VariableDiseased", "Pr(>|t|)"]
    }
    
    return(p_val)
  })
  
  # Format results
  disease.test = data.frame(t(data.frame(disease.test)))
  rownames(disease.test) = unique(div.combined$Disease)
  disease.test$FDR = p.adjust(disease.test[,1])
  disease.test$Sign = ifelse(disease.test$FDR < 0.05, "Significant", "NS")
  return(disease.test)
}

sign_both = alpha_sign(div_both)
sign_cultured = alpha_sign(div_cultured)
sign_uncultured = alpha_sign(div_uncultured)