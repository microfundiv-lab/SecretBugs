# load libraries
library(reshape2)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(lme4)
library(ggpubr)

# Optional parallel (auto-fallback to lapply if not installed)
future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
glmer_apply = function(X, FUN) future.apply::future_lapply(X, FUN, future.seed = TRUE)
options(contrasts = c("contr.treatment", "contr.poly"))

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")
top1.pos = scan("fastspar/top-pos-1_cag-170.txt", what="")
genome.metadata  = read.delim("metadata/genomes-nr_metadata.tsv")
cobrapy.res = read.delim("genofan/cag170_ecology/cobrapy_M3.tsv", header=FALSE)
cobrapy.uptake = cobrapy.res[which(cobrapy.res$V4 == "Uptake"),]
cobrapy.secretion = cobrapy.res[which(cobrapy.res$V4 == "Secretion"),]
cobrapy.dict = unique(rbind(cobrapy.secretion[,c("V2", "V3")], cobrapy.uptake[,c("V2", "V3")]))
rownames(cobrapy.dict) = cobrapy.dict$V2

# parse taxonomy
cag170.species = unique(genome.metadata[grep("g__CAG-170", genome.metadata$Lineage),"Species_rep"])
genome.selected = genome.metadata[genome.metadata$Species_rep %in% c(cag170.species, top1.pos), c("Genome","Species_rep", "Completeness")]
rownames(genome.selected) = genome.selected$Genome
selected = intersect(genome.selected$Genome, unique(cobrapy.res$V1))
genome.selected = genome.selected[selected,]
genome.selected$Classification = NA
genome.selected$Classification = ifelse(genome.selected$Species_rep %in% cag170.species, "CAG-170", genome.selected$Classification)
genome.selected$Classification = ifelse(genome.selected$Species_rep %in% top1.pos, "Co-colonizers", genome.selected$Classification)

# prepare model
cobrapy_analysis = function(x) {
  cobrapy.matrix = data.frame(acast(x, V2 ~ V1, length)) # secreted
  cobrapy.matrix = cobrapy.matrix[,selected]
  cobrapy.matrix = t(cobrapy.matrix[rowSums(cobrapy.matrix > 0) > ncol(cobrapy.matrix) * 0.01, , drop = FALSE])
  cobrapy.matrix[cobrapy.matrix > 0] = 1L
  storage.mode(cobrapy.matrix) = "integer"
  features = colnames(cobrapy.matrix)
  
  glmer.df = merge(cobrapy.matrix, genome.selected, by = "row.names")
  rownames(glmer.df) = glmer.df$Row.names
  glmer.df$Row.names = NULL
  
  # ensure factors are well-defined
  glmer.df$Classification = droplevels(factor(glmer.df$Classification))
  glmer.df$Species_rep  = factor(glmer.df$Species_rep)
  
  # glmer function
  fit_one = function(x) {
    
    # define new df
    dat = glmer.df
    y = dat[[x]]
    
    # Basic validity checks
    if (!all(y %in% c(0L, 1L)) || length(unique(y)) < 2L) {
      return(c(Estimate = NA_real_, Std.Error = NA_real_, Pvalue = NA_real_))
    }
    dat$resp = as.numeric(y)
    
    # Levels after dropping rows
    dat$Classification = droplevels(dat$Classification)
    dat$Species_rep = droplevels(dat$Species_rep)
    
    # check there is variation within the target genus; perfect separation will otherwise break the model
    tab_tg = table(dat$resp[dat$Classification == "CAG-170"])
    if (length(tab_tg) == 1L) {
      
      # all 0s or all 1s within target genus -> treat as non-estimable
      return(c(Estimate = NA_real_, Std.Error = NA_real_, Pvalue = NA_real_))
    }
    
    # fit glmer safely
    fit = tryCatch(
      suppressWarnings(
        glmer(resp ~ Classification + (1 | Species_rep),
              data = dat,
              family = binomial,
              nAGQ = 0,  # faster (Laplace would be nAGQ=1)
              control = glmerControl(optimizer = "bobyqa",
                                     calc.derivs = FALSE,
                                     optCtrl = list(maxfun = 1e5)))),
      error = function(e) NULL)
    if (is.null(fit)) {
      return(c(Estimate = NA_real_, Std.Error = NA_real_, Pvalue = NA_real_))
    }
    
    # convergence / singularity checks
    conv = fit@optinfo$conv
    conv_msgs = c(conv$lme4$messages, conv$messages)
    opt_code  = tryCatch(conv$opt, error = function(e) NULL)
    
    if (!is.null(conv_msgs) && length(conv_msgs) > 0L) {
      return(c(Estimate = NA_real_, Std.Error = NA_real_, Pvalue = NA_real_))
    }
    if (!is.null(opt_code) && !identical(opt_code, 0L)) {
      return(c(Estimate = NA_real_, Std.Error = NA_real_, Pvalue = NA_real_))
    }
    if (isTRUE(isSingular(fit, tol = 1e-4))) {
      return(c(Estimate = NA_real_, Std.Error = NA_real_, Pvalue = NA_real_))
    }
    
    # extract coefficient for target genus
    sm = summary(fit)
    cf = sm$coefficients
    row_name = "ClassificationCo-colonizers"
    if (!(row_name %in% rownames(cf))) {
      return(c(Estimate = NA_real_, Std.Error = NA_real_, Pvalue = NA_real_))
    }
    
    out = cf[row_name, c("Estimate", "Std. Error", "Pr(>|z|)")]
    names(out) = c("Estimate", "Std.Error", "Pvalue")
    as.numeric(out)
  }
  
  # run model in parallel
  glmer.summary = glmer_apply(features, fit_one)
  
  # collect results
  raw.output = do.call(rbind, glmer.summary) |> as.data.frame()
  rownames(raw.output) = features
  colnames(raw.output) = c("Estimate", "Std.Error", "Pvalue")
  raw.output = raw.output[which(!is.na(raw.output$Pvalue)),]
  raw.output$FDR = p.adjust(raw.output$Pvalue, method = "BH")
  raw.output$Result = ifelse(raw.output$FDR < 0.05, ifelse(raw.output$Estimate > 0, "Enriched in co-colonizers", "Enriched in CAG-170"), "Non-significant")
  sign.output = raw.output[which(raw.output$FDR < 0.05),]
  return(sign.output)
}

# results
uptake.results = cobrapy_analysis(cobrapy.uptake)
secretion.results = cobrapy_analysis(cobrapy.secretion)
common = intersect(rownames(uptake.results), rownames(secretion.results))

# prepare colors
annot = genome.selected[, "Classification", drop=FALSE]
class.colors = c("palegreen3", "steelblue")
names(class.colors) = c("CAG-170", "Co-colonizers")

cobrapy.uptake.color.df = unique(cobrapy.uptake[,c("V2", "V4")])
rownames(cobrapy.uptake.color.df) = paste0("Uptake_", cobrapy.uptake.color.df$V2)
cobrapy.uptake.color.df = cobrapy.uptake.color.df[,-1, drop=FALSE]

cobrapy.secretion.color.df = unique(cobrapy.secretion[,c("V2", "V4")])
rownames(cobrapy.secretion.color.df) = paste0("Secretion_", cobrapy.secretion.color.df$V2)
cobrapy.secretion.color.df = cobrapy.secretion.color.df[,-1, drop=FALSE]

cobrapy.color.df = rbind(cobrapy.secretion.color.df, cobrapy.uptake.color.df)
colnames(cobrapy.color.df) = "Type"
metab.colors = c("#d2c3ff", "grey40")
names(metab.colors) = c("Secretion", "Uptake")

heat.colors = list(Type=metab.colors, Classification=class.colors)

# prepare heatmap
cobrapy.uptake.res = data.frame(acast(cobrapy.uptake, V2 ~ V1, length))
cobrapy.uptake.res = t(cobrapy.uptake.res[rownames(uptake.results), selected])
colnames(cobrapy.uptake.res) = paste0("Uptake_", colnames(cobrapy.uptake.res))

cobrapy.secretion.res = data.frame(acast(cobrapy.secretion, V2 ~ V1, length))
cobrapy.secretion.res[cobrapy.secretion.res > 0] = 1
cobrapy.secretion.res = t(cobrapy.secretion.res[rownames(secretion.results), selected])
colnames(cobrapy.secretion.res) = paste0("Secretion_", colnames(cobrapy.secretion.res))

cobrapy.combined = merge(cobrapy.secretion.res, cobrapy.uptake.res, by="row.names")
rownames(cobrapy.combined) = cobrapy.combined$Row.names
cobrapy.combined = cobrapy.combined[,-1]
col.labels = cobrapy.dict[gsub("Uptake_", "", gsub("Secretion_", "", colnames(cobrapy.combined))),"V3"]
pheatmap(t(cobrapy.combined), show_colnames = FALSE, color=c("#fdfde4", "#f3756b"),
         annotation_col = annot, annotation_row = cobrapy.color.df, annotation_colors = heat.colors,
         labels_row = col.labels, clustering_method = "ward.D",
         filename = "figures/cag170_metabolites_heatmap.pdf", width=11, height=7)

# prepare data for barplot
cobrapy.combined$Classification = genome.selected[rownames(cobrapy.combined), "Classification"]
cobrapy.merged = data.frame(t(aggregate(. ~ Classification, data=cobrapy.combined, FUN=sum)))
colnames(cobrapy.merged) = c("CAG-170", "Co-colonizers")
cobrapy.merged = cobrapy.merged[-1,]
cobrapy.merged$`CAG-170` = as.numeric(cobrapy.merged$`CAG-170`)
cobrapy.merged$`Co-colonizers` = as.numeric(cobrapy.merged$`Co-colonizers`)
cobrapy.merged$`CAG-170_prop` = c(cobrapy.merged$`CAG-170` / length(which(genome.selected$Classification == "CAG-170")))*100
cobrapy.merged$`Co-colonizers_prop` = c(cobrapy.merged$`Co-colonizers` / length(which(genome.selected$Classification == "Co-colonizers")))*100
cobrapy.merged$Metabolite = rownames(cobrapy.merged)
cobrapy.melt = reshape2::melt(cobrapy.merged)
cobrapy.melt = separate(Metabolite, into = c("type", "metabolite"), sep = "_", extra = "merge", fill = "right", data=cobrapy.melt)
cobrapy.melt$metabolite_name = cobrapy.dict[cobrapy.melt$metabolite,"V3"]

# secretion barplot
cobrapy.melt.sec = cobrapy.melt[which(grepl("_prop", cobrapy.melt$variable) & cobrapy.melt$type == "Secretion"),]
cobrapy.melt.sec.cag170 = cobrapy.melt.sec[which(cobrapy.melt.sec$variable == "CAG-170_prop"),]
order.sec = cobrapy.melt.sec.cag170[order(cobrapy.melt.sec.cag170$value, decreasing=TRUE),"metabolite_name"]
sec.plot = ggplot(cobrapy.melt.sec, aes(x=metabolite_name, y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge", alpha=0.75) +
  ylab("Proportion of genomes (%)") +
  xlab("") +
  theme_classic() +
  theme(axis.title = element_text(size=14)) +
  theme(axis.text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) +
  scale_x_discrete(limits=order.sec) +
  scale_fill_manual(values=c("palegreen3", "steelblue"), labels=c("CAG-170", "Co-colonizers"), name="Classification")

# uptake barplot
cobrapy.melt.up = cobrapy.melt[which(grepl("_prop", cobrapy.melt$variable) & cobrapy.melt$type == "Uptake"),]
cobrapy.melt.up.cag170 = cobrapy.melt.up[which(cobrapy.melt.up$variable == "CAG-170_prop"),]
order.up = cobrapy.melt.up.cag170[order(cobrapy.melt.up.cag170$value, decreasing=TRUE),"metabolite_name"]
up.plot = ggplot(cobrapy.melt.up, aes(x=metabolite_name, y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge", alpha=0.75) +
  ylab("Proportion of genomes (%)") +
  xlab("") +
  theme_classic() +
  theme(axis.title = element_text(size=14)) +
  theme(axis.text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1)) +
  scale_x_discrete(limits=order.up) +
  scale_fill_manual(values=c("palegreen3", "steelblue"), labels=c("CAG-170", "Co-colonizers"), name="")

# combine plots
ggarrange(sec.plot, up.plot, nrow=1, align="h", common.legend=TRUE, widths=c(1,0.3))
ggsave(filename = "figures/cag170_metabolites_sign.pdf", width=15, height=5)
