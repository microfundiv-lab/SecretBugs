# load libraries
library(reshape2)
library(pheatmap)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(lme4)

# Optional parallel (auto-fallback to lapply if not installed)
future::plan(future::multisession, workers = max(1, future::availableCores() - 1))
glmer_apply = function(X, FUN) future.apply::future_lapply(X, FUN, future.seed = TRUE)
options(contrasts = c("contr.treatment", "contr.poly"))

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")

species.metadata = read.delim("metadata/species_uhgg_v1.2.tsv")
species.oscillo  = species.metadata[grepl("f__Oscillospiraceae", species.metadata$Lineage), ]

genome.metadata  = read.delim("metadata/genomes-nr_metadata.tsv")

# parse taxonomy
ranks = c("Domain","Phylum","Class","Order","Family","Genus","Species")
species.oscillo  = separate(species.oscillo, col = Lineage, into = ranks, sep = ";")
rownames(species.oscillo) = species.oscillo$Genome

genome.oscillo = genome.metadata[genome.metadata$Species_rep %in% species.oscillo$Species_rep, c("Genome","Species_rep", "Completeness")]
genome.oscillo$Genus = species.oscillo[genome.oscillo$Species_rep, "Genus"]
genome.oscillo = genome.oscillo[genome.oscillo$Genus != "g__", ]
rownames(genome.oscillo) = genome.oscillo$Genome

# load kegg presence/absence data
keggor.df = read.delim("oscillo/kegg_orthologs.tsv", sep = "", header = FALSE)
keggor.df$V1 = gsub("_\\d+$", "", keggor.df$V1)
keggor.matrix = data.frame(acast(keggor.df, V2 ~ V1, length))

# filter dataset
selected.genomes = intersect(colnames(keggor.matrix), rownames(genome.oscillo))
genome.oscillo   = genome.oscillo[selected.genomes, , drop = FALSE]
keggor.matrix    = keggor.matrix[, selected.genomes, drop = FALSE]

# keep KO features present in >10% of genomes and binarize
keggor.matrix = t(keggor.matrix[rowSums(keggor.matrix > 0) > ncol(keggor.matrix) * 0.01, , drop = FALSE])
keggor.matrix[keggor.matrix > 0] = 1L
storage.mode(keggor.matrix) = "integer"
kegg.features = colnames(keggor.matrix)

# prepare model
genus.df = merge(keggor.matrix, genome.oscillo, by = "row.names")
rownames(genus.df) = genus.df$Row.names
genus.df$Row.names = NULL

# ensure factors are well-defined
genus.df$Genus        = droplevels(factor(genus.df$Genus))
genus.df$Species_rep  = factor(genus.df$Species_rep)

# function for glmm per feature
target_genus = "g__CAG-170"

fit_one = function(x) {
  
  # define new df
  dat = genus.df
  y = dat[[x]]
 
  # Basic validity checks
  if (!all(y %in% c(0L, 1L)) || length(unique(y)) < 2L) {
    return(c(Estimate = NA_real_, Std.Error = NA_real_, Pvalue = NA_real_))
  }
  dat$resp = as.numeric(y)
  
  # Levels after dropping rows
  dat$Genus       = droplevels(dat$Genus)
  dat$Species_rep = droplevels(dat$Species_rep)
  
  # check there is variation within the target genus; perfect separation will otherwise break the model
  tab_tg = table(dat$resp[dat$Genus == target_genus])
  if (length(tab_tg) == 1L) {
    
    # all 0s or all 1s within target genus -> treat as non-estimable
    return(c(Estimate = NA_real_, Std.Error = NA_real_, Pvalue = NA_real_))
  }
  
  # fit glmer safely
  fit = tryCatch(
    suppressWarnings(
      glmer(resp ~ Genus + (1 | Species_rep),
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
  row_name = paste0("Genus", target_genus)
  if (!(row_name %in% rownames(cf))) {
    return(c(Estimate = NA_real_, Std.Error = NA_real_, Pvalue = NA_real_))
  }
  
  out = cf[row_name, c("Estimate", "Std. Error", "Pr(>|z|)")]
  names(out) = c("Estimate", "Std.Error", "Pvalue")
  as.numeric(out)
}

# run model in parallel
glmer.summary = glmer_apply(kegg.features, fit_one)

# collect results
raw.output = do.call(rbind, glmer.summary) |> as.data.frame()
rownames(raw.output) = kegg.features
colnames(raw.output) = c("Estimate", "Std.Error", "Pvalue")
raw.output = raw.output[which(!is.na(raw.output$Pvalue)),]
raw.output$FDR = p.adjust(raw.output$Pvalue, method = "BH")
raw.output = cbind(KEGG_Ortholog = rownames(raw.output), raw.output)
raw.output$Result = ifelse(raw.output$FDR < 0.05, ifelse(raw.output$Estimate > 0, "Enriched in CAG-170", "Depleted in CAG-170"), "Non-significant")

# add kegg data
ko.classes.ori = read.delim("genofan/KO_Orthology_ko00001.txt", header = FALSE)
ko_map = transform(
  ko.classes.ori,
  KO   = sub("^\\s*(K\\d+).*", "\\1", V4),
  Desc = sub("^\\s*K\\d+\\s*", "", V4))
ko_map = unique(ko_map[, c("KO", "Desc")])
rownames(ko_map) = ko_map$KO
raw.output$Description = ko_map[raw.output$KEGG_Ortholog, "Desc"]

# find significant genes
sign.output = raw.output[!is.na(raw.output$FDR) & raw.output$FDR < 0.05, , drop = FALSE]
sign.output = sign.output[order(abs(sign.output$Estimate), decreasing=TRUE), , drop = FALSE]

# volcano plot
volcano.plot = ggplot(raw.output, aes(x = Estimate, y = -log(FDR), colour = Result)) +
  geom_point(alpha=0.8, size=1) +
  geom_vline(xintercept=0, linetype="dashed", colour="darkgrey") +
  geom_hline(yintercept=-log(0.05), linetype="dashed", linewidth=0.5) +
  xlab("log(Odds Ratio)") +
  ylab("-log(FDR)") +
  scale_colour_manual(values=c("tomato", "steelblue", "darkgrey")) +
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(axis.title = element_text(size=14)) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_blank())
ggsave(file="figures/oscillo_genes.pdf", heigh=5, width=6)

# save output
write.table(sign.output, file = "oscillo/ko-glmer_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(selected.genomes, file = "oscillo/selected_genomes.txt", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
