#!/usr/bin/env Rscript

# load libraries
library(tidyverse)
library(ALDEx2)
library(data.table)

# arguments to run in HPC
args <- commandArgs(trailingOnly=TRUE)

# define input / output and other variables
comparison <- args[1] #  All, or each of the 14 diseases
outfile <- paste0("aldex2_", comparison)

# load input file for transformed count reads
abund_data <- fread("/rds/project/rds-aFEMMKDjWlo/aalmeida/project_secretbugs/bwa/exp06/bwa_counts-filtered_samples.csv")

# load metadata
metadata <- fread("/rds/project/rds-aFEMMKDjWlo/rfs_data/metadata/metagenomes_03-2024.tsv")

# load input file to define taxa
infoSpecies <- fread("/rds/project/rds-aFEMMKDjWlo/rfs_data/metadata/species_uhgg_v1.2.tsv")

# clean dataframe 
abund_data <- column_to_rownames(abund_data, var = "Genome")

# make data frame useful by cleaning up a bit
metadata1 <- metadata %>% rename_with( ~ (gsub(".", "_", .x, fixed=T))) #%>% # I prefer underscores instead of points

# subset your metadata according to the disease that is input as args[1]
which_studies <- metadata1 %>% filter(Disease_name==comparison) %>% distinct(Study)
  
# filter metadata for those studies
metadata_disease <- metadata1 %>% filter(Study %in% which_studies$Study)
  
# now filter so that only the specific disease + healthy are included
metadata2 <- metadata_disease %>% filter(!is.na(Disease_name) &!is.na(Age_group) & !is.na(Continent)) %>%
  filter(Disease_name==ifelse((Health_state=="Diseased" & Disease_name!=comparison),NA,Disease_name))

# exclude studies that don't have health+disease
disease_studies = unique(metadata2[which(metadata2$Disease_name == comparison),"Study"])
healthy_studies = unique(metadata2[which(metadata2$Disease_name == "Healthy"),"Study"])
studies = intersect(disease_studies, healthy_studies)

# filter studies that have both healthy + disease samples
metadata3 <- metadata2 %>% filter(Study %in% studies$Study)

# update row names
metadata_ready <- metadata3 %>% mutate(Run1=Run) %>% column_to_rownames(var = "Run1")

# deduplicate samples
read_counts_samples = aggregate(Read_count ~ Sample, data=metadata_ready, FUN=sum)
rownames(read_counts_samples) = read_counts_samples$Sample
metadata_ready = unique(metadata_ready[,which(!colnames(metadata_ready) %in% c("Run", "Read_count", "Timepoint"))])
rownames(metadata_ready) = metadata_ready$Sample
metadata_ready$Read_count = read_counts_samples[rownames(metadata_ready),"Read_count"]

# filter abundance data using the metadata
abund_data_filtered <- abund_data[,metadata_ready$Sample]

# filter data based on prevalence 
thresh = 0.01 # prevalence_value #ratio 
filtered_genomes <- names(which(rowSums(abund_data_filtered > 0)/ncol(abund_data_filtered) > thresh))
abund_data_ready <- abund_data_filtered[filtered_genomes,]

# define what variables are going to be used as references:
age.group_ref <- dplyr::count(metadata_ready, Age_group) %>% arrange(desc(n)) %>% pull(Age_group) %>% first()
continent_ref <- dplyr::count(metadata_ready, Continent) %>% arrange(desc(n)) %>% pull(Continent) %>% first()
study_ref <- dplyr::count(metadata_ready, Study) %>% arrange(desc(n)) %>% pull(Study) %>% first()

# and print it for your records
cat(paste("We're using as reference the values of", age.group_ref, continent_ref, study_ref, "\n"))

# setup model matrix:
metadata_ready$Disease_name = relevel(factor(metadata_ready$Disease_name), "Healthy")
metadata_ready$Age_group = relevel(factor(metadata_ready$Age_group), age.group_ref)
metadata_ready$Continent = relevel(factor(metadata_ready$Continent), continent_ref)
metadata_ready$Study = relevel(factor(metadata_ready$Study), study_ref)

# define all variables of interest for the matrix
variables <- c("Disease_name", "Age_group", "Continent", "Study")
new_variables <- variables
for (v in variables) {
  if (length(levels(metadata_ready[[v]])) <= 1 & v!="Disease_name") { # check if they have more than one value
    new_variables <- new_variables[new_variables != v] # update the variables
  }
}

# transform those variables in the matrix expression
new_variables = c(new_variables, "Read_count")
formula_str <- paste0("~ ", paste0(new_variables, collapse = " + "))

# and finally create the model matrix
mm = model.matrix(as.formula(formula_str), metadata_ready[colnames(abund_data_ready),])

# run ALDEx2 
cat("Running ALDEx2 with", ncol(abund_data_ready), "samples and", nrow(abund_data_ready), "taxa ...")
clr.values = aldex.clr(abund_data_ready, mm, mc.samples = 128, denom = "all", useMC = TRUE) # genomes in rows, samples in columns 
glm.test = aldex.glm(clr.values, verbose = TRUE)
  
# make the column name with target disease
res_column <- paste0("Disease_name",comparison)

# filter results
all_results <- glm.test[,c(paste(res_column,":Est", sep=""), paste(res_column, ":pval", sep=""))]
  
# apply FDR correction
all_results$FDR <- p.adjust(all_results[,2], method="fdr")

# save table
results_final <- cbind(rownames(all_results), all_results)
colnames(results_final)[1] = "feature"
write.table(results_final, file = paste0(outfile, ".tsv"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep="\t")
  
