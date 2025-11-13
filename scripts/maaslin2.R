#!/usr/bin/env Rscript

# load libraries
library(tidyverse)
library(data.table)
library(Maaslin2)

args <- commandArgs(trailingOnly=TRUE)

# define input / output and other variables
comparison <- args[1] #  All, or each of the 14 diseases
cpus <- args[2] # in the cluster we can ask lots to make it faster!
outfile <- paste0("maaslin2_", comparison)

# load input file for transformed count reads
abund_data <- fread("/rds/project/rds-aFEMMKDjWlo/aalmeida/project_secretbugs/bwa/bwa_counts-filtered_samples.csv")

# load metadata
metadata <- fread("/rds/project/rds-aFEMMKDjWlo/rfs_data/metadata/metagenomes_03-2024.tsv")

# load input file to calculate genome lengths
infoSpecies <- fread("/rds/project/rds-aFEMMKDjWlo/rfs_data/metadata/species_uhgg_v1.2.tsv")

# filter genome lengths
genome_length <- infoSpecies %>% dplyr::select(Genome, Length)

# calculate relative abundance using transformed counts / genome length
relative_abund <- inner_join(abund_data,genome_length, by="Genome") %>% # joins both df
  column_to_rownames(var = "Genome") %>% # names the rows with Genome
  mutate(across(everything(), ~ .x/Length)) # divides Length by all columns of each row
  
relative_abund1 <- scale(relative_abund, center = FALSE, scale = colSums(relative_abund))
relative_abundance <- as.data.frame(relative_abund1) %>% dplyr::select(!Length) %>% # deletes the Length column
  mutate(across(everything(), ~ .x*100)) # makes it percentage

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
metadata_ready <- metadata3 %>% mutate(Run1=Run) %>% column_to_rownames(var = "Run1") %>%
  mutate(Disease_name = ifelse(Disease_name == "Healthy", "a_Healthy", Disease_name)) # see note below

# deduplicate samples
read_counts_samples = aggregate(Read_count ~ Sample, data=metadata_ready, FUN=sum)
rownames(read_counts_samples) = read_counts_samples$Sample
metadata_ready = unique(metadata_ready[,which(!colnames(metadata_ready) %in% c("Run", "Read_count", "Timepoint"))])
rownames(metadata_ready) = metadata_ready$Sample
metadata_ready$Read_count = read_counts_samples[rownames(metadata_ready),"Read_count"]

# filter abundance data using the metadata
abund_data_ready <- relative_abundance[,metadata_ready$Sample]
abund_data_ready[is.na(abund_data_ready)] = 0

# define the minimum percent of samples for which a feature is detected at minimum abundance:
prevalence <- 0.01 # here we're using most common value of 1%

# define what variables are going to be used as references:
age.group_ref <- (dplyr::count(metadata_ready,Age_group) %>% arrange(-n) %>% dplyr::select(Age_group))[[1]][1]
continent_ref <- (dplyr::count(metadata_ready,Continent) %>% arrange(-n) %>% dplyr::select(Continent))[[1]][1]
study_ref <- (dplyr::count(metadata_ready,Study) %>% arrange(-n) %>% dplyr::select(Study))[[1]][1]

# print it to have a reference
cat(paste("We're using as reference the values of", age.group_ref, continent_ref, study_ref, "\n"))
cat("Running Maaslin2 with", ncol(abund_data_ready), "samples","\n")

# define all variables of interest for the fixed effects
variables <- c("Disease_name","Age_group", "Continent")#, "Study") # not going to use study as a fixed effect now!
new_variables <- c("Disease_name")

# check and add variables based on conditions
for (v in variables) {
  if (v != "Disease_name" && length(unique(metadata_ready[[v]])) > 1) { # check if they have more than one value
    new_variables <- c(new_variables, v) # update the variables
  }
}

# transform those variables in to an expression
fixed <- c(new_variables, "Read_count")

# determine what variable to add as random effects
min_rep = as.numeric(sort(table(metadata_ready$Subject), decreasing=FALSE)[1])
n_studies = length(unique(metadata_ready$Study))

if (min_rep > 1 & n_studies > 1) {
  random <- c("Subject")
  fixed <- c(fixed, "Study")
} else if (min_rep == 1 & n_studies > 1) {
  random <- c("Study")
} else if (min_rep > 1 & n_studies == 1) {
  random <- c("Subject")
} else {
  random <- c("")
}

# print this information to have a reference
cat("Running Maaslin2 with fixed effects", fixed, "and random effects", random,".")

# call Maaslin2 function to fit data
fit_data <- Maaslin2(
  cores = as.numeric(cpus),
  normalization = "NONE", # already normalized into relab
  transform = "LOG", #  other options: LOG, LOGIT, AST
  min_prevalence = prevalence,
  max_significance = 0.05, # the q-value threshold for significance
  correction = "BH", # the correction method for computing the q-value 
  input_data = abund_data_ready, 
  input_metadata = metadata_ready,
  output = outfile,
  fixed_effects = c(fixed), # variables that we want to study
  reference = c("Disease_name,a_Healthy", paste0("Age_group,",age.group_ref), paste0("Continent,",continent_ref), paste0("Study,",study_ref)),
  random_effects = c(random),
  plot_heatmap = FALSE,
  plot_scatter = FALSE)
