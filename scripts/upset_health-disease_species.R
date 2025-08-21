# load libraries
library(data.table)
library(matrixStats)
library(ComplexUpset)
library(reshape2)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")

# load da data
da.results = read.csv("overlap/species_ds.csv")
da.sign = da.results[da.results$DS != 0,]

# rename diseases
colnames(da.sign) = gsub("CRC", "Colorectal cancer", colnames(da.sign))
colnames(da.sign) = gsub("CD", "Crohn's disease", colnames(da.sign))
colnames(da.sign) = gsub("UC", "Ulcerative colitis", colnames(da.sign))
colnames(da.sign) = gsub("T1D", "Type 1 diabetes", colnames(da.sign))
colnames(da.sign) = gsub("T2D", "Type 2 diabetes", colnames(da.sign))
colnames(da.sign) = gsub("MS", "Multiple sclerosis", colnames(da.sign))
colnames(da.sign) = gsub("ME", "Myalgic encephalomyelitis", colnames(da.sign))
colnames(da.sign) = gsub("T2D", "Type 2 diabetes", colnames(da.sign))
colnames(da.sign) = gsub("AS", "Ankylosing spondylitis", colnames(da.sign))
colnames(da.sign) = gsub("RA", "Rheumatoid arthritis", colnames(da.sign))
colnames(da.sign) = gsub("Parkinson", "Parkinson's disease", colnames(da.sign))
colnames(da.sign) = gsub("BoneDisease", "Bone disease", colnames(da.sign))

# process data
end.col = which(colnames(da.sign) == "DS")-1
da.health = da.sign
da.health[,2:end.col][da.health[,2:end.col] > 0] = 0
da.health[,2:end.col][da.health[,2:end.col] != 0] = 1

da.disease = da.sign
da.disease[,2:end.col][da.disease[,2:end.col] < 0] = 0
da.disease[,2:end.col][da.disease[,2:end.col] != 0] = 1

# define orders colours
orders =  c("Actinomycetales", "Bacteroidales", "Christensenellales", "Coriobacteriales", "Enterobacterales", "Erysipelotrichales", "Gastranaerophilales",
            "Lachnospirales", "Lactobacillales", "Oscillospirales", "Peptostreptococcales", "RF39", "TANB77", "Veillonellales", "Verrucomicrobiales")
colours = c(colorRampPalette(brewer.pal(12, "Set3"))(length(orders)), "darkgrey")
names(colours) = c(orders, "Other")
colours["Erysipelotrichales"] = "lightblue"
colours["TANB77"] = "steelblue1"

# function to generate upset
generate_upset = function(input_df, output_df) {
  # filter taxa
  upset.df = input_df
  upset.df = upset.df[rowSums(upset.df[,2:end.col]) > 0,]
  upset.df[,"Order"] = gsub("o__", "", upset.df[,"Order"])
  upset.df[,"Order"] = gsub("_A", "", upset.df[,"Order"])
  keep.taxa = names(sort(table(upset.df[,"Order"]), decreasing=TRUE)[1:10])
  upset.df[,"Order"] = ifelse(upset.df[,"Order"] %in% keep.taxa, upset.df[,"Order"], "Other")
  upset.df[,"Order"] = factor(upset.df[,"Order"], levels=c(sort(keep.taxa), "Other"))
  
  # upset plot
  upset.plot = upset(upset.df, colnames(upset.df)[2:end.col], n_intersections=10, name="",
                     base_annotations=list('Intersection size'=intersection_size(
                       mapping=aes(fill=Status)) 
                       + scale_fill_manual(values=c('Uncultured'='palegreen4', 'Cultured'='steelblue'), name="Status")),
                     annotations = list(
                       taxon=(
                         ggplot(mapping=aes(fill=Order))
                         + geom_bar(stat='count', position='fill')
                         + scale_y_continuous(labels=scales::percent_format())
                         + scale_fill_manual(values=colours, name="Order")
                         + ylab('% of species'))),
                     set_sizes = upset_set_size() + ylab('Number of species'),
                     width_ratio=0.25,
                     themes=upset_default_themes(text=element_text(size=14)))
  ggsave(upset.plot, filename=output_df, height=8, width=9, dpi=300)
}

# calculate proportions
unique.health = length(which(rowSums(da.health[,2:end.col]) == 1))/length(which(rowSums(da.health[,2:end.col]) > 0))*100
shared.health = length(which(rowSums(da.health[,2:end.col]) > 1))/length(which(rowSums(da.health[,2:end.col]) > 0))*100
unique.disease = length(which(rowSums(da.disease[,2:end.col]) == 1))/length(which(rowSums(da.disease[,2:end.col]) > 0))*100
shared.disease = length(which(rowSums(da.disease[,2:end.col]) > 1))/length(which(rowSums(da.disease[,2:end.col]) > 0))*100
discordant.species = da.sign[as.vector(which(rowSums(da.disease[,2:end.col]) > 0 & rowSums(da.health[,2:end.col]) > 0)),"Species_rep"]

generate_upset(da.health, "figures/upset_health_species.pdf")
generate_upset(da.disease, "figures/upset_disease_species.pdf")