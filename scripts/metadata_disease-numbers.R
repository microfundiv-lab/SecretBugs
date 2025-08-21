#load libraries
library(ggplot2)
library(pheatmap)

#load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/")
metadata = read.delim("metadata/metagenomes_03-2024_samples.tsv")
rownames(metadata) = metadata$Sample
metadata = metadata[!is.na(metadata$Disease.name) & !is.na(metadata$Age.group) & !is.na(metadata$Continent),]

#metadata by disease
metadata$Disease.clean = metadata$Disease.name
metadata$Disease.clean = gsub("CRC", "Colorectal cancer", metadata$Disease.clean)
metadata$Disease.clean = gsub("CD", "Crohn's disease", metadata$Disease.clean)
metadata$Disease.clean = gsub("UC", "Ulcerative colitis", metadata$Disease.clean)
metadata$Disease.clean = gsub("T1D", "Type 1 diabetes", metadata$Disease.clean)
metadata$Disease.clean = gsub("T2D", "Type 2 diabetes", metadata$Disease.clean)
metadata$Disease.clean = gsub("MS", "Multiple sclerosis", metadata$Disease.clean)
metadata$Disease.clean = gsub("ME", "Myalgic encephalomyelitis", metadata$Disease.clean)
metadata$Disease.clean = gsub("T2D", "Type 2 diabetes", metadata$Disease.clean)
metadata$Disease.clean = gsub("AS", "Ankylosing spondylitis", metadata$Disease.clean)
metadata$Disease.clean = gsub("RA", "Rheumatoid arthritis", metadata$Disease.clean)
metadata$Disease.clean = gsub("Parkinson", "Parkinson's disease", metadata$Disease.clean)
metadata$Disease.clean = gsub("BoneDisease", "Bone disease", metadata$Disease.clean)
selected.diseases = c("Adenoma", "Ankylosing spondylitis", "Bone disease", "Crohn's disease", "Colorectal cancer", "Myalgic encephalomyelitis", "Multiple sclerosis", 
                      "Obesity", "Parkinson's disease", "Rheumatoid arthritis", "Type 1 diabetes", "Type 2 diabetes", "Ulcerative colitis")

# save supplementary table
write.table(metadata, file="tables/supptable1.tsv", row.names=FALSE, quote=FALSE, sep="\t")

#prepare df by disease and continent
disease.df = data.frame(matrix(nrow = length(selected.diseases), ncol=3))
rownames(disease.df) = selected.diseases
colnames(disease.df) = c("Healthy", "Diseased", "Studies")
continents = c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")
disease.cont = data.frame(matrix(nrow = length(selected.diseases), ncol=length(continents)))
rownames(disease.cont) = selected.diseases
colnames(disease.cont) = continents

studies.all = c()
for (disease in selected.diseases) {
  disease.studies = unique(metadata[which(metadata$Disease.clean == disease),"Study"])
  healthy.studies = unique(metadata[which(metadata$Disease.clean == "Healthy"),"Study"])
  studies = intersect(disease.studies, healthy.studies)
  studies.all = unique(c(studies.all, studies))
  disease.df[disease,"Healthy"] = nrow(metadata[which(metadata$Disease.clean == "Healthy" & metadata$Study %in% studies),])
  disease.df[disease,"Diseased"] = nrow(metadata[which(metadata$Disease.clean == disease & metadata$Study %in% studies),])
  disease.df[disease,"Studies"] = length(studies)
  for (c in continents) {
    disease.cont[disease,c] = nrow(metadata[which(metadata$Continent == c & metadata$Study %in% studies),])
  }
}
disease.df$Name = metadata[match(rownames(disease.df), metadata$Disease.clean),"Disease.name"]
disease.counts = rowSums(disease.df[,c("Healthy", "Diseased")])
names(disease.counts) = disease.df$Name

# plot heatmap
disease.cont.melt = reshape2::melt(as.matrix(log10(disease.cont+1)))
cont.heat = ggplot(disease.cont.melt, aes(x=Var2, y=Var1, fill=value)) +
  geom_tile() +
  theme_classic() +
  scale_fill_gradient2(high="steelblue", name="Number of samples (log10)") +
  scale_y_discrete(limits = rev(levels(disease.cont.melt$Var1))) +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(size=12, angle=45, vjust=1, hjust=1)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title = element_blank())

# plot barplot
disease.order = names(sort(rowSums(disease.df[,c("Healthy", "Diseased")]), decreasing=TRUE))
disease.df.melt = reshape2::melt(as.matrix(disease.df[,c("Healthy", "Diseased")]))
meta.plot = ggplot(disease.df.melt, aes(x=Var1, y=log10(value), fill=Var2)) +
  geom_bar(stat="identity", alpha=0.85, width=0.5, position="dodge") +
  theme_classic() +
  scale_x_discrete(limits=disease.order) +
  scale_fill_manual(values=c("steelblue", "tomato"), name="Sample type") +
  ylab("log10(Number of samples)") +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1)) +
  theme(axis.text.y = element_text(size=12))
ggsave(file="figures/metadata_disease-numbers.png", height=4, width=10, dpi=100)

# count by disease group
studies.all = unique(metadata[which(metadata$Disease.clean %in% selected.diseases),"Study"])
metadata.disease = metadata[which(metadata$Study %in% studies.all),]
gastro = c("Adenoma", "Crohn's disease", "Ulcerative colitis", "Colorectal cancer")
gastro.prop = sum(disease.df[gastro,c(2)])

metabo = c("Obesity", "Type 1 diabetes", "Type 2 diabetes")
metabo.prop = sum(disease.df[metabo,c(2)])

neuro = c("Myalgic encephalomyelitis", "Multiple sclerosis", "Parkinson's disease")
neuro.prop = sum(disease.df[neuro,c(2)])

other = c("Ankylosing spondylitis", "Bone disease", "Rheumatoid arthritis")
other.prop = sum(disease.df[other,c(2)])

# calculate healthy samples
disease.studies = unique(metadata[which(metadata$Disease.clean != "Healthy"),"Study"])
metadata.healthy = metadata[which(!metadata$Study %in% disease.studies),]
