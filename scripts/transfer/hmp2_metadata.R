# load data
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/metadata")
dysb.scores = read.delim("hmp2_ibd/dysbiosis_scores.tsv", header=FALSE)
colnames(dysb.scores) = c("External.ID", "Dysb_score", "Dysb")
metadata.run = read.delim("hmp2_ibd/SRP115494_run-metadata.tsv")
metadata.run$sample_alias = gsub("_.*", "", metadata.run$sample_alias)

metadata.all = read.delim("metagenomes_09-2023.tsv")
metadata.subs = metadata.all[which(metadata.all$Study == "SRP115494"),]
metadata.subs$Sample_alias = metadata.run[match(metadata.subs$Sample, metadata.run$secondary_sample_accession),"sample_alias"]

metadata.hmp2 = read.csv("hmp2_ibd/hmp2_metadata_2018-08-20.csv")
metadata.hmp2 = unique(metadata.hmp2[which(metadata.hmp2$data_type == "metagenomics"),c("External.ID", "site_sub_coll", "week_num", "visit_num")])
metadata.hmp2 = metadata.hmp2[which(metadata.hmp2$External.ID %in% dysb.scores$External.ID),]
metadata.hmp2 = metadata.hmp2[!grepl("_", metadata.hmp2$External.ID),]
colnames(metadata.hmp2)[2] = "Sample_alias"

metadata.merged = merge(metadata.subs, metadata.hmp2, by="Sample_alias")
metadata.merged = merge(metadata.merged, dysb.scores, by="External.ID")
metadata.merged = metadata.merged[,c("Run", "Sample", "Subject", "week_num", "visit_num", "Health.state", "Disease.name", "Dysb_score", "Age", "Age.group", "Gender")]

# save metadata
write.table(metadata.merged, file="hmp2_ibd/metagenomes_dysb-score.tsv", sep="\t", row.names=FALSE, quote=FALSE)
