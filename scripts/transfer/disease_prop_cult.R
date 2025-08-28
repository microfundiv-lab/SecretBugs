# load libraries
library(data.table)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# set working directory
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/")

# load da data
da.results = read.csv("overlap/species_ds.csv")
da.sign = da.results[da.results$DS != 0,]
da.sign$Classification = ifelse(da.sign$DS > 0, "Disease", "Health")
da.status = as.data.frame(as.matrix(table(da.sign[,c("Classification", "Status")])))

# plot barplot
bar.plot = ggplot(da.status, aes(x=Classification, y=Freq, fill=Status)) +
  geom_bar(stat='identity', colour="grey50", linewidth=0.2) +
  theme_classic() +
  ylab("Number of significant species") +
  scale_fill_manual(values=c("#b4c8ff", "#d9ffc1")) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(size=12)) +
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size=12))
ggsave(filename = "figures/prop_cult-status_disease.pdf", dpi=300, height=5, width=4)

# check significance
chisq.res = chisq.test(as.matrix(table(da.sign[,c("Classification", "Status")])))
print(chisq.res)
