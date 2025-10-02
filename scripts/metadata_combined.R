# load libraries
library(RColorBrewer)
library(readxl)
library(data.table)
library(ggplot2)
library(reshape2)
library(grid)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggpubr)

# load input
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data/")
metadata = read.delim("metadata/metagenomes_03-2024_samples.tsv")
rownames(metadata) = metadata$Sample
metadata = metadata[!is.na(metadata$Disease.name) & !is.na(metadata$Age.group) & !is.na(metadata$Continent),]
metadata$Country = gsub("United Republic of Tanzania", "Tanzania", metadata$Country)
countries = unique(metadata$Country)

# count metagenomes per country
ddf = data.frame(matrix(NA, ncol=2, nrow=length(countries)), row.names=countries)
colnames(ddf) = c("country", "samples")
ddf$country = rownames(ddf)

for (c in ddf$country){
  samples = rownames(metadata)[which(metadata$Country == c)]
  ddf[c,2] = length(samples)
}
ddf$samples_class = NA
ddf$samples_class[which(ddf$samples >=100)] = "100-1000"
ddf$samples_class[which(ddf$samples > 1000)] = ">1000"
ddf$samples_class[which(is.na(ddf$samples_class))] = "<100"
ddf$samples_class = factor(ddf$samples_class, levels=c("<100", "100-1000", ">1000"))

# load world map using sf
world = ne_countries(scale = "medium", returnclass = "sf")
world = world %>% filter(admin != "Antarctica")

# merge sample data with world map
world$country = world$admin  # rename for clarity
world_merged = left_join(world, ddf, by = c("country"))

# plot map
map.plot = ggplot(data = world_merged, aes(fill=log10(samples+1))) +
  geom_sf(alpha = 0.8, color = "grey20", linewidth = 0.05) +
  scale_fill_gradient2(high="steelblue", name="Number of samples (log10)", na.value="white") +
  coord_sf(crs = "+proj=merc") +
  theme(panel.background = element_rect(fill = "grey95", color = "grey95", linewidth = 0.2) , panel.grid = element_blank(), panel.border = element_blank(),
        plot.margin = margin(2, 4, 2, 2),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom", legend.text = element_text(size=12), legend.title = element_text(size=14))

# save plot
source("../scripts/alex/metadata_disease-numbers.R")
source("../scripts/alex/uhgg_donut.R")
ggarrange(ggarrange(map.plot, cont.heat, common.legend=TRUE, labels=c("A", "B"), legend = "bottom", font.label = list(size=20)), 
          ggarrange(meta.plot, donut.plot, widths=c(2.5,1), labels=c("C","D"), font.label = list(size=20)), ncol=1, heights=c(1.5,1))
ggsave("figures/metadata_combined.pdf", height = 8, width = 12, dpi = 300, bg="white")
