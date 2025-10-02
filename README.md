# Project SecretBugs
Repository for microbiome project investigating the clinical relevance of the uncultured microbiome

# Workflows
These are the custom workflows that were used in this study. Each individual repo has detailed instructions on how to install and run the pipelines.

* [MetaGen-Fetch - Processing public metagenomes](https://github.com/alexmsalmeida/metagen-fetch)
* [metaMap - Quantifying genomes in metagenomes](https://github.com/alexmsalmeida/metamap)
* [ML-Microbiome - Machine learning classification of microbiome data](https://github.com/alexmsalmeida/ml-microbiome)

# Scripts
The `scripts/` folder contains R scripts that were used to process downstream results and plot figures.

### Metadata
* [metadata_combined.R](scripts/metadata_combined.R)
* [metadata_disease-numbers.R](scripts/metadata_disease-numbers.R)

### Diversity analysis
* [alpha-div_plot.R](scripts/alpha-div_plot.R)
* [alpha-div_rscl.R](scripts/alpha-div_rscl.R)
* [mapping_rates.R](scripts/mapping_rates.R)
* [top-prev_uncult.R](scripts/top-prev_uncult.R)

### Machine learning
* [ml_summary_best-method.R](scripts/ml_summary_best-method.R)
* [ml_summary_cross-vs-pool.R](scripts/ml_summary_cross-vs-pool.R)
* [ml_summary_pairwise.R](scripts/ml_summary_pairwise.R)
* [ml_summary_pairwise_cult-vs-all.R](scripts/ml_summary_pairwise_cult-vs-all.R)
* [ml_summary_pool-study_both-vs-cult.R](scripts/ml_summary_pool-study_both-vs-cult.R)
* [ml_summary_combined.R](scripts/ml_summary_combined.R)

### Differential abundance
* [maaslin2.R](scripts/maaslin2.R)
* [aldex2.R](scripts/aldex2.R)
* [disease_gen-overlap.R](scripts/disease_gen-overlap.R)

### Meta-analysis across diseases
* [disease_process-overlap.R](scripts/disease_process-overlap.R)
* [disease_prop_cult.R](scripts/disease_prop_cult.R)
* [disease_summarize_genera.R](scripts/disease_summarize_genera.R)
* [disease_tax-fisher.R](scripts/disease_tax-fisher.R)
* [meta-analysis_cross-disease.R](scripts/meta-analysis_cross-disease.R)
* [meta-analysis_genus-effect.R](scripts/meta-analysis_genus-effect.R)
* [upset_health-disease_species.R](scripts/upset_health-disease_species.R)
* [cag17_prev-abund.R](scripts/cag17_prev-abund.R)

### Network analysis in health
* [healthy_hub-centr.R](scripts/healthy_hub-centr.R)
* [healthy_hub-centr_strict.R](scripts/healthy_hub-centr_strict.R)
* [healthy_hub-centr_tax-status.R](scripts/healthy_hub-centr_tax-status.R)
* [healthy_mean-prev.R](scripts/healthy_mean-prev.R)
* [healthy_prev-hub_corr.R](scripts/healthy_prev-hub_corr.R)
* [healthy_prev-hub_heatmaps.R](scripts/healthy_prev-hub_heatmaps.R)

### Functional comparisons
* [ko_glmer_oscillo_genomes.R](scripts/ko_glmer_oscillo_genomes.R)
* [ko_glm_species-disease.R](scripts/ko_glm_species-disease.R)

### HMP2-IBD analyses
* [hmp2_cag170-dysb.R](scripts/hmp2_cag170-dysb.R)
* [hmp2_cag170-genes_plots.R](scripts/hmp2_cag170-genes_plots.R)
* [hmp2_cag170-genes.R](scripts/hmp2_cag170-genes.R)
