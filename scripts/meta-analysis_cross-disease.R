# libraries
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)
library(tibble)

# input files (one per disease)
setwd("~/OneDrive - University of Cambridge/MFD_shared/Projects/2023_AnaSilva_SecretBugs/data")
in_dir  = "aldex2"
pattern = "^aldex2_(.+?)\\.tsv$"
files = list.files(in_dir, pattern = pattern, full.names = TRUE)
if (!length(files)) stop("No files matched pattern in: ", in_dir)

# get sample sizes per disease to use as weights
source("../scripts/alex/metadata_disease-numbers.R")
sample_sizes = disease.counts

# custom functions
signed_z = function(effect, pval) {
  pval = pmin(pmax(pval, .Machine$double.xmin), 1 - .Machine$double.eps)
  z = qnorm(pval / 2, lower.tail = FALSE)  # two-sided
  z * sign(effect)
}

stouffer = function(z, w) {
  z = z[is.finite(z)]
  w = w[is.finite(w)]
  if (!length(z) || length(z) != length(w)) return(NA_real_)
  sum(w * z) / sqrt(sum(w^2))
}

# reader for aldex2_[disease].tsv (species col1, effect col2, pvalue col3)
read_one = function(fp) {
  nm = basename(fp)
  disease = stringr::str_match(nm, pattern)[,2]
  
  # read as tab-separated with no assumptions about headers
  df = readr::read_tsv(fp, col_names = FALSE, show_col_types = FALSE, progress = FALSE)
  
  if (ncol(df) < 3) stop("File has fewer than 3 columns: ", nm)
  
  # name the first 3 columns so we can reference them safely
  names(df)[1:3] <- c("species", "combined_effect", "combined_p")
  
  df %>%
    dplyr::transmute(
      disease = disease,
      species = as.character(.data$species),
      combined_effect = suppressWarnings(as.numeric(.data$combined_effect)),
      combined_p      = suppressWarnings(as.numeric(.data$combined_p))
    )
}

long = map_dfr(files, read_one) %>% distinct(disease, species, .keep_all = TRUE)

# check sample sizes
found_diseases = sort(unique(long$disease))
if (length(sample_sizes) == 0L) {
  stop("Please fill `sample_sizes` with named entries for these diseases: ",
       paste(found_diseases, collapse = ", "))
}
if (is.null(names(sample_sizes)) || any(!nzchar(names(sample_sizes)))) {
  stop("`sample_sizes` must be a NAMED numeric vector (names = disease identifiers).")
}
missing = setdiff(found_diseases, names(sample_sizes))
if (length(missing)) {
  stop("`sample_sizes` is missing these diseases: ",
       paste(missing, collapse = ", "),
       ". Provided names: ", paste(names(sample_sizes), collapse = ", "))
}

long = long %>%
  mutate(weight = as.numeric(sample_sizes[disease]))

# meta-analysis across diseases
meta = long %>%
  mutate(z = signed_z(combined_effect, combined_p)) %>%
  group_by(species) %>%
  summarise(
    n_diseases  = sum(is.finite(z) & is.finite(weight)),
    meta_z      = stouffer(z, w = weight),
    meta_p      = ifelse(is.finite(meta_z), 2 * pnorm(abs(meta_z), lower.tail = FALSE), NA_real_),
    
    # Weighted mean of per-disease combined_effect using the same weights
    meta_effect = {
      ok = is.finite(combined_effect) & is.finite(weight)
      if (any(ok)) sum(combined_effect[ok] * weight[ok]) / sum(weight[ok]) else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(meta_p, method = "BH")) %>%
  arrange(FDR, meta_p)

# write outputs
out_meta = file.path(in_dir, "combined_across_diseases_meta.csv")
write_csv(meta, out_meta)
write_csv(long, file.path(in_dir, "combined_across_diseases_long.csv"))
message("Wrote final output files")