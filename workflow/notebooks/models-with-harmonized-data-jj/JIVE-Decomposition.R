library(r.jive)

# Paths
assays <- list(
  plasma_cytokine_concentrations='results/processed/harmonized/training_plasma_cytokine_concentrations.tsv',
  pbmc_cell_frequency='results/processed/harmonized/training_pbmc_cell_frequency.tsv',
  plasma_antibody_levels='results/processed/harmonized/training_plasma_antibody_levels.tsv',
  pbmc_gene_expression='results/processed/harmonized/training_pbmc_gene_expression.tsv'
)

common_samples <- intersect_all(lapply(assays, function(fn) {
  read.table(fn, sep='\t', header=TRUE)$subject_id
}))

jive_data <- lapply(assays, function(fn) {
  data <- read.table(fn, sep='\t', header=TRUE)
  rownames(data) <- data$subject_id
  data[common_samples, -1]
})

factorizations_jive <- jive(jive_data, rankJ=10, rankA=rep(10, length(jive_data)), method="given")
saveRDS(factorizations_jive, "results/jive/harmonized/jive_decomposition.rds")
