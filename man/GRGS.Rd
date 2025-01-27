% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRGS.R
\name{GRGS}
\alias{GRGS}
\title{Run Genetic Algorithm for Feature Selection}
\usage{
GRGS(
  data,
  labels,
  lr_data,
  pop_size,
  generations,
  mutation_rate,
  variance_threshold,
  early_stop_limit,
  dynamic_mutation = TRUE
)
}
\arguments{
\item{data}{A matrix or data frame of gene expression data (samples in rows, genes in columns).}

\item{labels}{A factor vector of labels for each sample (e.g., "diseased" vs "healthy").}

\item{lr_data}{A data frame of ligand-receptor pairs with columns `ligand_gene_symbol` and `receptor_gene_symbol`.
This dataset is built into the package and can be accessed using `data("lr_data")`.}

\item{pop_size}{Integer, population size for the genetic algorithm.}

\item{generations}{Integer, number of generations to run.}

\item{mutation_rate}{Numeric, initial mutation probability (0-1).}

\item{variance_threshold}{Numeric, proportion of genes to retain based on variance (0-1).}

\item{early_stop_limit}{Integer, early stopping criteria (generations without improvement).}

\item{dynamic_mutation}{Logical, enable adaptive mutation rates during stagnation.}
}
\value{
A list containing:
\itemize{
  \item results - Detailed results per generation
  \item ga_genes - Final selected genes after pruning
}
}
\description{
This function implements a genetic algorithm for feature selection in biological data,
integrating ligand-receptor interaction analysis and differential gene expression.
}
\examples{
\donttest{
# Load built-in datasets
data("thca_gse_154763")
data("thca_gse_154763_Labels")
data("lr_data")

thca_gse_154763 <- as.data.frame(thca_gse_154763)
thca_gse_154763_Labels <- as.factor(thca_gse_154763_Labels)
lr_data <- as.data.frame(lr_data)
# Run genetic algorithm
results <- GRGS(
  data = thca_gse_154763,
  labels = thca_gse_154763_Labels,
  lr_data = lr_data,
  pop_size = 100,
  generations = 500,
  mutation_rate = 0.03,
  variance_threshold = 0.8,
  early_stop_limit = 30,
  dynamic_mutation = TRUE
)

# View best genes
print(results$ga_genes)

# Access full results
str(results$results[[1]])  # First generation results
}
}
