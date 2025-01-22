---
title: "Optimizing Gene Selection in High-Dimensional Data Using the GRGS Framework: A Genetic Algorithm Approach"
output: github_document
---

## Authors
Name1, Name2, Name3  
Department of Bioinformatics, University/Institute  
Corresponding Author Email: example@university.edu  

## Abstract
The identification of significant genes is critical for understanding disease mechanisms and advancing biomarker discovery for diagnostics and therapeutics. In this study, we present the Genetic Receptor-Gene Selector (GRGS), a novel framework that leverages Genetic Algorithms (GAs) to overcome the challenges posed by high-dimensional gene expression data. GRGS integrates variance filtering, differential expression analysis, and ligand-receptor (LR) enrichment to enhance both predictive accuracy and biological relevance. The framework employs edgeR normalization for precise differential expression identification, variance filtering to mitigate noise, and a multi-objective fitness function that optimizes decision tree performance, ligand-receptor pair enrichment, and the median log2 fold change of selected genes. This comprehensive approach ensures robust feature selection while maintaining computational efficiency. The performance of GRGS was evaluated using single-cell RNA sequencing data from Familial Medullary Thyroid Carcinoma, where it identified a biologically meaningful subset of 124 genes. GRGS outperformed the full dataset and three widely used feature selection methods—Lasso, Boruta, and standard Genetic Algorithm (GA)—in predictive accuracy and biological interpretability across multiple machine learning classifiers. These results highlight the potential of GRGS as an advanced tool for biomarker discovery and precision medicine applications.

## Keywords
Genetic Algorithm, Feature Selection, Gene Expression, Ligand-Receptor Interaction, Differential Expression, Biomarker Discovery, Precision Medicine

---

## GRGS {GRGS} R Documentation

### Run Genetic Algorithm for Feature Selection

#### Description
This function implements a genetic algorithm for feature selection in biological data, integrating ligand-receptor interaction analysis and differential gene expression.

#### Usage
```r
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

Βεβαίως! Ακολουθεί το περιεχόμενο σε R Markdown μορφή για να μπορείς να το κάνεις απλά αντιγραφή και επικόλληση:

```rmd
---
title: "Optimizing Gene Selection in High-Dimensional Data Using the GRGS Framework: A Genetic Algorithm Approach"
output: github_document
---

## Authors
Name1, Name2, Name3  
Department of Bioinformatics, University/Institute  
Corresponding Author Email: example@university.edu  

## Abstract
The identification of significant genes is critical for understanding disease mechanisms and advancing biomarker discovery for diagnostics and therapeutics. In this study, we present the Genetic Receptor-Gene Selector (GRGS), a novel framework that leverages Genetic Algorithms (GAs) to overcome the challenges posed by high-dimensional gene expression data. GRGS integrates variance filtering, differential expression analysis, and ligand-receptor (LR) enrichment to enhance both predictive accuracy and biological relevance. The framework employs edgeR normalization for precise differential expression identification, variance filtering to mitigate noise, and a multi-objective fitness function that optimizes decision tree performance, ligand-receptor pair enrichment, and the median log2 fold change of selected genes. This comprehensive approach ensures robust feature selection while maintaining computational efficiency. The performance of GRGS was evaluated using single-cell RNA sequencing data from Familial Medullary Thyroid Carcinoma, where it identified a biologically meaningful subset of 124 genes. GRGS outperformed the full dataset and three widely used feature selection methods—Lasso, Boruta, and standard Genetic Algorithm (GA)—in predictive accuracy and biological interpretability across multiple machine learning classifiers. These results highlight the potential of GRGS as an advanced tool for biomarker discovery and precision medicine applications.

## Keywords
Genetic Algorithm, Feature Selection, Gene Expression, Ligand-Receptor Interaction, Differential Expression, Biomarker Discovery, Precision Medicine

---

## GRGS {GRGS} R Documentation

### Run Genetic Algorithm for Feature Selection

#### Description
This function implements a genetic algorithm for feature selection in biological data, integrating ligand-receptor interaction analysis and differential gene expression.

#### Usage
```r
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
```

#### Arguments
- **data**: A matrix or data frame of gene expression data (samples in rows, genes in columns).
- **labels**: A factor vector of labels for each sample (e.g., "diseased" vs "healthy").
- **lr_data**: A data frame of ligand-receptor pairs with columns 'ligand_gene_symbol' and 'receptor_gene_symbol'. This dataset is built into the package and can be accessed using `data("lr_data")`.
- **pop_size**: Integer, population size for the genetic algorithm.
- **generations**: Integer, number of generations to run.
- **mutation_rate**: Numeric, initial mutation probability (0-1).
- **variance_threshold**: Numeric, proportion of genes to retain based on variance (0-1).
- **early_stop_limit**: Integer, early stopping criteria (generations without improvement).
- **dynamic_mutation**: Logical, enable adaptive mutation rates during stagnation.

#### Value
A list containing:
- **results**: Detailed results per generation
- **ga_genes**: Final selected genes after pruning

#### Examples
```r
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
```

```r
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

```
