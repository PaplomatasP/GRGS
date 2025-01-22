#' Run Genetic Algorithm for Feature Selection
#'
#' This function implements a genetic algorithm for feature selection in biological data,
#' integrating ligand-receptor interaction analysis and differential gene expression.
#'
#' @param data A matrix or data frame of gene expression data (samples in rows, genes in columns).
#' @param labels A factor vector of labels for each sample (e.g., "diseased" vs "healthy").
#' @param lr_data A data frame of ligand-receptor pairs with columns `ligand_gene_symbol` and `receptor_gene_symbol`.
#' This dataset is built into the package and can be accessed using `data("lr_data")`.
#' @param pop_size Integer, population size for the genetic algorithm.
#' @param generations Integer, number of generations to run.
#' @param mutation_rate Numeric, initial mutation probability (0-1).
#' @param variance_threshold Numeric, proportion of genes to retain based on variance (0-1).
#' @param early_stop_limit Integer, early stopping criteria (generations without improvement).
#' @param dynamic_mutation Logical, enable adaptive mutation rates during stagnation.
#'
#' @return A list containing:
#' \itemize{
#'   \item results - Detailed results per generation
#'   \item ga_genes - Final selected genes after pruning
#' }
#'
#' @import parallel
#' @export
#'
#' @examples
#' \donttest{
#' # Load built-in datasets
#' data("thca_gse_154763")
#' data("thca_gse_154763_Labels")
#' data("lr_data")
#'
#' thca_gse_154763 <- as.data.frame(thca_gse_154763)
#' thca_gse_154763_Labels <- as.factor(thca_gse_154763_Labels)
#' lr_data <- as.data.frame(lr_data)
#' # Run genetic algorithm
#' results <- GRGS(
#'   data = thca_gse_154763,
#'   labels = thca_gse_154763_Labels,
#'   lr_data = lr_data,
#'   pop_size = 100,
#'   generations = 500,
#'   mutation_rate = 0.03,
#'   variance_threshold = 0.8,
#'   early_stop_limit = 30,
#'   dynamic_mutation = TRUE
#' )
#'
#' # View best genes
#' print(results$ga_genes)
#'
#' # Access full results
#' str(results$results[[1]])  # First generation results
#' }
GRGS <- function(data, labels, lr_data,
                 pop_size, generations, mutation_rate,
                 variance_threshold, early_stop_limit,
                 dynamic_mutation = TRUE) {
  set.seed(1987)
  library(parallel)
  preprocessed   <- preprocess_data(data, labels, variance_threshold)
  filtered_data  <- preprocessed$filtered_data
  deg_results    <- preprocessed$deg_results

  num_features <- ncol(filtered_data) - 1
  population <- generate_population(pop_size, num_features)

  best_fitness <- 0
  early_stop_count <- 0
  results <- list()

  cl <- makeCluster(detectCores() - 2)
  clusterExport(cl, c("evaluate", "lr_data", "crossover_and_mutate", "deg_results",
                      "filtered_data"), envir = environment())

  ## Κύριο loop γενεών
  for (gen in seq_len(generations)) {  # Εδώ ξεκινάει το loop γενεών
    evaluations <- parLapply(cl, seq_len(nrow(population)), function(i) {
      res <- evaluate(
        features      = population[i, ],
        original_data = filtered_data,
        lr_data       = lr_data,
        deg_results   = deg_results,
        prune_rate    = prune_rate
      )

      if (res$Q_LR < 0.2) {
        res$fitness <- res$fitness * 0.8
      }

      return(res)
    })

    fitness_values <- sapply(evaluations, `[[`, "fitness")
    max_idx        <- which.max(fitness_values)
    best_fitness_g <- fitness_values[max_idx]

    results[[gen]] <- list(
      generation                       = gen,
      best_fitness                     = best_fitness_g,
      mean_fitness                     = mean(fitness_values, na.rm = TRUE),
      train_error                      = evaluations[[max_idx]]$train_error,
      Q_LR                             = evaluations[[max_idx]]$Q_LR,
      logFC                            = evaluations[[max_idx]]$logFC,
      best_individual                  = population[max_idx, ],
      selected_features_before_pruning = evaluations[[max_idx]]$selected_features_before_pruning,
      selected_features_after_pruning  = evaluations[[max_idx]]$selected_features_after_pruning,
      ga_genes_before_pruning          = colnames(filtered_data)[evaluations[[max_idx]]$selected_features_before_pruning],
      ga_genes_after_pruning           = colnames(filtered_data)[evaluations[[max_idx]]$selected_features_after_pruning]
    )

    cat(sprintf(
      "Generation %d: Best Fitness = %.4f, Mean Fitness = %.4f, Q_LR = %.4f, logFC = %.4f, Train Error = %.4f, Selected Features = %d\n",
      gen,
      results[[gen]]$best_fitness,
      results[[gen]]$mean_fitness,
      results[[gen]]$Q_LR,
      results[[gen]]$logFC,
      results[[gen]]$train_error,
      length(results[[gen]]$selected_features_after_pruning)
    ))


    if (best_fitness_g > best_fitness) {
      best_fitness <- best_fitness_g
      early_stop_count <- 0
      elite <- population[max_idx, ]
    } else {
      early_stop_count <- early_stop_count + 1
    }

    if (early_stop_count >= early_stop_limit) {
      cat(sprintf("Early stopping triggered at generation %d.\n", gen))
      break
    }

    if (dynamic_mutation) {
      if (early_stop_count > 25) {
        mutation_rate <- 0.95
      } else if (early_stop_count > 20) {
        mutation_rate <- min(0.9, mutation_rate * 1.5)
      } else if (early_stop_count > 15) {
        mutation_rate <- min(0.8, mutation_rate * 1.4)
      } else if (early_stop_count > 10) {
        mutation_rate <- min(0.8, mutation_rate * 1.3)
      } else if (early_stop_count > 5) {
        mutation_rate <- min(0.8, mutation_rate * 1.2)
      }
    }

    next_population <- matrix(0, nrow = pop_size, ncol = num_features)
    if (exists("elite")) {
      next_population[1, ] <- elite
    }
    for (i in 2:pop_size) {
      parents <- sample(seq_len(pop_size), 2)
      next_population[i, ] <- crossover_and_mutate(
        parent1       = population[parents[1], ],
        parent2       = population[parents[2], ],
        mutation_rate = mutation_rate,lr_data
      )
    }
    population <- next_population
  }  # Εδώ κλείνει το loop γενεών

  stopCluster(cl)

  best_gen <- which.max(sapply(results, `[[`, "best_fitness"))
  best_selected_genes <- results[[best_gen]]$selected_features_after_pruning
  best_selected_genes1 <- results[[best_gen]]$selected_features_before_pruning

  if (length(best_selected_genes) > 0 && !all(is.na(best_selected_genes))) {
    best_gene_names <- colnames(filtered_data)[best_selected_genes]
  } else {
    best_gene_names <- character(0)
  }

  return(list(
    results  = results,
    ga_genes = best_gene_names
  ))
}

