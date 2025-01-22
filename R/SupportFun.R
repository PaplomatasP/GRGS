#' Preprocess Data for Genetic Algorithm
#'
#' This function preprocesses the input data by filtering genes based on variance and performing differential expression analysis.
#'
#' @param data A matrix or data frame of gene expression data.
#' @param labels A vector of labels for each sample.
#' @param variance_threshold A numeric value between 0 and 1 indicating the proportion of genes to retain based on variance.
#' @return A list containing the filtered data and DEG results.
#' @import edgeR
preprocess_data <- function(data, labels, variance_threshold = 0.8) {
  # Μετατροπή σε data.frame
  data <- as.data.frame(data)
  data$labels <- labels
  cat("The preprocessing of GRGS is starting. Please wait...\n")

  ## 1.1 Variance filtering
  num_genes <- ncol(data) - 1  # Αφαιρούμε τη στήλη των labels
  top_genes_count <- ceiling(num_genes * variance_threshold)
  variances <- apply(data[, -ncol(data)], 2, var)
  top_genes <- order(variances, decreasing = TRUE)[1:top_genes_count]

  filtered_data <- data[, c(top_genes, ncol(data))]

  ## 1.2 edgeR DEG analysis
  library(edgeR)
  group <- factor(filtered_data$labels)
  counts <- as.matrix(filtered_data[, -ncol(filtered_data)])

  # Μεταφορά των μετρήσεων γιατί η edgeR θέλει τις σειρές ως γονίδια
  y <- DGEList(counts = t(counts), group = group)
  y <- calcNormFactors(y)
  design <- model.matrix(~ group)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit)

  deg_results <- topTags(qlf, n = Inf)$table
  significant_genes <- which(deg_results$PValue < 0.05)
  if (length(significant_genes) < 2) {
    stop("Not enough significant genes remain after DEG analysis.")
  }

  ## 1.3 Τελικό φιλτράρισμα
  filtered_data <- filtered_data[, c(significant_genes, ncol(filtered_data))]

  return(list(
    filtered_data = filtered_data,
    deg_results   = deg_results
  ))
}

#' Evaluate Fitness of a Solution
#'
#' This function evaluates the fitness of a solution based on classification error, ligand-receptor interactions, and differential expression.
#'
#' @param features A binary vector indicating which features are active.
#' @param original_data The original data frame containing gene expression data and labels.
#' @param lr_data A data frame containing ligand-receptor pairs.
#' @param deg_results A data frame containing differential expression results.
#' @param prune_rate A numeric value indicating the proportion of features to prune.
#' @return A list containing fitness, train error, Q_LR, logFC, and selected features.
#' @import rpart
evaluate <- function(features, original_data, lr_data, deg_results, prune_rate) {
  library(rpart)

  #############################
  ## 1. Ποια features είναι ενεργά;
  #############################
  selected_features_before_pruning <- which(features == 1)

  if (length(selected_features_before_pruning) == 0) {
    # Κανένα feature ενεργό: fitness = 0
    return(list(
      fitness = 0,
      train_error = NA,
      Q_LR = NA,
      logFC = NA,
      selected_features_before_pruning = integer(0),
      selected_features_after_pruning = integer(0),
      feature_scores = numeric(0)
    ))
  }

  #############################
  ## 2. Φτιάχνουμε local_data μόνο για αυτό το evaluation
  #############################
  local_data <- original_data[, c(selected_features_before_pruning, ncol(original_data)), drop = FALSE]

  features_data <- local_data[, -ncol(local_data), drop = FALSE]
  labels <- as.factor(local_data$labels)

  #############################
  ## 3. Decision tree training
  #############################
  tree_model <- rpart(labels ~ ., data = cbind(features_data, labels), method = "class")
  predictions <- predict(tree_model, features_data, type = "class")
  train_error <- 1 - sum(diag(table(Predicted = predictions, Actual = labels))) / length(labels)

  #############################
  ## 4. Υπολογισμός Q_LR
  #############################
  selected_gene_symbols <- colnames(local_data)[seq_along(selected_features_before_pruning)]

  lr_pairs <- lr_data[
    lr_data$ligand_gene_symbol   %in% selected_gene_symbols &
      lr_data$receptor_gene_symbol %in% selected_gene_symbols,
  ]

  Q_LR <- max(nrow(lr_pairs) / nrow(lr_data), 1e-6)

  #############################
  ## 5. Υπολογισμός M_DEG (median(|logFC|)), με fallback = 0
  #############################

  # 5.1. Ποια γονίδια υπάρχουν και στο deg_results;
  common_genes <- intersect(selected_gene_symbols, rownames(deg_results))

  # 5.2. Αν δεν υπάρχει τίποτα κοινό, δίνουμε M_DEG=0
  if (length(common_genes) > 0) {
    M_DEG <- median(abs(deg_results$logFC[common_genes]), na.rm = TRUE)
    if (is.na(M_DEG)) {
      M_DEG <- 0
    }
  } else {
    M_DEG <- 0
  }

  #############################
  ## 6. Κανονικοποίηση M_DEG
  #############################
  all_logfc <- deg_results$logFC

  M_DEG <- median(abs(deg_results$logFC[selected_features_before_pruning]), na.rm = TRUE)
  M_DEG_normalized <- log1p(M_DEG) / max(log1p(deg_results$logFC), na.rm = TRUE)
  #############################
  ## 7. Fitness = w1*Q_LR + w2*M_DEG_norm + w3*(1 - train_error)
  #############################
  w1 <- 0.4
  w2 <- 0.3
  w3 <- 0.3


  raw_fitness <-     w1 * Q_LR + w2 * M_DEG_normalized + w3 * (1 - train_error)

  # Σε περίπτωση που κάποιος όρος είναι NA, μετατρέπουμε το fitness σε 0
  if (is.na(raw_fitness)) {
    fitness <- 0
  } else {
    fitness <- max(0, raw_fitness)
  }

  #############################
  ## 8. Υπολογισμός feature_scores για pruning
  #############################
  feature_scores <- numeric(length(selected_gene_symbols))
  names(feature_scores) <- selected_gene_symbols

  for (i in seq_along(selected_gene_symbols)) {
    gene <- selected_gene_symbols[i]

    # Μετατροπή σε κεφαλαία για συμβατότητα (αν χρειάζεται)
    gene_upper <- toupper(gene)

    # Q_LR local: 1 αν το γονίδιο είναι σε ligand/receptor, αλλιώς 0
    q_lr <- as.integer(
      gene_upper %in% toupper(lr_data$ligand_gene_symbol) |
        gene_upper %in% toupper(lr_data$receptor_gene_symbol)
    )

    # Υπολογισμός του normalized logFC για το συγκεκριμένο gene
    # (Βεβαιωθείτε ότι οι rownames του deg_results ταιριάζουν με τα gene names)
    if (gene %in% rownames(deg_results) && !is.na(deg_results$logFC[gene])) {
      raw_logfc_value <- log1p(abs(deg_results$logFC[gene]))
      m_deg_i_normalized <- ifelse(max_logfc > 0, raw_logfc_value / max_logfc, 0)
    } else {
      m_deg_i_normalized <- 0
    }

    # Υπολογισμός του feature score
    feature_scores[i] <- w1 * q_lr + w2 * m_deg_i_normalized
  }


  #############################
  ## 9. Pruning
  #############################
  prune_rate <- 0.05

  lr_genes <- which(feature_scores > w1)
  non_lr_genes <- which(feature_scores < w1)

  if (length(lr_genes) > 0) {
    selected_features_after_pruning <- lr_genes
  } else {
    # Ορίζουμε ποια είναι τα top indices με βάση το feature_scores
    top_indices <- order(feature_scores[non_lr_genes], decreasing = TRUE)

    num_to_prune <- ceiling(length(top_indices) * prune_rate)
    # Αν το num_to_prune > length(non_lr_genes), το pmin το διορθώνει
    num_to_prune <- pmin(num_to_prune, length(top_indices))

    # Τώρα παίρνουμε τους πρώτους num_to_prune
    selected_features_after_pruning <- non_lr_genes[ top_indices[seq_len(num_to_prune)] ]
  }



  #############################
  ## 10. Επιστροφή
  #############################
  return(list(
    fitness                           = fitness,
    train_error                       = train_error,
    Q_LR                              = Q_LR,
    logFC                             = M_DEG_normalized,
    selected_features_before_pruning  = selected_features_before_pruning,
    selected_features_after_pruning   = selected_features_after_pruning,
    feature_scores                    = feature_scores
  ))
}



#' Generate Initial Population
#'
#' This function generates an initial population of solutions for the genetic algorithm.
#'
#' @param pop_size The size of the population.
#' @param num_features The number of features in the data.
#' @param prob The probability of a feature being active in the initial population.
#' @return A matrix representing the initial population.

generate_population <- function(pop_size, num_features, prob = 0.1) {
  if (pop_size <= 0 || num_features <= 0) {
    stop("pop_size and num_features must be positive integers.")
  }
  # Δημιουργία matrix pop_size x num_features, με πιθανότητα prob να είναι 1
  population <- matrix(
    rbinom(num_features * pop_size, size = 1, prob = prob),
    nrow = pop_size,
    ncol = num_features
  )
  return(population)
}


#' Perform Crossover and Mutation
#'
#' This function performs crossover and mutation to generate offspring.
#'
#' @param parent1 A binary vector representing the first parent.
#' @param parent2 A binary vector representing the second parent.
#' @param lr_data A data frame containing ligand-receptor pairs.
#' @param mutation_rate The probability of mutation.
#' @param lr_mutation_boost The increased mutation probability for ligand-receptor genes.
#' @return A binary vector representing the offspring.
crossover_and_mutate <- function(parent1, parent2, lr_data, mutation_rate = 0.1, lr_mutation_boost = 0.3) {
  if (length(parent1) != length(parent2)) {
    stop("Both parents must have the same length.")
  }
  if (mutation_rate < 0 || mutation_rate > 1) {
    stop("Mutation rate must be between 0 and 1.")
  }

  # Σημεία crossover
  crossover_points <- sort(sample(1:(length(parent1) - 1), 2))

  # Δημιουργία offspring
  offspring <- numeric(length(parent1))
  offspring[1:crossover_points[1]] <- parent1[1:crossover_points[1]]
  offspring[(crossover_points[1] + 1):crossover_points[2]] <- parent2[(crossover_points[1] + 1):crossover_points[2]]
  offspring[(crossover_points[2] + 1):length(parent1)] <- parent1[(crossover_points[2] + 1):length(parent1)]

  # Εύρεση LR γονιδίων στις θέσεις του offspring
  lr_genes <- which(names(parent1) %in% c(lr_data$ligand_gene_symbol, lr_data$receptor_gene_symbol))

  # Διατήρηση LR γονιδίων με μεγαλύτερη πιθανότητα να είναι ενεργά (1)
  offspring[lr_genes] <- sample(c(0, 1), length(lr_genes),
                                prob = c(0.2, 0.8), replace = TRUE)

  # Κλασσικό mutation με αυξημένη πιθανότητα για LR γονίδια
  mutation_indices <- which(runif(length(offspring)) < mutation_rate)

  # Ειδική μεταχείριση για τα LR γονίδια - αύξηση πιθανότητας μετάλλαξης
  lr_mutation_indices <- intersect(mutation_indices, lr_genes)
  offspring[lr_mutation_indices] <- sample(c(0, 1), length(lr_mutation_indices),
                                           prob = c(1 - lr_mutation_boost, lr_mutation_boost), replace = TRUE)

  # Κανονική μετάλλαξη για τα υπόλοιπα γονίδια
  non_lr_mutation_indices <- setdiff(mutation_indices, lr_genes)
  offspring[non_lr_mutation_indices] <- 1 - offspring[non_lr_mutation_indices]

  return(offspring)
}


