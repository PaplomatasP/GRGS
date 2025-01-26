feature_selection <- function(full_data, label_col, method = c("RFE", "edgeR", "GA", "Boruta", "LASSO")) {
  library(caret)
  library(edgeR)
  library(GA)
  library(Boruta)
  library(class)
  library(glmnet)  # Προσθήκη βιβλιοθήκης για LASSO
  
  # Εξασφάλιση σωστής επιλογής μεθόδου
  method <- match.arg(method)
  
  ## 1.1 Variance filtering
  num_genes <- ncol(full_data) - 1  # Αφαιρούμε τη στήλη των labels
  top_genes_count <- ceiling(num_genes * 0.8)
  variances <- apply(full_data[, -ncol(full_data)], 2, var)
  top_genes <- order(variances, decreasing = TRUE)[1:top_genes_count]
  
  full_data <- full_data[, c(top_genes, ncol(full_data))]
  
  # RFE Method
  run_rfe <- function() {
    cat("Τώρα εκτελείται η μέθοδος RFE...\n")
    control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
    rfe_results <- rfe(
      x = full_data[, -which(colnames(full_data) == label_col)],
      y = as.factor(full_data[[label_col]]),
      sizes = c(10, 20, 50, 100),
      rfeControl = control
    )
    selected_genes <- predictors(rfe_results)
    return(selected_genes)
  }
  
  # edgeR Method
  run_edgeR <- function() {
    group <- factor(full_data[[label_col]])
    counts <- as.matrix(full_data[, -which(colnames(full_data) == label_col)])
    y <- DGEList(counts = t(counts), group = group)
    y <- calcNormFactors(y)
    design <- model.matrix(~ group)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit)
    deg_results <- topTags(qlf, n = Inf)$table
    significant_genes <- rownames(deg_results[deg_results$FDR < 0.05 & abs(deg_results$logFC) > 1, ])
    return(significant_genes)
  }
  
  # GA Method (with kNN)
  run_ga <- function() {
    cat("Τώρα εκτελείται η μέθοδος GA...\n")
    set.seed(123)
    fitness_function <- function(selected_genes) {
      selected_genes <- as.logical(selected_genes)
      if (sum(selected_genes) == 0) return(0)
      trainIndex <- createDataPartition(full_data[[label_col]], p = 0.7, list = FALSE)
      train_data <- full_data[trainIndex, ]
      test_data <- full_data[-trainIndex, ]
      train_subset <- train_data[, c(which(selected_genes), which(colnames(full_data) == label_col))]
      test_subset <- test_data[, c(which(selected_genes), which(colnames(full_data) == label_col))]
      train_x <- train_subset[, -ncol(train_subset)]
      train_y <- train_subset[[label_col]]
      test_x <- test_subset[, -ncol(test_subset)]
      test_y <- test_subset[[label_col]]
      preds <- knn(train = train_x, test = test_x, cl = train_y, k = 3)
      accuracy <- sum(preds == test_y) / length(test_y)
      penalty_weight <- 0.0005
      max_genes_allowed <- 600
      penalty <- sum(selected_genes) / ncol(full_data) * penalty_weight
      if (sum(selected_genes) > max_genes_allowed) {
        penalty <- penalty + (sum(selected_genes) - max_genes_allowed) * 0.0001
      }
      reward_for_fewer_genes <- (1 / (sum(selected_genes) + 1)) * 0.8
      return(accuracy - penalty + reward_for_fewer_genes)
    }
    ga_result <- ga(type = "binary", fitness = fitness_function, nBits = ncol(full_data) - 1, 
                                     maxiter = 10, popSize = 10, run = 10)
    selected_genes <- colnames(full_data)[which(ga_result@solution[1, ] == 1)]
    return(selected_genes)
  }
  
  # Boruta Method
  run_boruta <- function() {
    cat("Τώρα εκτελείται η μέθοδος Boruta...\n")
    set.seed(123)
    formula <- as.formula(paste(label_col, "~ ."))
    boruta_result <- Boruta(formula, data = full_data, doTrace = 2)
    selected_genes <- getSelectedAttributes(boruta_result, withTentative = TRUE)
    return(selected_genes)
  }
  
  # LASSO Method
  run_lasso <- function() {
    cat("Τώρα εκτελείται η μέθοδος LASSO...\n")
    x <- as.matrix(full_data[, -which(colnames(full_data) == label_col)])
    y <- as.factor(full_data[[label_col]])
    lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)
    selected <- as.matrix(coef(lasso_model, s = "lambda.min"))
    nonzero_indices <- which(selected != 0)
    selected_genes <- rownames(selected)[nonzero_indices]
    selected_genes <- selected_genes[selected_genes != "(Intercept)"]
    return(selected_genes)
  }
  
  # Επιλογή μεθόδου
  selected_genes <- switch(
    method,
    "RFE" = run_rfe(),
    "edgeR" = run_edgeR(),
    "GA" = run_ga(),
    "Boruta" = run_boruta(),
    "LASSO" = run_lasso()
  )
  
  return(list(method = method, genes = selected_genes))
}




full_data = filtered_data
full_data$labels = as.factor(thca_gse_154763_Labels)


# Δημιουργία λίστας για αποθήκευση αποτελεσμάτων
methods_to_run <- c( "GA", "Boruta","RFE")
results <- list()

# Τρέχουμε τις επιλεγμένες μεθόδους
for (method in methods_to_run) {
  cat("Running method:", method, "\n")
  results[[method]] <- feature_selection(full_data = full_data, label_col = "labels", method = method)
}

# Εμφάνιση αποτελεσμάτων
for (method in names(results)) {
  cat("\nMethod:", method, "\n")
  print(results[[method]])
}

result_edgeR <- feature_selection(full_data = full_data, label_col = "labels", method = "edgeR")

result_LASSO <- feature_selection(full_data, label_col = "labels", method = "LASSO")

 
# Χρήση της συνάρτησης
result_rfe <- feature_selection(full_data = full_data, label_col = "labels", method = "RFE")
#result_edger <- feature_selection(full_data = full_data, label_col = "labels", method = "edgeR")
result_ga <- feature_selection(full_data = filtered_data, label_col = "labels", method = "GA")
Boruta_result <- feature_selection(full_data, label_col = "labels", method = "Boruta")

reduce_genes_randomly <- function(genes, top_n = 500) {
  sampled_genes <- sample(genes, top_n)
  cat("Αρχικός αριθμός γονιδίων:", length(genes), "\n")
  cat("Γονίδια μετά την τυχαία επιλογή:", length(sampled_genes), "\n")
  return(sampled_genes)
}

# Εφαρμογή:
reduced_genes <- reduce_genes_randomly(
  genes = results$GA$genes,
  top_n = 836
)
reduced_genes
reduced_genes
results$GA$genes= reduced_genes
