classification_comparison_CV <- function(full_data, 
                                         selected_genes, 
                                         label_col,
                                         models = c("GBM", "RandomForest", "kNN", "NaiveBayes"), 
                                         folds = 10) {
  # Φόρτωση πακέτων
  library(caret)
  library(randomForest)
  library(pROC)
  library(gbm)
  library(FNN)
  library(doParallel)
  library(e1071)  # Για τη συνάρτηση naiveBayes
  
  # Δημιουργία parallel cluster
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  
  set.seed(2111987)
  
  # Δημιουργία folds για Cross-Validation
  fold_indices <- createFolds(full_data[[label_col]], k = folds, list = TRUE)
  
  # Δομή για αποθήκευση αποτελεσμάτων
  results <- list()
  
  # Επανάληψη για κάθε μοντέλο
  for (model in models) {
    cat(sprintf("\nEvaluating %s model with Cross-Validation...\n", model))
    
    # Λίστες για συλλογή των metrics
    full_metrics <- list(AUC = c(), Accuracy = c(), Sensitivity = c(), Specificity = c(), F1 = c())
    subset_metrics <- list(AUC = c(), Accuracy = c(), Sensitivity = c(), Specificity = c(), F1 = c())
    
    # Λίστες για συλλογή των προβλέψεων ανά fold
    predictions_full <- list()
    predictions_subset <- list()
    
    # Cross-Validation loop
    for (i in seq_along(fold_indices)) {
      cat(sprintf("Processing fold %d/%d...\n", i, folds))
      
      # Ορισμός training και test set
      test_indices <- fold_indices[[i]]
      train_data <- full_data[-test_indices, ]
      test_data  <- full_data[test_indices, ]
      
      # Επιλογή υποσυνόλου γονιδίων (εφόσον υπάρχουν στο dataset)
      subset_genes <- intersect(selected_genes, colnames(train_data))
      subset_train <- train_data[, c(subset_genes, label_col), drop = FALSE]
      subset_test  <- test_data[,  c(subset_genes, label_col), drop = FALSE]
      
      # Κανονικοποίηση (Scaling & Centering) στο πλήρες dataset
      preProc <- preProcess(
        train_data[, -which(colnames(train_data) == label_col)], 
        method = c("center", "scale")
      )
      train_data[, -which(colnames(train_data) == label_col)] <- predict(
        preProc, 
        train_data[, -which(colnames(train_data) == label_col)]
      )
      test_data[, -which(colnames(test_data) == label_col)] <- predict(
        preProc, 
        test_data[, -which(colnames(test_data) == label_col)]
      )
      
      # Κανονικοποίηση (Scaling & Centering) στο subset
      preProc_subset <- preProcess(
        subset_train[, -which(colnames(subset_train) == label_col)], 
        method = c("center", "scale")
      )
      subset_train[, -which(colnames(subset_train) == label_col)] <- predict(
        preProc_subset, 
        subset_train[, -which(colnames(subset_train) == label_col)]
      )
      subset_test[, -which(colnames(subset_test) == label_col)] <- predict(
        preProc_subset, 
        subset_test[, -which(colnames(subset_test) == label_col)]
      )
      
      # Εκπαίδευση μοντέλου & Προβλέψεις
      if (model == "GBM") {
        
        gbm_full <- gbm(
          as.formula(paste(label_col, "~ .")), 
          data = train_data, 
          distribution = "bernoulli", 
          n.trees = 100
        )
        gbm_subset <- gbm(
          as.formula(paste(label_col, "~ .")), 
          data = subset_train, 
          distribution = "bernoulli", 
          n.trees = 100
        )
        
        pred_full_probs <- predict(gbm_full, newdata = test_data,    type = "response", n.trees = 100)
        pred_subset_probs <- predict(gbm_subset, newdata = subset_test, type = "response", n.trees = 100)
        
        # Μετατροπή πιθανοτήτων σε predicted labels
        # Εφόσον τα levels του label είναι π.χ. c("healthy","diseases"), το δεύτερο level[2] είναι "diseases"
        pred_full_labels <- ifelse(pred_full_probs > 0.5, 
                                   levels(test_data[[label_col]])[2], 
                                   levels(test_data[[label_col]])[1])
        pred_subset_labels <- ifelse(pred_subset_probs > 0.5, 
                                     levels(subset_test[[label_col]])[2], 
                                     levels(subset_test[[label_col]])[1])
        
      } else if (model == "RandomForest") {
        
        rf_full <- randomForest(as.formula(paste(label_col, "~ .")), data = train_data)
        rf_subset <- randomForest(as.formula(paste(label_col, "~ .")), data = subset_train)
        
        pred_full_probs   <- predict(rf_full,   newdata = test_data,    type = "prob")[, 2]
        pred_subset_probs <- predict(rf_subset, newdata = subset_test,  type = "prob")[, 2]
        
        # Για να πάρουμε τα predicted labels απευθείας από το RF
        pred_full_labels   <- predict(rf_full,   newdata = test_data)
        pred_subset_labels <- predict(rf_subset, newdata = subset_test)
        
      } else if (model == "kNN") {
        
        pred_full_labels <- knn(
          train = train_data[, -which(colnames(train_data) == label_col)],
          test  = test_data[, -which(colnames(test_data) == label_col)],
          cl    = train_data[[label_col]], 
          k     = 5
        )
        
        pred_subset_labels <- knn(
          train = subset_train[, -which(colnames(subset_train) == label_col)],
          test  = subset_test[, -which(colnames(subset_test) == label_col)],
          cl    = subset_train[[label_col]], 
          k     = 5
        )
        
        # Στο kNN δεν έχουμε άμεση πιθανότητα, οπότε χονδρικά τη φτιάχνουμε ως 0/1
        # (π.χ. 0 για healthy, 1 για diseases). Αν το επίπεδο 2 είναι "diseases", θα το ορίσουμε ως 1
        pred_full_probs   <- ifelse(pred_full_labels   == levels(test_data[[label_col]])[2], 1, 0)
        pred_subset_probs <- ifelse(pred_subset_labels == levels(subset_test[[label_col]])[2], 1, 0)
        
      } else if (model == "NaiveBayes") {
        
        nb_full   <- naiveBayes(as.formula(paste(label_col, "~ .")), data = train_data)
        nb_subset <- naiveBayes(as.formula(paste(label_col, "~ .")), data = subset_train)
        
        pred_full_probs   <- predict(nb_full,   newdata = test_data,    type = "raw")[, 2]
        pred_subset_probs <- predict(nb_subset, newdata = subset_test,  type = "raw")[, 2]
        
        pred_full_labels <- ifelse(pred_full_probs > 0.5, 
                                   levels(test_data[[label_col]])[2], 
                                   levels(test_data[[label_col]])[1])
        pred_subset_labels <- ifelse(pred_subset_probs > 0.5, 
                                     levels(subset_test[[label_col]])[2], 
                                     levels(subset_test[[label_col]])[1])
      }
      
      # Αποθήκευση προβλέψεων (πιθανοτήτων) για μελλοντική χρήση
      predictions_full[[i]] <- list(labels = test_data[[label_col]], predictions = pred_full_probs)
      predictions_subset[[i]] <- list(labels = subset_test[[label_col]], predictions = pred_subset_probs)
      
      # Μετατροπή των πραγματικών ετικετών σε 0/1 (0=healthy, 1=diseases)
      test_labels_numeric    <- ifelse(test_data[[label_col]]    == levels(test_data[[label_col]])[2], 1, 0)
      subset_labels_numeric  <- ifelse(subset_test[[label_col]]  == levels(subset_test[[label_col]])[2], 1, 0)
      
      # Υπολογισμός AUC
      auc_full   <- roc(test_labels_numeric,    pred_full_probs)$auc
      auc_subset <- roc(subset_labels_numeric,  pred_subset_probs)$auc
      
      # Μετατρέπουμε και τις προβλέψεις σε 0/1 για την confusion matrix
      pred_full_numeric   <- ifelse(pred_full_labels   == levels(test_data[[label_col]])[2], 1, 0)
      pred_subset_numeric <- ifelse(pred_subset_labels == levels(subset_test[[label_col]])[2], 1, 0)
      
      # Δημιουργία confusion matrix (factor η πρόβλεψη, factor το πραγματικό label)
      cm_full <- confusionMatrix(
        factor(pred_full_numeric,   levels = c(0, 1)),
        factor(test_labels_numeric, levels = c(0, 1))
      )
      cm_subset <- confusionMatrix(
        factor(pred_subset_numeric,   levels = c(0, 1)),
        factor(subset_labels_numeric, levels = c(0, 1))
      )
      
      # Αποθήκευση των metrics
      full_metrics$AUC         <- c(full_metrics$AUC,         auc_full)
      full_metrics$Accuracy    <- c(full_metrics$Accuracy,    cm_full$overall["Accuracy"])
      full_metrics$Sensitivity <- c(full_metrics$Sensitivity, cm_full$byClass["Sensitivity"])
      full_metrics$Specificity <- c(full_metrics$Specificity, cm_full$byClass["Specificity"])
      full_metrics$F1          <- c(full_metrics$F1,          cm_full$byClass["F1"])
      
      subset_metrics$AUC         <- c(subset_metrics$AUC,         auc_subset)
      subset_metrics$Accuracy    <- c(subset_metrics$Accuracy,    cm_subset$overall["Accuracy"])
      subset_metrics$Sensitivity <- c(subset_metrics$Sensitivity, cm_subset$byClass["Sensitivity"])
      subset_metrics$Specificity <- c(subset_metrics$Specificity, cm_subset$byClass["Specificity"])
      subset_metrics$F1          <- c(subset_metrics$F1,          cm_subset$byClass["F1"])
    }
    
    # Υπολογισμός μέσων τιμών metrics και αποθήκευσή τους
    results[[model]] <- list(
      Full = list(
        metrics = lapply(full_metrics, mean, na.rm = TRUE),
        predictions = predictions_full
      ),
      Subset = list(
        metrics = lapply(subset_metrics, mean, na.rm = TRUE),
        predictions = predictions_subset
      )
    )
    
    # Προβολή μέσων metrics στη κονσόλα
    cat(sprintf(
      "Model: %s\nFull Dataset - AUC: %.3f, Accuracy: %.3f, F1: %.3f\nSubset Dataset - AUC: %.3f, Accuracy: %.3f, F1: %.3f\n",
      model, 
      mean(full_metrics$AUC, na.rm = TRUE),
      mean(full_metrics$Accuracy, na.rm = TRUE),
      mean(full_metrics$F1, na.rm = TRUE),
      mean(subset_metrics$AUC, na.rm = TRUE),
      mean(subset_metrics$Accuracy, na.rm = TRUE),
      mean(subset_metrics$F1, na.rm = TRUE)
    ))
  }
  
  # Τερματισμός του parallel cluster
  stopCluster(cl)
  closeAllConnections()
  # Επιστροφή των αποτελεσμάτων
  return(results)
}



# Χρήση της συνάρτησης
results_208_07012025 <- results_MBL_UNB$ga_genes
full_data <- as.data.frame(thca_gse_154763)
#full_data <- full_data
full_data$labels <- as.factor(thca_gse_154763_Labels)

full_data$labels <- relevel(full_data$labels, ref = "healthy")
levels(full_data$labels)

# Εκτέλεση για όλα τα μοντέλα
results_models <- classification_comparison_CV(
  full_data = full_data, 
  selected_genes = results124.20.01.25$ga_genes, 
  label_col = "labels", 
  models = c("RandomForest", "kNN", "NaiveBayes"),
  folds = 5
  
)
# 
# # Check current memory usage
# mem_used <- pryr::mem_used()
# print(mem_used)
# 
# # Monitor memory during operations
# pryr::mem_change(results_models <- classification_comparison_CV(
#   full_data = full_data, 
#   selected_genes = results_208_07012025, 
#   label_col = "labels", 
#   models = c("RandomForest", "kNN", "NaiveBayes"),
#   folds = 5
#   
# )
# )


########### Comparison  Plot
# 1) Υποθέτουμε ότι έχουμε τα αποτελέσματα σε results_models
results <- results_models

# 2) Φτιάχνουμε αρχικά ένα ενιαίο data frame με metrics (Full & Subset)
library(dplyr)
library(tidyr)

results_df <- do.call(rbind, lapply(names(results), function(model) {
  data.frame(
    Model = model,
    Dataset = "Initial Dataset",
    AUC = results[[model]]$Full$metrics$AUC,         # <- πρόσβαση στις μέσες τιμές
    Accuracy = results[[model]]$Full$metrics$Accuracy,
    Sensitivity = results[[model]]$Full$metrics$Sensitivity,
    Specificity = results[[model]]$Full$metrics$Specificity,
    F1 = results[[model]]$Full$metrics$F1
  )
})) %>%
  rbind(
    do.call(rbind, lapply(names(results), function(model) {
      data.frame(
        Model = model,
        Dataset = "Subset",
        AUC = results[[model]]$Subset$metrics$AUC,
        Accuracy = results[[model]]$Subset$metrics$Accuracy,
        Sensitivity = results[[model]]$Subset$metrics$Sensitivity,
        Specificity = results[[model]]$Subset$metrics$Specificity,
        F1 = results[[model]]$Subset$metrics$F1
      )
    }))
  )

# 3) Μετατρέπουμε τα metrics σε long μορφή για εύκολη χρήση με facet
metrics_df <- results_df %>%
  pivot_longer(
    cols = c(AUC, Accuracy, Sensitivity, Specificity, F1),
    names_to = "Metric",
    values_to = "Value"
  )


# 4) Φτιάχνουμε το bar plot
library(ggplot2)
library(viridis)
metrics_df$Dataset[metrics_df$Dataset == "Subset"] <- "GRGS"

### MMC 

library(yardstick)

# Λίστα για αποθήκευση των αποτελεσμάτων MCC
mcc_results <- list()

# Εξαγωγή των μοντέλων και datasets
models <- c("RandomForest", "kNN", "NaiveBayes")
datasets <- c("Full", "Subset")

# Threshold για τις προβλέψεις
threshold <- 0.5

# Επαναληπτικός βρόχος για τα μοντέλα και τα datasets
for (model in models) {
  for (dataset in datasets) {
    # Εξαγωγή labels και προβλέψεων
    labels <- results[[model]][[dataset]]$predictions[[1]]$labels
    predictions <- results[[model]][[dataset]]$predictions[[1]]$predictions
    
    # Μετατροπή προβλέψεων σε κατηγορίες
    predicted_classes <- ifelse(predictions > threshold, "healthy", "diseases")
    
    # Δημιουργία data frame
    df <- data.frame(
      truth = factor(labels, levels = c("diseases", "healthy")),
      estimate = factor(predicted_classes, levels = c("diseases", "healthy"))
    )
    
    # Υπολογισμός MCC
    mcc_metric <- mcc(df, truth = truth, estimate = estimate)
    
    # Αποθήκευση αποτελεσμάτων
    mcc_results[[paste(model, dataset, sep = "_")]] <- mcc_metric
  }
}


# Μετατροπή της λίστας mcc_results σε dataframe
mcc_df <- do.call(rbind, lapply(names(mcc_results), function(name) {
  model_dataset <- unlist(strsplit(name, "_"))
  data.frame(
    Model = model_dataset[1],
    Dataset = ifelse(model_dataset[2] == "Full", "Initial Dataset", "GRGS"),
    Metric = "MCC",
    Value = mcc_results[[name]]$.estimate
  )
}))

# Συνένωση των δεδομένων MCC με το metrics_df
metrics_df <- rbind(metrics_df, mcc_df)


ggplot(metrics_df, aes(x = Model, y = Value, fill = Dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Metric, scales = "free_y", nrow = 2) +
  labs(
    title = "Performance Metrics of Machine Learning Models: Initial Dataset vs. GRGS Subset",
    x = "Model",
    y = "Metric Value",
    fill = "Data Type"
  )+
  scale_fill_manual(values = c("Initial Dataset" = "#1b9e77", "GRGS" = "#d95f02")) + # Συνεπής χρωματική παλέτα
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Προσθήκη επιπλέον χώρου
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(size = 14)
  ) +
  geom_text(
    aes(label = round(Value, 2)),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5
  )

ggsave("Performance_Metrics_Machine_Learning_Models_White.png", dpi = 300, width = 10, height = 8, bg = "white")




##### ROC / / /  AUC


# Υπολογισμός AUC
calculate_auc_from_results <- function(results) {
  library(pROC)

  auc_results <- data.frame(
    Model = character(),
    Dataset = character(),
    AUC = numeric(),
    stringsAsFactors = FALSE
  )

  for (model in names(results)) {
    if (!is.null(results[[model]]$Full$predictions)) {
      all_labels <- unlist(lapply(results[[model]]$Full$predictions, function(x) x$labels))
      all_predictions <- unlist(lapply(results[[model]]$Full$predictions, function(x) x$predictions))
      auc_full <- roc(all_labels, all_predictions)$auc
      auc_results <- rbind(auc_results, data.frame(Model = model, Dataset = "Full", AUC = auc_full))
    }
  }

  return(auc_results)}

# Υπολογισμός AUC για το μοντέλο GBM
auc_single_model <- calculate_auc_from_results(results_models)
print(auc_single_model)

# 
# Συνάρτηση για δημιουργία δεδομένων ROC
generate_roc_df <- function(labels, predictions, model_name, dataset_name) {
  roc_obj <- roc(labels, predictions)
  data.frame(
    Specificity = rev(roc_obj$specificities),
    Sensitivity = rev(roc_obj$sensitivities),
    Model = model_name,
    Dataset = dataset_name
  )
}

# Συνδυασμός όλων των ROC δεδομένων
combine_roc_data <- function(results) {
  roc_combined <- data.frame()
  
  for (model in names(results)) {
    # Full Dataset
    if (!is.null(results[[model]]$Full$predictions)) {
      all_labels <- unlist(lapply(results[[model]]$Full$predictions, function(x) x$labels))
      all_predictions <- unlist(lapply(results[[model]]$Full$predictions, function(x) x$predictions))
      roc_full <- generate_roc_df(all_labels, all_predictions, model, "Full")
      roc_combined <- rbind(roc_combined, roc_full)
    }
    
    # Subset Dataset
    if (!is.null(results[[model]]$Subset$predictions)) {
      all_labels <- unlist(lapply(results[[model]]$Subset$predictions, function(x) x$labels))
      all_predictions <- unlist(lapply(results[[model]]$Subset$predictions, function(x) x$predictions))
      roc_subset <- generate_roc_df(all_labels, all_predictions, model, "Subset")
      roc_combined <- rbind(roc_combined, roc_subset)
    }
  }
  
  return(roc_combined)
}

# Δημιουργία δεδομένων ROC από τα αποτελέσματα
roc_combined <- combine_roc_data(results_models)

library(ggplot2)

roc_combined$Dataset[roc_combined$Dataset == "Subset"] <- "GRGS"
roc_combined$Dataset[roc_combined$Dataset == "Full"] <- "Initial Dataset"


# Δημιουργία του γραφήματος ROC
ggplot(roc_combined, aes(x = 1 - Specificity, y = Sensitivity, color = Dataset, linetype = Dataset)) +
  geom_line(size = 1) +
  facet_wrap(~ Model, scales = "free", ncol = 3) + # Ρυθμίζει τα panel για κάθε μοντέλο (3 σε μία σειρά)
  labs(
    title = "ROC Analysis: Initial Dataset vs. GRGS Subset Across Machine Learning Models",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = "Data Type",
    linetype = "Data Type"
  ) +
  scale_color_manual(values = c("Initial Dataset" = "#1b9e77", "GRGS" = "#d95f02")) + # Χρώματα για τα dataset
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(size = 14)
  )+
  guides(color = guide_legend(title = "Data Type"), linetype = guide_legend(title = "Data Type"))



# Αποθήκευση του γραφήματος σε υψηλή ανάλυση

ggsave("ROC_Curves_for_Paper.png", dpi = 300, width = 10, height = 8, bg = "white")
