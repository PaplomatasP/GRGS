plot_evaluation <- function(results) {
  
  
  # 1. Βγάζεις όσες γενιές είναι πράγματι λίστα
  valid_idx <- sapply(results$results, is.list)
  good_list <- results$results[valid_idx]
  
  # 2. Φτιάχνεις το data frame
  results_df <- do.call(rbind, lapply(good_list, function(x) {
    data.frame(
      Generation   = x$generation,
      Best_Fitness = x$best_fitness,
      Mean_Fitness = x$mean_fitness,
      Train_Error  = x$train_error,
      Q_LR         = x$Q_LR,
      logFC        = x$logFC
    )
  }))
  
  # 3. Βιβλιοθήκες
  library(ggplot2)
  library(patchwork)
  cutoff_generation <-  which.max(results_df$Best_Fitness)
  cutoff_fitness <- results_df$Best_Fitness[cutoff_generation]
  
  # 4. Δημιουργία plots
  # Best Fitness Plot
  p1 <- ggplot(results_df, aes(x = Generation, y = Best_Fitness)) +
    geom_line(color = "blue", size = 1.2) +
    labs(title = "Best Fitness", x = "Generation", y = "Fitness") +
    theme_minimal() +
    geom_vline(
      xintercept = cutoff_generation,
      color = "red",
      linetype = "twodash",
      size = 1
    ) +
    # Προσθήκη κειμένου στο cutoff
    annotate(
      "text",
      x = cutoff_generation,
      y = cutoff_fitness,
      label =paste("Optimal Generation:", cutoff_generation),
      color = "red",
      hjust = -0.1,
      vjust = 1.2,
      size = 4
    ) +
    theme(plot.title = element_text(hjust = 0.5)) # Κεντραρισμένος τίτλος 
  
  # Mean Fitness Plot
  p2 <- ggplot(results_df, aes(x = Generation, y = Mean_Fitness)) +
    geom_line(color = "orange", size = 1.2) +
    labs(title = "Mean Fitness", x = "Generation", y = "Fitness") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) # Κεντραρισμένος τίτλος 
  
  # Train Error Plot
  p3 <- ggplot(results_df, aes(x = Generation, y = Train_Error)) +
    geom_line(color = "red", size = 1.2) +
    labs(title = "Train Error", x = "Generation", y = "Train Error") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) # Κεντραρισμένος τίτλος 
  
  # Q_LR Evolution Plot
  p4 <- ggplot(results_df, aes(x = Generation, y = Q_LR)) +
    geom_line(color = "purple", size = 1.2) +
    labs(title = "Q_LR Evolution", x = "Generation", y = "Q_LR") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) # Κεντραρισμένος τίτλος 
  
  p5 <- ggplot(results_df, aes(x = Generation, y = logFC)) +
    geom_line(color = "darkgreen", size = 1.2) +
    #geom_point(color = "darkgreen", size = 2) + # Προσθήκη σημείων
    #geom_smooth(method = "loess", color = "black", linetype = "dashed") + # Εξομάλυνση
    labs(title = "logFC Evolution",
         x = "Generation", y = "Median logFC (normalized)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5)) #+
  # geom_vline(
  #   xintercept = cutoff_generation,
  #   color = "red",
  #   linetype = "dotted",
  #   size = 1
  # )
  
  # 5. Συνδυασμός plots με patchwork
  combined_plot <- p1 / (p2 | p3) / (p4 | p5)
  
  # Εμφάνιση
  print(combined_plot)
}
results = results124.20.01.25

# Πρόσθεσε το flag evaluationPlot για να κάνεις την οπτικοποίηση
plot_evaluation(results)

ggsave("FitnessEvolution.png", dpi = 300, width = 10, height = 8, bg = "white")
