#==================================================
# LIBRARIES
#==================================================

library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(hillR)
library(tibble)
library(mgcv)
library(glue)

#==================================================
# MAIN PIPELINE FUNCTION
#==================================================


run_gamma_diversity_pipeline <- function(taxa_name,
                                         base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2",
                                         save_figs = TRUE) {
  
  #-----------------------------
  # FILE PATHS
  #-----------------------------
  
  taxa_path <- file.path(base_path, "taxa_lists", "gamma_diversity",
                         paste0(taxa_name, "_gamma_relative_abundance.tsv"))
  
  meta_path <- file.path(base_path, "taxa_lists", "gamma_diversity",
                         paste0(taxa_name, "_gamma_metadata.tsv"))
  
  output_dir  <- file.path(base_path, "figures/gamma_diversity", taxa_name)
  
  #-----------------------------
  # LOAD DATA
  #-----------------------------
  
  taxa <- fread(taxa_path)
  meta <- fread(meta_path)
  
  #-----------------------------
  # HILL NUMBERS
  #-----------------------------
  
  # Extract abundance matrix
  comm <- taxa %>%
    column_to_rownames("Genus") %>%
    as.matrix()
  
  # Convert percentages back to proportions
  comm <- comm / 100
  
  # Transpose so samples are rows
  comm <- t(comm)
  
  # Remove samples with zero total abundance
  empty_samples <- rownames(comm)[rowSums(comm) == 0]
  
  if(length(empty_samples) > 0){
    
    message(
      "Removing empty samples:\n",
      paste(empty_samples, collapse = "\n")
    )
    
    comm <- comm[rowSums(comm) > 0, ]
  }
  
  # Calculate Hill numbers
  hill_results <- data.frame(
    sample_name = rownames(comm),
    q0 = hill_taxa(comm, q = 0),
    q1 = hill_taxa(comm, q = 1),
    q2 = hill_taxa(comm, q = 2)
  )
  
  # Combine with metadata
  hill_results <- meta %>%
    left_join(hill_results, by = "sample_name")
  
  # Add day columngamma_diversity <- gamma_diversity %>%
  hill_results <- hill_results %>%
    mutate(
      day = case_when(
        time == "T2" ~ 28,
        time == "T3" ~ 49,
        time == "T4" ~ 73,
        time == "T5" ~ 91,
        TRUE ~ NA_real_
      )
    )
  
  # Put all hill numbers into one column
  hill_long <- hill_results %>%
    pivot_longer(
      cols = c(q0, q1, q2),
      names_to = "q",
      values_to = "gamma"
    )
  
  #-----------------------------
  # PLOT GAMMA OVERTIME
  #-----------------------------
  
  # Colour scale
  nutrient_cols <- c(
    "0" = "#440154FF",
    "25" = "#39568CFF",
    "50" = "#20A387FF",
    "100" = "#73D055FF",
    "200" = "#FDE725FF"
  )
  
  # Average biological replicates
  gamma_mean <- hill_long %>%
    group_by(time, day, nutrient_addition, flow_rate, q) %>%
    summarise(
      gamma = mean(gamma, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      nutrient_addition = factor(
        nutrient_addition,
        levels = c(0, 25, 50, 100, 200)
      )
    )
  
  # Plot function
  plot_gamma <- function(q_label, output_file) {
    
    df <- gamma_mean %>%
      filter(q == q_label)
    
    p <- ggplot(
      df,
      aes(
        x = day,
        y = gamma,
        color = nutrient_addition,
        group = nutrient_addition
      )
    ) +
      
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      
      facet_wrap(~ flow_rate, nrow = 1, scales = "fixed") +
      
      scale_color_manual(
        values = nutrient_cols,
        name = "Nutrient Addition"
      ) +
      
      labs(
        x = "Day",
        y = paste("Gamma diversity", q_label)
      ) +
      
      theme_classic() +
      
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      )
    
    ggsave(
      filename = output_file,
      plot = p,
      width = 14,
      height = 4,
      units = "in",
      dpi = 600
    )
    
    return(p)
  }
  
  # Create output 
  if(save_figs){
    
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    output_plot_q0 <- file.path(
      output_dir,
      paste0(taxa_name, "_gamma_q0.png")
    )
    
    output_plot_q1 <- file.path(
      output_dir,
      paste0(taxa_name, "_gamma_q1.png")
    )
    
    output_plot_q2 <- file.path(
      output_dir,
      paste0(taxa_name, "_gamma_q2.png")
    )
  }
  
  # Run and save plots
  p0 <- plot_gamma("q0", output_plot_q0)
  p1 <- plot_gamma("q1", output_plot_q1)
  p2 <- plot_gamma("q2", output_plot_q2)
  
  #-----------------------------
  # RUN GAMM
  #-----------------------------
  
  # Clean Data for GAMM
  gamma <- hill_results %>%
    mutate(
      nutrient_addition = as.numeric(nutrient_addition),
      flow_rate = as.numeric(flow_rate),
      day = as.numeric(day),
      
      # GAMM variables
      Movement = flow_rate,
      NutrientP = nutrient_addition,
      Day_num = day,
      
      # random effect
      Mesocosm = as.factor(tank)
    )
  
  #Remove empty samples
  empty_samples <- gamma %>%
    filter(
      q0 == 0 | is.na(q1) | is.na(q2)
    ) %>%
    select(sample_name, time, flow_rate, nutrient_addition)
  
  
  print(empty_samples)
  
  gamma <- gamma %>%
    filter(q0 > 0, !is.na(q1), !is.na(q2))
  
  #Run GAM
  model_list <- list()
  
  deviance_table <- data.frame()
  
  for (hill in c("q0", "q1", "q2")) {
    
    k_day <- 4
    k <- 4
    
    formula_text <- glue(
      "{hill} ~
        ti(Day_num, k = {k_day}) +
        ti(Movement, k = {k}) +
        ti(NutrientP, k = {k}) +
        ti(Movement, Day_num, k = c({k},{k_day})) +
        ti(NutrientP, Day_num, k = c({k},{k_day})) +
        ti(Movement, NutrientP, k = c({k},{k})) +
        ti(Movement, NutrientP, Day_num, k = c({k},{k},{k_day})) +
        s(Mesocosm, bs = 're')"
    )
    
    
    model <- gam(
      as.formula(formula_text),
      data = gamma,
      method = "REML",
      family = Gamma(link = "log")
    )
    
    
    model_list[[hill]] <- model
    
    # -----------------------------
    # GAM CHECK (diagnostics)
    # -----------------------------
    
    # Save the 4-panel diagnostic plots
    png(
      file.path(output_dir, paste0(taxa_name, "_gamma_", hill, "_gamcheck.png")),
      width = 8, height = 8, units = "in", res = 300
    )
    gam.check(model)
    dev.off()
    
    # Save the printed convergence + k-index text
    check_text <- capture.output(gam.check(model))
    writeLines(
      check_text,
      file.path(output_dir, paste0(taxa_name, "_gamma_", hill, "_gamcheck.txt"))
    )
    
    # Save the k-index table itself as a clean data frame
    k_table <- as.data.frame(k.check(model))
    k_table$Term <- rownames(k_table)
    rownames(k_table) <- NULL
    
    fwrite(
      k_table,
      file.path(output_dir, paste0(taxa_name, "_gamma_", hill, "_kcheck.tsv")),
      sep = "\t"
    )
    
    # Deviance explained
    deviance_table <- rbind(
      deviance_table,
      data.frame(
        Hill = hill,
        Deviance_explained = summary(model)$dev.expl
      )
    )
    
    print(summary(model))
    
  }
  
  #-----------------------------
  # SAVE GAM TABLES
  #-----------------------------
  
  # Function to format p-values
  format_p <- function(x){
    
    if(x < 0.001){
      return("<0.001")
    } else if(x < 0.01){
      return("<0.01")
    } else {
      return(format(round(x, 2), nsmall = 2))
    }
  }
  
  
  # Desired term names and order
  term_names <- c(
    "ti(Day_num)" = "ti(Day)",
    "ti(Movement)" = "ti(Connectivity)",
    "ti(NutrientP)" = "ti(Nutrient Enrichment)",
    "ti(Movement,Day_num)" = "ti(Connectivity, Day)",
    "ti(NutrientP,Day_num)" = "ti(Nutrient Enrichment, Day)",
    "ti(Movement,NutrientP)" = "ti(Connectivity, Nutrient Enrichment)",
    "ti(Movement,NutrientP,Day_num)" = "ti(Connectivity, Nutrient Enrichment, Day)",
    "s(Mesocosm)" = "s(PairID)"
  )
  
  
  # Extract each GAM table
  for(hill in c("q0", "q1", "q2")){
    
    model <- model_list[[hill]]
    
    gam_table <- as.data.frame(summary(model)$s.table)
    
    # Move rownames into column
    gam_table$Term <- rownames(gam_table)
    rownames(gam_table) <- NULL
    
    # Rename columns
    gam_table <- gam_table %>%
      select(
        Term,
        edf = edf,
        F = `F`,
        p_value = `p-value`
      )
    
    # Rename terms
    gam_table$Term <- term_names[gam_table$Term]
    
    # Reorder terms
    gam_table$Term <- factor(
      gam_table$Term,
      levels = unname(term_names)
    )
    
    gam_table <- gam_table %>%
      arrange(Term)
    
    # Format numbers
    gam_table <- gam_table %>%
      mutate(
        edf = round(edf, 2),
        F = round(F, 2),
        p_value = sapply(p_value, format_p)
      ) %>%
      mutate(
        Term = as.character(Term)
      )
    
    
    # Save table
    fwrite(
      gam_table,
      file.path(
        output_dir,
        paste0(taxa_name, "_gamma_", hill, "_GAM_table.tsv")
      ),
      sep = "\t"
    )
  }
  
  
  #-----------------------------
  # SAVE DEVIANCE EXPLAINED
  #-----------------------------
  
  deviance_table <- deviance_table %>%
    mutate(
      Deviance_explained = round(Deviance_explained * 100, 2)
    )
  
  
  fwrite(
    deviance_table,
    file.path(
      output_dir,
      paste0(taxa_name, "_gamma_deviance_explained.tsv")
    ),
    sep = "\t"
  )
  
  #-----------------------------
  # SAVE GAMs Themselves
  #-----------------------------
  
  saveRDS(
    model_list[["q0"]],
    file.path(output_dir, paste0(taxa_name, "_gamma_q0.rds"))
  )
  
  saveRDS(
    model_list[["q1"]],
    file.path(output_dir, paste0(taxa_name, "_gamma_q1.rds"))
  )
  
  saveRDS(
    model_list[["q2"]],
    file.path(output_dir, paste0(taxa_name, "_gamma_q2.rds"))
  )
  
  #-----------------------------
  # PREDICTION GRAPHS
  #-----------------------------
  
  time_points <- c(28, 49, 73, 91)
  
  plot_list <- list()
  
  
  for (hill in c("q0", "q1", "q2")) {
    
    
    #-----------------------------------------------
    # Prediction grid
    #-----------------------------------------------
    
    pred_grid <- expand.grid(
      
      Movement = seq(min(gamma$Movement),
                     max(gamma$Movement),
                     length.out = 150),
      
      NutrientP = seq(min(gamma$NutrientP),
                      max(gamma$NutrientP),
                      length.out = 150),
      
      Day_num = time_points,
      
      Mesocosm = levels(gamma$Mesocosm)[1]
    )
    
    
    pred_grid$Time <- factor(
      pred_grid$Day_num,
      levels = time_points
    )
    
    
    #-----------------------------------------------
    # Predictions
    #-----------------------------------------------
    
    pred_grid$Prediction <- predict(
      model_list[[hill]],
      newdata = pred_grid,
      type = "response",
      exclude = "s(Mesocosm)"
    )
    
    #-----------------------------------------------
    # Contour bins
    #-----------------------------------------------
    
    contour_breaks <- seq(
      min(pred_grid$Prediction, na.rm = TRUE),
      max(pred_grid$Prediction, na.rm = TRUE),
      length.out = 11
    )
    
    
    pred_grid$PredictionBin <- cut(
      pred_grid$Prediction,
      breaks = contour_breaks,
      include.lowest = TRUE,
      ordered_result = TRUE
    )
    
    
    #-----------------------------------------------
    # Plot
    #-----------------------------------------------
    
    p <- ggplot(pred_grid,
                aes(x = NutrientP,
                    y = Movement)) +
      
      geom_raster(
        aes(fill = PredictionBin)
      ) +
      
      geom_contour(
        aes(z = Prediction),
        breaks = contour_breaks,
        colour = "black",
        linewidth = 0.4
      ) +
      
      facet_wrap(~Time,
                 nrow = 1) +
      
      scale_fill_viridis_d(
        option = "viridis",
        name = "Predicted Gamma Diversity"
      ) +
      
      labs(
        x = "Nutrient Enrichment (µg L⁻¹ week⁻¹)",
        y = "Connectivity (% week⁻¹)"
      ) +
      
      theme_bw() +
      
      theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold")
      )
    
    
    plot_list[[hill]] <- p
    
  }
  
  
  #==================================================
  # SAVE ALL PLOTS
  #==================================================
  
  for (hill in names(plot_list)) {
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0(taxa_name, "_", hill, "_surface.png")
      ),
      plot = plot_list[[hill]],
      width = 12,
      height = 5,
      dpi = 600
    )
    
  }
  
}

#==================================================
# RUN PIPELINE FUNCTION
#==================================================

result_23S <- run_gamma_diversity_pipeline("23S_total")
result_23S_green_algae <- run_gamma_diversity_pipeline("23S_green_algae")
result_23S_cyanobacteria <- run_gamma_diversity_pipeline("23S_cyanobacteria")

result_18S <- run_gamma_diversity_pipeline("18S_total")
result_18S_phytoplankton <- run_gamma_diversity_pipeline("18S_phytoplankton")
result_18S_green_algae <- run_gamma_diversity_pipeline("18S_green_algae")
result_18S_zooplankton <- run_gamma_diversity_pipeline("18S_zooplankton")
result_18S_micro_heterotrophs <- run_gamma_diversity_pipeline("18S_micro_heterotrophs")
result_18S_fungi <- run_gamma_diversity_pipeline("18S_fungi")

result_COI <- run_gamma_diversity_pipeline("COI_total")
result_COI_zooplankton <- run_gamma_diversity_pipeline("COI_zooplankton")




