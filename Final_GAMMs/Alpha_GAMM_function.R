#==================================================
# LIBRARIES
#==================================================

library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(mgcv)
library(ggplot2)
library(glue)
library(patchwork)
library(gratia)
library(hillR)
library(tibble)
library(glue)

#==================================================
# MAIN PIPELINE FUNCTION
#==================================================

run_alpha_diversity_pipeline <- function(taxa_name,
                                         base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2",
                                         save_figs = TRUE) {
  
  #-----------------------------
  # FILE PATHS
  #-----------------------------
  
  taxa_path <- file.path(base_path, "taxa_lists", "cleaned",
                         paste0(taxa_name, "_genus_relative_abundance.tsv"))
  
  meta_path <- file.path(base_path, "Metadata", "23S_LEAP22_metadata_NR.tsv")
  
  output_dir  <- file.path(base_path, "figures/alpha_diversity", paste0(taxa_name))
  
  #-----------------------------
  # LOAD DATA
  #-----------------------------
  
  taxa <- fread(taxa_path)
  meta <- fread(meta_path)
  
  taxa[is.na(taxa)] <- 0
  
  #-----------------------------
  # HILL NUMBERS
  #-----------------------------
  
  # Remove T1
  meta <- meta %>%
    filter(time != "T1")
  
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
  # RUN ALPHA OVERTIME
  #-----------------------------
  
  # Colour scale
  nutrient_cols <- c(
    "0" = "#440154FF",
    "25" = "#39568CFF",
    "50" = "#20A387FF",
    "100" = "#73D055FF",
    "200" = "#FDE725FF"
  )
  
  
  #-----------------------------
  # Prepare plotting data
  #-----------------------------
  
  alpha_plot_data <- hill_results %>%
    
    mutate(
      
      Movement = factor(
        as.character(flow_rate),
        levels = c("0","10","20","30","40")
      ),
      
      Node = factor(
        upstream_downstream,
        levels = c("Up","Down"),
        labels = c("Upstream","Downstream")
      ),
      
      NutrientP = factor(
        nutrient_addition,
        levels = c(0,25,50,100,200)
      )
      
    ) %>%
    
    group_by(
      Movement,
      NutrientP,
      Node,
      day
    ) %>%
    
    summarise(
      
      q0 = mean(q0, na.rm = TRUE),
      q1 = mean(q1, na.rm = TRUE),
      q2 = mean(q2, na.rm = TRUE),
      
      .groups = "drop"
      
    )
  
  
  #-----------------------------
  # Plot function
  #-----------------------------
  
  plot_alpha <- function(response, ylab){
    
    ggplot(
      alpha_plot_data,
      aes(
        x = day,
        y = .data[[response]],
        colour = NutrientP,
        group = NutrientP
      )
    ) +
      
      geom_line(
        linewidth = 1,
        na.rm = TRUE
      ) +
      
      geom_point(
        size = 2,
        na.rm = TRUE
      ) +
      
      facet_grid(
        rows = vars(Node),
        cols = vars(Movement),
        switch = "y"
      ) +
      
      scale_colour_manual(
        values = nutrient_cols,
        name = "Nutrient Enrichment"
      ) +
      
      scale_x_continuous(
        breaks = c(0,28,49,73,91)
      ) +
      
      labs(
        x = "Day",
        y = ylab
      ) +
      
      theme_classic() +
      
      theme(
        
        legend.position = "bottom",
        
        strip.background = element_blank(),
        
        strip.placement = "outside",
        
        strip.text.x = element_text(
          face = "bold",
          size = 13
        ),
        
        strip.text.y.right = element_text(
          face = "bold",
          size = 13,
          angle = -90
        ),
        
        panel.spacing = unit(1,"lines")
        
      )
    
  }
  
  
  #-----------------------------
  # Create plots
  #-----------------------------
  
  p_q0 <- plot_alpha(
    "q0",
    "Hill diversity (q0)"
  )
  
  p_q1 <- plot_alpha(
    "q1",
    "Hill diversity (q1)"
  )
  
  p_q2 <- plot_alpha(
    "q2",
    "Hill diversity (q2)"
  )
  
  
  #-----------------------------
  # Save plots
  #-----------------------------
  
  ggsave(
    filename = file.path(output_dir,
                         paste0(taxa_name,"_alpha_q0.png")),
    plot = p_q0,
    width = 16,
    height = 7,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(output_dir,
                         paste0(taxa_name,"_alpha_q1.png")),
    plot = p_q1,
    width = 16,
    height = 7,
    dpi = 600
  )
  
  ggsave(
    filename = file.path(output_dir,
                         paste0(taxa_name,"_alpha_q2.png")),
    plot = p_q2,
    width = 16,
    height = 7,
    dpi = 600
  )
  
  #-----------------------------
  # PROCESS FOR GAMMS
  #-----------------------------
  
  colnames(meta)
  
  colnames(hill_results)
  
  
  # Clean Data for GAMM
  alpha <- hill_results %>%
    mutate(
      nutrient_addition = as.numeric(nutrient_addition),
      flow_rate = as.numeric(flow_rate),
      day = as.numeric(day),
      upstream_downstream = as.factor(upstream_downstream),
      
      # GAMM variables
      Movement = flow_rate,
      NutrientP = nutrient_addition,
      Day_num = day,
      Node = upstream_downstream,
      
      # random effect
      Mesocosm = as.factor(tank)
    )
  
  #Remove empty samples
  empty_samples <- alpha %>%
    filter(
      q0 == 0 | is.na(q1) | is.na(q2)
    ) %>%
    select(sample_name, time, flow_rate, nutrient_addition, upstream_downstream)
  
  
  print(empty_samples)
  
  #-----------------------------
  # RUN GAMMS
  #-----------------------------
  
  if (grepl("^18S", taxa_name)) {
    
    alpha <- alpha %>%
      mutate(
        q0 = log1p(q0),
        q1 = log1p(q1),
        q2 = log1p(q2)
      )
    
  } else {
    
    alpha <- alpha %>%
      mutate(
        q0 = sqrt(q0),
        q1 = log1p(q1),
        q2 = log1p(q2)
      )
  }
  
  model_list <- list()
  deviance_table <- data.frame()
  
  for (hill in c("q0", "q1", "q2")) {
    
    k_day <- 4
    k <- 4
    
    formula_text <- glue(
      "{hill} ~
        Node +
        ti(Day_num, by = Node, k = {k_day}) +
        ti(Movement, by = Node, k = {k}) +
        ti(NutrientP, by = Node, k = {k}) +
        ti(Movement, Day_num, by = Node, k = c({k},{k_day})) +
        ti(NutrientP, Day_num, by = Node, k = c({k},{k_day})) +
        ti(Movement, NutrientP, by = Node, k = c({k},{k})) +
        ti(Movement, NutrientP, Day_num, by = Node, k = c({k},{k},{k_day})) +
        s(Mesocosm, bs = 're')"
    )
    
    model <- gam(as.formula(formula_text),
                 data = alpha,
                 method = "REML")
    
    model_list[[hill]] <- model
    
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
      "Node" = "Node",
      "ti(Day_num):NodeUp" =
        "ti(Day):Upstream",
      "ti(Day_num):NodeDown" =
        "ti(Day):Downstream",
      "ti(Movement):NodeUp" =
        "ti(Connectivity):Upstream",
      "ti(Movement):NodeDown" =
        "ti(Connectivity):Downstream",
      "ti(NutrientP):NodeUp" =
        "ti(Nutrient Enrichment):Upstream",
      "ti(NutrientP):NodeDown" =
        "ti(Nutrient Enrichment):Downstream",
      "ti(Movement,Day_num):NodeUp" =
        "ti(Connectivity, Day):Upstream",
      "ti(Movement,Day_num):NodeDown" =
        "ti(Connectivity, Day):Downstream",
      "ti(NutrientP,Day_num):NodeUp" =
        "ti(Nutrient Enrichment, Day):Upstream",
      "ti(NutrientP,Day_num):NodeDown" =
        "ti(Nutrient Enrichment, Day):Downstream",
      "ti(Movement,NutrientP):NodeUp" =
        "ti(Connectivity, Nutrient Enrichment):Upstream",
      "ti(Movement,NutrientP):NodeDown" =
        "ti(Connectivity, Nutrient Enrichment):Downstream",
      "ti(Movement,NutrientP,Day_num):NodeUp" =
        "ti(Connectivity, Nutrient Enrichment, Day):Upstream",
      "ti(Movement,NutrientP,Day_num):NodeDown" =
        "ti(Connectivity, Nutrient Enrichment, Day):Downstream",
      "s(Mesocosm)" =
        "s(Mesocosm)"
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
          paste0(taxa_name, "_alpha_", hill, "_GAM_table.tsv")
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
        paste0(taxa_name, "_alpha_deviance_explained.tsv")
      ),
      sep = "\t"
    )
    
    #-----------------------------
    # SAVE GAMs Themselves
    #-----------------------------
    
    saveRDS(
      model_list[["q0"]],
      file.path(output_dir, paste0(taxa_name, "_alpha_q0.rds"))
    )
    
    saveRDS(
      model_list[["q1"]],
      file.path(output_dir, paste0(taxa_name, "_alpha_q1.rds"))
    )
    
    saveRDS(
      model_list[["q2"]],
      file.path(output_dir, paste0(taxa_name, "_alpha_q2.rds"))
    )
    
    #-----------------------------
    # PREDICTION GRAPHS
    #-----------------------------
    
    time_points <- c(28, 49, 73, 91)
    
    pred_grid <- expand.grid(
      Movement = seq(
        min(alpha$Movement),
        max(alpha$Movement),
        length.out = 200
      ),
      NutrientP = seq(
        min(alpha$NutrientP),
        max(alpha$NutrientP),
        length.out = 200
      ),
      Day_num = time_points,
      Node = levels(model.frame(model_list[[1]])$Node)
    )
    
    # Downstream ponds do not exist below 10% connectivity
    pred_grid <- pred_grid %>%
      filter(!(Node == "Down" & Movement < 10))
    
    pred_grid$Mesocosm <- levels(alpha$Mesocosm)[1]
    
    pred_grid$Time <- factor(
      pred_grid$Day_num,
      levels = time_points
    )
    
    #==================================================
    # SURFACE PLOTS
    #==================================================
    
    surface_plots <- list()
    
    for(hill in c("q0","q1","q2")){
      
      # Keep model input unchanged
      plot_df <- pred_grid
      
      #--------------------------------------------------
      # Predictions
      #--------------------------------------------------
      
      plot_df$Prediction <- predict(
        model_list[[hill]],
        newdata = plot_df,
        type = "response",
        exclude = "s(Mesocosm)"
      )
      
      #--------------------------------------------------
      # Undo transformations
      #--------------------------------------------------
      
      if(hill == "q0"){
        
        if(grepl("^18S", taxa_name)){
          
          plot_df$Prediction <- exp(plot_df$Prediction) - 1
          
        } else{
          
          plot_df$Prediction <- plot_df$Prediction^2
          
        }
        
      } else{
        
        plot_df$Prediction <- exp(plot_df$Prediction) - 1
        
      }
      
      
      #--------------------------------------------------
      # Create colour bins
      #--------------------------------------------------
      
      contour_breaks <- seq(
        min(plot_df$Prediction, na.rm = TRUE),
        max(plot_df$Prediction, na.rm = TRUE),
        length.out = 11
      )
      
      plot_df$PredictionBin <- cut(
        plot_df$Prediction,
        breaks = contour_breaks,
        include.lowest = TRUE,
        ordered_result = TRUE
      )
      
      
      #--------------------------------------------------
      # Create plotting copy ONLY
      #--------------------------------------------------
      
      plot_plot <- plot_df
      
      plot_plot$Node <- factor(
        plot_plot$Node,
        levels = c("Up","Down"),
        labels = c("Upstream","Downstream")
      )
      
      
      #--------------------------------------------------
      # Plot
      #--------------------------------------------------
      
      p <- ggplot(
        plot_plot,
        aes(x = NutrientP, y = Movement)
      ) +
        
        geom_raster(aes(fill = PredictionBin)) +
        
        stat_contour(
          aes(z = Prediction),
          breaks = contour_breaks,
          colour = "black",
          linewidth = 0.4
        ) +
        
        facet_grid(Node ~ Time) +
        
        scale_fill_viridis_d(
          option = "viridis",
          name = "Predicted Diversity",
          drop = FALSE
        ) +
        
        scale_x_continuous(
          breaks = scales::pretty_breaks(n = 3), 
          expand = expansion(mult = 0.02)    
        ) +
        
        scale_y_continuous(
          breaks = scales::pretty_breaks(n = 4),
          expand = expansion(mult = 0.02)
        ) +
        
        labs(
          x = expression("Nutrient Enrichment (" * mu * "g L"^{-1} * " week"^{-1} * ")"),
          y = expression("Connectivity (% week"^{-1} * ")")
        ) +
        
        theme_bw() +
        
        theme(
          panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey90"),
          strip.text = element_text(face = "bold"),
          panel.spacing.x = unit(0.6, "lines"),
          panel.spacing.y = unit(0.6, "lines"),
          axis.text.x = element_text(size = 8),
          axis.text.y = element_text(size = 8)
        )
      
      surface_plots[[hill]] <- p
      
      if(save_figs){
        
        ggsave(
          file.path(
            output_dir,
            paste0(hill, "_surface_plot.png")
          ),
          p,
          width = 12,
          height = 7,
          dpi = 600
        )
        
      }
    }
}

#==================================================
# RUN PIPELINE FUNCTION
#==================================================

result_23S <- run_alpha_diversity_pipeline("23S_total")
result_23S_green_algae <- run_alpha_diversity_pipeline("23S_green_algae")
result_23S_cyanobacteria <- run_alpha_diversity_pipeline("23S_cyanobacteria")

result_18S <- run_alpha_diversity_pipeline("18S_total")
result_18S_phytoplankton <- run_alpha_diversity_pipeline("18S_phytoplankton")
result_18S_green_algae <- run_alpha_diversity_pipeline("18S_green_algae")
result_18S_zooplankton <- run_alpha_diversity_pipeline("18S_zooplankton")
result_18S_micro_heterotrophs <- run_alpha_diversity_pipeline("18S_micro_heterotrophs")
result_18S_fungi <- run_alpha_diversity_pipeline("18S_fungi")

result_COI <- run_alpha_diversity_pipeline("COI_total")
result_COI_zooplankton <- run_alpha_diversity_pipeline("COI_zooplankton")
  
  
  
  
  
  
  
  
  
  





