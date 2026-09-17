#==================================================
# LIBRARIES
#==================================================

library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(ape)
library(ggplot2)
library(hillR)
library(tibble)
library(cowplot)
library(patchwork)
library(grid)
library(glue)
library(mgcv)

#==================================================
# MAIN PIPELINE FUNCTION
#==================================================

run_beta_diversity_pipeline <- function(
    taxa_name,
    base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2",
    save_figs = TRUE
) {
  
  
  #-----------------------------
  # FILE PATHS
  #-----------------------------
  
  taxa_path <- file.path(
    base_path,
    "taxa_lists",
    "cleaned",
    paste0(taxa_name, "_genus_relative_abundance.tsv")
  )
  
  meta_path <- file.path(
    base_path,
    "Metadata",
    "23S_LEAP22_metadata_NR.tsv"
  )
  
  output_dir <- file.path(
    base_path,
    "figures",
    "beta_diversity",
    taxa_name
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  #-----------------------------
  # LOAD DATA
  #-----------------------------
  
  taxa <- fread(taxa_path)
  meta <- fread(meta_path)
  
  # Remove T1
  meta <- meta %>%
    filter(time != "T1")
  
  # Save genus names
  genus_names <- taxa$Genus
  
  # Remove genus column
  taxa <- taxa[, -1, with = FALSE]
  
  # Transpose
  taxa <- as.data.frame(t(taxa))
  
  # Name the columns by genus
  colnames(taxa) <- genus_names
  
  # Sample names become a column
  taxa$sample_name <- rownames(taxa)
  rownames(taxa) <- NULL
  
  # Put sample_name first
  taxa <- taxa %>%
    relocate(sample_name)
  
  # Merge metadata
  data_full <- left_join(
    as.data.frame(meta),
    as.data.frame(taxa),
    by = "sample_name"
  )
  
  # Genus columns
  genus_cols <- setdiff(
    names(taxa),
    "sample_name"
  )
  
  #-----------------------------
  # BRAY CURTIS
  #-----------------------------
  
  #===========================================================
  # CREATE PAIR ID
  #===========================================================
  
  pair_data <- data_full %>%
    mutate(
      PairID = paste(
        pmin(tank, connected_tank),
        pmax(tank, connected_tank),
        sep = "_"
      )
    )
  
  
  #===========================================================
  # STORAGE
  #===========================================================
  
  pair_results <- list()
  result_counter <- 1
  
  processed_samples <- character(0)
  
  
  #===========================================================
  # PROCESS EACH TIME POINT
  #===========================================================
  
  for(current_time in unique(pair_data$time)) {
    
    cat("\nProcessing:", current_time, "\n")
    
    time_data <- pair_data %>%
      filter(time == current_time)
    
    
    #=========================================================
    # SEPARATE 0% CONNECTIVITY PONDS
    #=========================================================
    
    zero_conn <- time_data %>%
      filter(flow_rate == 0)
    
    zero_nutrient <- zero_conn %>%
      filter(nutrient_addition == 0)
    
    nutrient_50 <- zero_conn %>%
      filter(nutrient_addition == 50)
    
    nutrient_200 <- zero_conn %>%
      filter(nutrient_addition == 200)
    
    
    #=========================================================
    # FUNCTION TO CALCULATE BRAY-CURTIS BETWEEN TWO SAMPLES
    #=========================================================
    
    calculate_bray <- function(sample1, sample2) {
      
      comm <- rbind(
        as.numeric(data_full[data_full$sample_name == sample1, genus_cols]),
        as.numeric(data_full[data_full$sample_name == sample2, genus_cols])
      )
      
      # Check for NA
      if(any(is.na(comm)))
        return(NA_real_)
      
      # Check for zero-abundance samples
      if(any(rowSums(comm) == 0))
        return(NA_real_)
      
      # Calculate Bray-Curtis
      as.numeric(
        vegdist(comm, method = "bray")
      )
    }
    
    
    #=========================================================
    # 1. PAIR 0N PONDS WITH OTHER 0N PONDS
    #=========================================================
    
    if(nrow(zero_nutrient) >= 2) {
      
      zero_pairs <- combn(
        zero_nutrient$sample_name,
        2,
        simplify = FALSE
      )
      
      for(pair in zero_pairs) {
        
        sample1 <- pair[1]
        sample2 <- pair[2]
        
        bc <- calculate_bray(sample1, sample2)
        
        if(is.na(bc))
          next
        
        tanks <- sort(c(
          zero_nutrient$tank[
            zero_nutrient$sample_name == sample1
          ],
          zero_nutrient$tank[
            zero_nutrient$sample_name == sample2
          ]
        ))
        
        pair_results[[result_counter]] <- data.frame(
          PairID = paste0(
            tanks[1], "_", tanks[2]
          ),
          Time = current_time,
          Day_num = c(
            T2 = 28,
            T3 = 49,
            T4 = 73,
            T5 = 91
          )[current_time],
          Movement = 0,
          NutrientP = 0,
          Bray = bc,
          Sample1 = sample1,
          Sample2 = sample2,
          PairType = "0N_0N"
        )
        
        result_counter <- result_counter + 1
      }
    }
    
    
    #=========================================================
    # 2. RANDOMLY PAIR 50N WITH 0N
    #=========================================================
    
    if(nrow(nutrient_50) > 0 &&
       nrow(zero_nutrient) > 0) {
      
      n_pairs <- min(
        nrow(nutrient_50),
        nrow(zero_nutrient)
      )
      
      # Randomly select unique 0N ponds
      zero_50 <- sample(
        zero_nutrient$sample_name,
        n_pairs,
        replace = FALSE
      )
      
      # Randomize 50N ponds
      nutrient_50_samples <- sample(
        nutrient_50$sample_name,
        n_pairs,
        replace = FALSE
      )
      
      for(k in seq_len(n_pairs)) {
        
        sample1 <- zero_50[k]
        sample2 <- nutrient_50_samples[k]
        
        bc <- calculate_bray(sample1, sample2)
        
        if(is.na(bc))
          next
        
        tanks <- sort(c(
          zero_nutrient$tank[
            zero_nutrient$sample_name == sample1
          ],
          nutrient_50$tank[
            nutrient_50$sample_name == sample2
          ]
        ))
        
        pair_results[[result_counter]] <- data.frame(
          PairID = paste0(
            "ARTIFICIAL_50N_",
            tanks[1], "_", tanks[2]
          ),
          Time = current_time,
          Day_num = c(
            T2 = 28,
            T3 = 49,
            T4 = 73,
            T5 = 91
          )[current_time],
          Movement = 0,
          NutrientP = 50,
          Bray = bc,
          Sample1 = sample1,
          Sample2 = sample2,
          PairType = "0N_50N"
        )
        
        result_counter <- result_counter + 1
      }
    }
    
    
    #=========================================================
    # 3. RANDOMLY PAIR 200N WITH 0N
    #=========================================================
    
    if(nrow(nutrient_200) > 0 &&
       nrow(zero_nutrient) > 0) {
      
      n_pairs <- min(
        nrow(nutrient_200),
        nrow(zero_nutrient)
      )
      
      # Randomly select unique 0N ponds
      zero_200 <- sample(
        zero_nutrient$sample_name,
        n_pairs,
        replace = FALSE
      )
      
      # Randomize 200N ponds
      nutrient_200_samples <- sample(
        nutrient_200$sample_name,
        n_pairs,
        replace = FALSE
      )
      
      for(k in seq_len(n_pairs)) {
        
        sample1 <- zero_200[k]
        sample2 <- nutrient_200_samples[k]
        
        bc <- calculate_bray(sample1, sample2)
        
        if(is.na(bc))
          next
        
        tanks <- sort(c(
          zero_nutrient$tank[
            zero_nutrient$sample_name == sample1
          ],
          nutrient_200$tank[
            nutrient_200$sample_name == sample2
          ]
        ))
        
        pair_results[[result_counter]] <- data.frame(
          PairID = paste0(
            "ARTIFICIAL_200N_",
            tanks[1], "_", tanks[2]
          ),
          Time = current_time,
          Day_num = c(
            T2 = 28,
            T3 = 49,
            T4 = 73,
            T5 = 91
          )[current_time],
          Movement = 0,
          NutrientP = 200,
          Bray = bc,
          Sample1 = sample1,
          Sample2 = sample2,
          PairType = "0N_200N"
        )
        
        result_counter <- result_counter + 1
      }
    }
    
    
    #=========================================================
    # MARK ALL 0% CONNECTIVITY SAMPLES AS PROCESSED
    #=========================================================
    
    processed_samples <- c(
      processed_samples,
      zero_conn$sample_name
    )
    
    
    #=========================================================
    # PROCESS NORMAL CONNECTED PAIRS
    #=========================================================
    
    current_time_data <- time_data %>%
      filter(flow_rate > 0)
    
    for(i in seq_len(nrow(current_time_data))) {
      
      current <- current_time_data[i, ]
      
      sample_name <- current$sample_name
      
      
      #-----------------------------------------
      # Skip if already processed
      #-----------------------------------------
      
      if(sample_name %in% processed_samples)
        next
      
      
      #-----------------------------------------
      # Find connected pond
      #-----------------------------------------
      
      partner <- time_data %>%
        filter(
          tank == current$connected_tank,
          time == current_time,
          flow_rate > 0
        )
      
      
      #-----------------------------------------
      # No partner found
      #-----------------------------------------
      
      if(nrow(partner) == 0) {
        
        processed_samples <- c(
          processed_samples,
          sample_name
        )
        
        next
      }
      
      
      partner_name <- partner$sample_name[1]
      
      
      #-----------------------------------------
      # Partner already processed
      #-----------------------------------------
      
      if(partner_name %in% processed_samples)
        next
      
      
      #-----------------------------------------
      # Calculate Bray-Curtis
      #-----------------------------------------
      
      bc <- calculate_bray(
        sample_name,
        partner_name
      )
      
      if(is.na(bc)) {
        
        processed_samples <- c(
          processed_samples,
          sample_name,
          partner_name
        )
        
        next
      }
      
      
      #-----------------------------------------
      # Pair ID
      #-----------------------------------------
      
      tanks <- sort(c(
        current$tank,
        partner$tank[1]
      ))
      
      new_pair_id <- paste(
        tanks[1],
        tanks[2],
        sep = "_"
      )
      
      
      #-----------------------------------------
      # Store result
      #-----------------------------------------
      
      pair_results[[result_counter]] <- data.frame(
        PairID = new_pair_id,
        Time = current_time,
        Day_num = c(
          T2 = 28,
          T3 = 49,
          T4 = 73,
          T5 = 91
        )[current_time],
        Movement = current$flow_rate,
        NutrientP = current$nutrient_addition,
        Bray = bc,
        Sample1 = sample_name,
        Sample2 = partner_name,
        PairType = "CONNECTED"
      )
      
      result_counter <- result_counter + 1
      
      
      #-----------------------------------------
      # Mark both as processed
      #-----------------------------------------
      
      processed_samples <- c(
        processed_samples,
        sample_name,
        partner_name
      )
    }
  }
  
  
  #===========================================================
  # COMBINE ALL RESULTS
  #===========================================================
  
  div <- bind_rows(pair_results)
  
  
  #===========================================================
  # SAVE TABLE
  #===========================================================
  
  fwrite(
    div,
    file.path(
      output_dir,
      paste0(taxa_name, "_paired_bray.tsv")
    ),
    sep = "\t"
  )
  
  #-----------------------------
  # BRAY CURTIS OVERTIME
  #-----------------------------
  
  nutrient_cols <- c(
    "0"   = "#440154FF",
    "25"  = "#39568CFF",
    "50"  = "#20A387FF",
    "100" = "#73D055FF",
    "200" = "#FDE725FF"
  )
  
  bray_mean <- div %>%
    filter(Movement != 0) %>%
    group_by(Day_num, Movement, NutrientP) %>%
    summarise(
      Bray = mean(Bray, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      NutrientP = factor(
        NutrientP,
        levels = c(0, 25, 50, 100, 200)
      )
    )
  
  plot_bray <- function(output_file){
    
    p <- ggplot(
      bray_mean,
      aes(
        x = Day_num,
        y = Bray,
        colour = NutrientP,
        group = NutrientP
      )
    ) +
      
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      
      facet_wrap(
        ~ Movement,
        nrow = 1,
        scales = "fixed"
      ) +
      
      scale_color_manual(
        values = nutrient_cols,
        name = "Nutrient Addition"
      ) +
      
      labs(
        x = "Day",
        y = "Bray-Curtis Dissimilarity"
      ) +
      
      theme_classic() +
      
      theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom"
      )
    
    ggsave(
      output_file,
      p,
      width = 12,
      height = 4,
      dpi = 600
    )
    
    p
  }
  
  output_plot <- file.path(
    output_dir,
    paste0(taxa_name, "_bray_time.png")
  )
  
  p_bray <- plot_bray(output_plot)
  
  #-----------------------------
  # PCoA OVER TIME
  #-----------------------------
  
  output_pcoa <- file.path(
    output_dir,
    paste0(taxa_name, "_pcoa_time.png")
  )
  
  run_pcoa_time <- function(timepoint) {
    
    df <- data_full %>%
      filter(time == timepoint)
    
    if(nrow(df) < 3) return(NULL)
    
    comm <- as.matrix(df[, genus_cols, drop = FALSE])
    
    comm[is.na(comm)] <- 0
    
    bray <- vegdist(
      comm,
      method = "bray"
    )
    
    ord <- cmdscale(
      bray,
      eig = TRUE,
      k = 2
    )
    
    scores <- as.data.frame(ord$points)
    colnames(scores) <- c("PC1","PC2")
    
    df_plot <- cbind(df, scores)
    
    #------------------------------------------
    # REMOVE BLANKS FROM PLOT
    #------------------------------------------
    
    df_plot <- df_plot %>%
      filter(SampleType != "blank")
    
    #------------------------------------------
    # Nutrients
    #------------------------------------------
    
    df_plot$nutrient_addition[
      is.na(df_plot$nutrient_addition)
    ] <- "control"
    
    df_plot$nutrient_addition <- factor(
      as.character(df_plot$nutrient_addition),
      levels = c(
        "control",
        "0",
        "25",
        "50",
        "100",
        "200"
      )
    )
    
    #------------------------------------------
    # Shapes
    #------------------------------------------
    
    df_plot$shape_group <- case_when(
      df_plot$SampleType == "resv" ~ "Reservoir",
      df_plot$upstream_downstream == "Up" ~ "Upstream",
      df_plot$upstream_downstream == "Down" ~ "Downstream",
      TRUE ~ NA_character_
    )
    
    df_plot$shape_group <- factor(
      df_plot$shape_group,
      levels = c(
        "Upstream",
        "Downstream",
        "Reservoir"
      )
    )
    
    df_plot <- df_plot %>%
      filter(!is.na(shape_group))
    
    #------------------------------------------
    # Connectivity dataframe
    #------------------------------------------
    
    connections <- data.frame(
      x = numeric(),
      y = numeric(),
      xend = numeric(),
      yend = numeric(),
      flow = numeric()
    )
    
    for(i in seq_len(nrow(df_plot))) {
      
      row_i <- df_plot[i, ]
      
      if(
        is.na(row_i$flow_rate) ||
        row_i$flow_rate == 0
      ) next
      
      target <- df_plot %>%
        filter(
          tank == row_i$connected_tank
        )
      
      if(nrow(target) == 0) next
      
      connections <- rbind(
        connections,
        data.frame(
          x = row_i$PC1,
          y = row_i$PC2,
          xend = target$PC1[1],
          yend = target$PC2[1],
          flow = row_i$flow_rate
        )
      )
    }
    
    #------------------------------------------
    # Plot
    #------------------------------------------
    
    p <- ggplot(
      df_plot,
      aes(PC1, PC2)
    ) +
      
      geom_segment(
        data = connections,
        aes(
          x = x,
          y = y,
          xend = xend,
          yend = yend,
          size = flow
        ),
        arrow = arrow(
          length = unit(0.15, "cm")
        ),
        inherit.aes = FALSE
      ) +
      
      geom_point(
        aes(
          color = nutrient_addition,
          shape = shape_group
        ),
        size = 3,
        stroke = 1
      ) +
      
      scale_color_manual(
        values = c(
          "control" = "#000000",
          "0" = "#440154FF",
          "25" = "#39568CFF",
          "50" = "#20A387FF",
          "100" = "#73D055FF",
          "200" = "#FDE725FF"
        ),
        name = "Nutrient Enrichment"
      ) +
      
      scale_shape_manual(
        values = c(
          "Upstream" = 16,
          "Downstream" = 18,
          "Reservoir" = 1
        ),
        limits = c(
          "Upstream",
          "Downstream",
          "Reservoir"
        ),
        drop = FALSE,
        name = "Node Position"
      ) +
      
      scale_size(
        range = c(0.5, 2),
        name = "Connectivity (% week⁻¹)"
      ) +
      
      labs(
        x = NULL,
        y = NULL
      ) +
      
      ggtitle(
        paste(
          "Time",
          gsub("T", "", timepoint)
        )
      ) +
      
      theme_classic() +
      
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none"
      )
    
    return(p)
  }
  
  #==================================================
  # RUN PCoA PLOTS
  #==================================================
  
  pcoa_times <- c("T2", "T3", "T4", "T5")
  
  pcoa_plots <- lapply(
    pcoa_times,
    run_pcoa_time
  )
  
  #==================================================
  # GLOBAL AXIS LIMITS
  #==================================================
  
  all_scores <- data.frame()
  
  for(t in pcoa_times) {
    
    df <- data_full %>%
      filter(time == t)
    
    if(nrow(df) < 3) next
    
    comm <- as.matrix(
      df[, genus_cols, drop = FALSE]
    )
    
    comm[is.na(comm)] <- 0
    
    bray_t <- vegdist(
      comm,
      method = "bray"
    )
    
    ord <- cmdscale(
      bray_t,
      eig = TRUE,
      k = 2
    )
    
    scores <- as.data.frame(
      ord$points
    )
    
    all_scores <- rbind(
      all_scores,
      scores
    )
  }
  
  pcoa_xlim <- range(all_scores$V1)
  pcoa_ylim <- range(all_scores$V2)
  
  pcoa_plots <- lapply(
    pcoa_plots,
    function(p) {
      p +
        coord_cartesian(
          xlim = pcoa_xlim,
          ylim = pcoa_ylim
        )
    }
  )
  
  #==================================================
  # CREATE THE UNIVERSAL LEGEND
  #==================================================
  
  legend_points <- expand.grid(
    nutrient_addition = c("control", "0", "25", "50", "100", "200"),
    shape_group = c("Upstream", "Downstream", "Reservoir")
  )
  
  legend_points$PC1 <- seq_len(nrow(legend_points))
  legend_points$PC2 <- 1
  
  legend_lines <- data.frame(
    x = 1,
    y = 1,
    xend = 2,
    yend = 1,
    flow = c(0, 25, 50, 100, 200)
  )
  
  legend_plot <- ggplot() +
    
    geom_segment(
      data = legend_lines,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        size = flow
      ),
      arrow = arrow(length = unit(0.15, "cm"))
    ) +
    
    geom_point(
      data = legend_points,
      aes(
        PC1,
        PC2,
        colour = nutrient_addition,
        shape = shape_group
      )
    ) +
    
    scale_color_manual(
      values = c(
        "control" = "#000000",
        "0" = "#440154FF",
        "25" = "#39568CFF",
        "50" = "#20A387FF",
        "100" = "#73D055FF",
        "200" = "#FDE725FF"
      ),
      limits = c("control","0","25","50","100","200"),
      drop = FALSE,
      name = "Nutrient Enrichment"
    ) +
    
    scale_shape_manual(
      values = c(
        "Upstream" = 16,
        "Downstream" = 18,
        "Reservoir" = 1
      ),
      limits = c(
        "Upstream",
        "Downstream",
        "Reservoir"
      ),
      drop = FALSE,
      name = "Pond Position"
    ) +
    
    scale_size_continuous(
      range = c(0.5, 2),
      name = "Connectivity"
    ) +
    
    theme_void() +
    theme(
      legend.position = "bottom"
    )
  
  pcoa_legend <- get_legend(legend_plot)
  
  #==================================================
  # COMBINE PCoA PLOTS
  #==================================================
  
  pcoa_final_plot <- plot_grid(
    wrap_plots(pcoa_plots, nrow = 1),
    pcoa_legend,
    ncol = 1,
    rel_heights = c(1, 0.15)
  )
  
  ggsave(
    filename = output_pcoa,
    plot = pcoa_final_plot,
    width = 14,
    height = 4,
    units = "in",
    dpi = 600
  )
  
  #-----------------------------
  # GAMM
  #-----------------------------
  
  # Clean beta diversity for GAMM
  bray <- div %>%
    mutate(
      
      NutrientP = as.numeric(NutrientP),
      Movement = as.numeric(Movement),
      Day_num = as.numeric(Day_num),
      
      # Random effect
      PairID = as.factor(PairID)
      
    )
  
  # Smithson and Verkuilen adjustment
  n <- nrow(bray)
  
  bray <- bray %>%
    mutate(
      Bray_beta = (Bray * (n - 1) + 0.5) / n
    )
  
  #GAMM
  k_day <- 4
  k <- 4
  
  
  formula_text <- glue(
    "Bray_beta ~

ti(Day_num, k = {k_day}) +

ti(Movement, k = {k}) +

ti(NutrientP, k = {k}) +

ti(Movement, Day_num, 
   k = c({k}, {k_day})) +

ti(NutrientP, Day_num,
   k = c({k}, {k_day})) +

ti(Movement, NutrientP,
   k = c({k}, {k})) +

ti(Day_num, Movement, NutrientP,
   k = c({k_day}, {k}, {k})) +

s(PairID, bs='re')"
  )
  
  
  bray_model <- gam(
    as.formula(formula_text),
    data = bray,
    family = betar(link="logit"),
    method = "REML"
  )
  
  summary(bray_model)
  
  # -----------------------------
  # GAM CHECK (diagnostics)
  # -----------------------------
  
  # Save the 4-panel diagnostic plots
  png(
    file.path(output_dir, paste0(taxa_name, "_beta_gamcheck.png")),
    width = 8, height = 8, units = "in", res = 300
  )
  gam.check(bray_model)
  dev.off()
  
  # Save the printed convergence + k-index text
  check_text <- capture.output(gam.check(bray_model))
  writeLines(
    check_text,
    file.path(output_dir, paste0(taxa_name, "_beta_gamcheck.txt"))
  )
  
  # Save the k-index table itself as a clean data frame
  k_table <- as.data.frame(k.check(bray_model))
  k_table$Term <- rownames(k_table)
  rownames(k_table) <- NULL
  
  fwrite(
    k_table,
    file.path(output_dir, paste0(taxa_name, "_beta_kcheck.tsv")),
    sep = "\t"
  )
  
  p1 <- summary(bray_model)$dev.expl
  
  p2 <- deviance(bray_model)
  
  p3 <- bray_model$deviance
  
  p4 <- bray_model$null.deviance
  
  p5 <-  range(bray$Bray, na.rm = TRUE)
  
  p6 <- any(bray$Bray <= 0 | bray$Bray >= 1, na.rm = TRUE)
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  
  # Save model
  saveRDS(
    bray_model,
    file.path(
      output_dir,
      paste0(taxa_name, "_bray_GAMM.rds")
    )
  )
  
  #-----------------------------
  # SAVE GAMM TABLES
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
    "ti(Day_num,Movement,NutrientP)" = "ti(Connectivity, Nutrient Enrichment, Day)",
    "s(PairID)" = "s(PairID)"
  )
  
  # Save GAMM table
  gam_table <- as.data.frame(
    summary(bray_model)$s.table
  )
  
  # Move rownames into column
  gam_table$Term <- rownames(gam_table)
  rownames(gam_table) <- NULL
  
  # Test
  print(colnames(gam_table))
  
  # Rename columns
  gam_table <- gam_table %>%
    select(
      Term,
      edf = edf,
      χ2 = `Chi.sq`,
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
      χ2 = round(χ2, 2),
      p_value = sapply(p_value, format_p)
    ) %>%
    mutate(
      Term = as.character(Term)
    )
  
  # Save GAMM table
  fwrite(
    gam_table,
    file.path(
      output_dir,
      paste0(taxa_name,"_bray_GAM_table.tsv")
    ),
    sep="\t"
  )
  
  # Save Deviance Explained
  deviance_explained <- summary(bray_model)$dev.expl * 100
  
  deviance_table <- data.frame(
    Deviance_explained = round(deviance_explained, 2)
  )
  
  fwrite(
    deviance_table,
    file.path(
      output_dir,
      paste0(taxa_name, "_bray_deviance_explained.tsv")
    ),
    sep = "\t"
  )
  
  #-----------------------------
  # PREDICTION GRAPHS
  #-----------------------------
  
  time_points <- c(28, 49, 73, 91)
  
  # Prediction grid
  pred_grid <- expand.grid(
    
    Movement = seq(
      min(bray$Movement),
      max(bray$Movement),
      length.out = 150
    ),
    
    NutrientP = seq(
      min(bray$NutrientP),
      max(bray$NutrientP),
      length.out = 150
    ),
    
    Day_num = time_points,
    
    # random effect placeholder
    PairID = levels(bray$PairID)[1]
  )
  
  pred_grid$Time <- factor(
    pred_grid$Day_num,
    levels = time_points
  )

  # Predictions
  pred_grid$Prediction <- predict(
    bray_model,
    newdata = pred_grid,
    type = "response",
    exclude = "s(PairID)"
  )
  
  # Contour bins
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
  
  # Plot
  p_bray_surface <- ggplot(
    pred_grid,
    aes(
      x = NutrientP,
      y = Movement
    )
  ) +
    
    geom_raster(
      aes(fill = PredictionBin)
    ) +
    
    geom_contour(
      aes(z = Prediction),
      breaks = contour_breaks,
      colour = "black",
      linewidth = 0.4
    ) +
    
    facet_wrap(
      ~Time,
      nrow = 1
    ) +
    
    scale_fill_viridis_d(
      option = "viridis",
      name = "Predicted\nBray-Curtis\nDissimilarity"
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
  
  # Save plot
  ggsave(
    filename = file.path(
      output_dir,
      paste0(taxa_name, "_bray_surface.png")
    ),
    plot = p_bray_surface,
    width = 12,
    height = 5,
    dpi = 600
  )
}
  

#==================================================
# RUN PIPELINE FUNCTION
#==================================================

result_23S <- run_beta_diversity_pipeline("23S_total")
result_23S_green_algae <- run_beta_diversity_pipeline("23S_green_algae")
result_23S_cyanobacteria <- run_beta_diversity_pipeline("23S_cyanobacteria")

result_18S <- run_beta_diversity_pipeline("18S_total")
result_18S_phytoplankton <- run_beta_diversity_pipeline("18S_phytoplankton")
result_18S_green_algae <- run_beta_diversity_pipeline("18S_green_algae")
result_18S_zooplankton <- run_beta_diversity_pipeline("18S_zooplankton")
result_18S_micro_heterotrophs <- run_beta_diversity_pipeline("18S_micro_heterotrophs")
result_18S_fungi <- run_beta_diversity_pipeline("18S_fungi")

result_COI <- run_beta_diversity_pipeline("COI_total")
result_COI_zooplankton <- run_beta_diversity_pipeline("COI_zooplankton")
  