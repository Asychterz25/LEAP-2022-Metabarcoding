library(data.table)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(patchwork)
library(mgcv)
library(lme4)
library(cowplot)

#==================================================
# 1) LOAD FILES
#==================================================

#23S
taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/23S/23S_combined_taxa_aggregated.tsv"
metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/23S_LEAP22_metadata.tsv"
output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_total/23S_BD_plot.png"
output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_total/23S_centroid_plot.png"
output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_total/23S_BD_permanova.csv"
centroid_name = "23S_total_centorids"


  #==================================================
  # 1) LOAD DATA
  #==================================================
  
  taxa <- fread(taxa_file)
  meta <- fread(metadata_file)
  
  #==================================================
  # 2) CLEAN TO GENUS LEVEL
  #==================================================
  
  sample_cols <- setdiff(
    names(taxa),
    c(
      "domain","supergroup","division","subdivision",
      "class","order","family","species"
    )
  )
  
  taxa_genus <- taxa
  
  taxa_genus <- taxa_genus[
    !genus %in% c(
      "Unknown",
      "unknown",
      "Unassigned",
      "unassigned"
    )
  ]
  
  taxa_genus <- taxa_genus %>%
    group_by(genus) %>%
    summarise(
      across(
        where(is.numeric),
        sum,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  sample_matrix <- taxa_genus[, -1, with = FALSE]
  col_totals <- colSums(sample_matrix)
  
  empty_samples <- names(col_totals[col_totals == 0])
  
  if(length(empty_samples) > 0) {
    
    message(
      "Removing empty samples: ",
      paste(empty_samples, collapse = ", ")
    )
    
    taxa_genus <- taxa_genus[
      ,
      !names(taxa_genus) %in% empty_samples,
      with = FALSE
    ]
  }
  
  taxa_genus_rel <- taxa_genus
  
  taxa_genus_rel[, -1] <- sweep(
    taxa_genus[, -1],
    2,
    colSums(taxa_genus[, -1]),
    "/"
  ) * 100
  
  #==================================================
  # TRANSPOSE TO SAMPLE x GENUS AND COLLAPSE REPLICATES
  #==================================================
  
  taxa_mat <- as.data.frame(taxa_genus_rel)
  
  rownames(taxa_mat) <- taxa_mat$genus
  taxa_mat$genus <- NULL
  
  taxa_mat <- t(taxa_mat)
  taxa_mat <- as.data.frame(taxa_mat)
  
  taxa_mat$sample_name <- rownames(taxa_mat)
  
  #==================================================
  # COLLAPSE TECHNICAL REPLICATES
  #==================================================
  
  taxa_mat$sample_base <- sub("_R[0-9]+$", "", taxa_mat$sample_name)
  
  taxa_mat_collapsed <- taxa_mat %>%
    group_by(sample_base) %>%
    summarise(
      across(
        where(is.numeric),
        mean,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  names(taxa_mat_collapsed)[1] <- "sample_name"
  
  # Metadata collapse
  
  meta$sample_base <- sub("_R[0-9]+$", "", meta$sample_name)
  
  meta_collapsed <- meta %>%
    group_by(sample_base) %>%
    slice(1) %>%
    ungroup()
  
  meta_collapsed$sample_name <- meta_collapsed$sample_base
  
  meta_collapsed$sample_base <- NULL
  
  # Merge collapsed data
  
  data_full <- merge(
    meta_collapsed,
    taxa_mat_collapsed,
    by = "sample_name"
  )
  
  genus_cols <- setdiff(
    colnames(taxa_mat_collapsed),
    "sample_name"
  )
  
  #==================================================
  # 3) BRAY CURTIS AND PCOA FUNCTION
  #==================================================
  
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
          tank == row_i$connected_tank,
          replicate == row_i$replicate
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
          "0" = "#0072B2",
          "25" = "#56B4E9",
          "50" = "#E69F00",
          "100" = "#F0E442",
          "200" = "#009E73"
        ),
        name = "Nutrient Addition"
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
      
      scale_size(
        range = c(0.5, 2),
        name = "Flow Rate"
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
    
      df_plot <- df_plot %>%
        filter(!is.na(shape_group))
    
    return(p)
  }
  
  #==================================================
  # RUN PLOTS
  #==================================================
  
  times <- c(
    "T1",
    "T2",
    "T3",
    "T4",
    "T5"
  )
  
  plots <- lapply(
    times,
    run_pcoa_time
  )
  
  #==================================================
  # GLOBAL AXIS LIMITS
  #==================================================
  
  all_scores <- data.frame()
  
  for(t in times) {
    
    df <- data_full %>%
      filter(time == t)
    
    if(nrow(df) < 3) next
    
    comm <- as.matrix(
      df[, genus_cols, drop = FALSE]
    )
    
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
    
    scores <- as.data.frame(
      ord$points
    )
    
    all_scores <- rbind(
      all_scores,
      scores
    )
  }
  
  xlim <- range(all_scores$V1)
  ylim <- range(all_scores$V2)
  
  plots <- lapply(
    plots,
    function(p) {
      p +
        coord_cartesian(
          xlim = xlim,
          ylim = ylim
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
  
  #Create the legend
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
        "0" = "#0072B2",
        "25" = "#56B4E9",
        "50" = "#E69F00",
        "100" = "#F0E442",
        "200" = "#009E73"
      ),
      limits = c("control","0","25","50","100","200"),
      drop = FALSE,
      name = "Nutrient Addition"
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
      name = "Flow Rate"
    ) +
    
    theme_void() +
    theme(
      legend.position = "bottom"
    )
  
  legend <- get_legend(legend_plot)
  
  #==================================================
  # COMBINE PLOTS
  #==================================================
  
  final_plot <- plot_grid(
    wrap_plots(plots, nrow = 1),
    legend,
    ncol = 1,
    rel_heights = c(1, 0.15)
  )
  
  print(final_plot)
  
  ggsave(
    filename = output_plot,
    plot = final_plot,
    width = 14,
    height = 4,
    units = "in",
    dpi = 600
  )
  
  #==================================================
  # 4) MULTIVARIATE STATISTICS
  #==================================================
  
  stats_file <- output_permanova
  
  df_perm <- data_full %>%
    filter(
      SampleType != "blank",
      SampleType != "resv"
    ) %>%
    filter(
      !is.na(upstream_downstream),
      !is.na(time)
    )
  
  #--------------------------------------------------
  # Clean treatment variables
  #--------------------------------------------------
  
  #Keep only samples
  df_perm <- df_perm[
    df_perm$SampleType == "sample",
  ]
  
  df_perm$upstream_downstream <- factor(
    df_perm$upstream_downstream
  )
  
  df_perm$time <- factor(
    df_perm$time
  )
  
  df_perm$nutrient_addition[
    is.na(df_perm$nutrient_addition)
  ] <- "control"
  
  df_perm$nutrient_addition <- factor(
    as.character(df_perm$nutrient_addition)
  )
  
  # Remove 0 Flow Rate
  df_perm <- df_perm[df_perm$flow_rate != 0, ]
  
  #Remove any NA
  df_perm <- df_perm[!is.na(df_perm$sample_name), ]
  
  #--------------------------------------------------
  # Bray matrix
  #--------------------------------------------------
  
  comm <- as.matrix(df_perm[, genus_cols, drop = FALSE])
  rownames(comm) <- df_perm$sample
  
  comm[is.na(comm)] <- 0
  
  bray <- vegdist(comm, method = "bray")
  
  #==================================================
  # PERMANOVA
  #==================================================
  
  perm <- tryCatch({
    
    adonis2(
      bray ~
        upstream_downstream +
        nutrient_addition +
        flow_rate +
        time,
      data = df_perm,
      permutations = 999
    )
    
  }, error = function(e){
    
    e
    
  })
  
  #==================================================
  # BETADISPER
  #==================================================
  
  disp_anova <- NULL
  disp_perm <- NULL
  
  #-----------------------------
  # 1. CLEAN METADATA
  #-----------------------------
  df_perm <- df_perm[
    !is.na(df_perm$sample) &
      df_perm$sample != "",
  ]
  
  df_perm$sample <- as.character(df_perm$sample)
  
  #-----------------------------
  # 2. BUILD COMMON IDS
  #-----------------------------
  bray_ids <- rownames(as.matrix(bray))
  meta_ids <- df_perm$sample
  
  common <- intersect(bray_ids, meta_ids)
  
  #-----------------------------
  # 3. HARD SUBSET BRAY (ORDERED)
  #-----------------------------
  bray_mat <- as.matrix(bray)
  bray_mat <- bray_mat[common, common, drop = FALSE]
  bray <- as.dist(bray_mat)
  
  #-----------------------------
  # 4. HARD SUBSET + ORDER METADATA
  #-----------------------------
  df_perm <- df_perm[match(common, df_perm$sample), ]
  
  # safety check
  stopifnot(all(df_perm$sample == common))
  
  #-----------------------------
  # 5. ENSURE VALID GROUP
  #-----------------------------
  df_perm$upstream_downstream <- droplevels(
    factor(df_perm$upstream_downstream)
  )
  
  #-----------------------------
  # 6. RUN BETADISPER SAFELY
  #-----------------------------
  if (nlevels(df_perm$upstream_downstream) > 1 &&
      nrow(df_perm) == attr(bray, "Size")) {
    
    disp <- betadisper(
      bray,
      df_perm$upstream_downstream
    )
    
    disp_anova <- anova(disp)
    
    disp_perm <- permutest(
      disp,
      permutations = 999
    )
  }
  #==================================================
  # CENTROID DISTANCES
  #==================================================
  
  centroid_distances <- list()
  
  for(tp in unique(df_perm$time)){
    
    df_tp <- df_perm %>%
      filter(time == tp)
    
    if(nrow(df_tp) < 4) next
    
    comm_tp <- as.matrix(
      df_tp[, genus_cols, drop = FALSE]
    )
    
    comm_tp[is.na(comm_tp)] <- 0
    
    bray_tp <- vegdist(
      comm_tp,
      method = "bray"
    )
    
    ord <- cmdscale(
      bray_tp,
      eig = TRUE,
      k = 2
    )
    
    scores <- as.data.frame(
      ord$points
    )
    
    colnames(scores) <- c(
      "PC1",
      "PC2"
    )
    
    df_tp <- cbind(
      df_tp,
      scores
    )
    
    centroids <- df_tp %>%
      group_by(
        nutrient_addition,
        flow_rate,
        upstream_downstream
      ) %>%
      summarise(
        PC1 = mean(PC1),
        PC2 = mean(PC2),
        .groups = "drop"
      )
    
    treatments <- unique(
      centroids[, c(
        "nutrient_addition",
        "flow_rate"
      )]
    )
    
    for(i in seq_len(nrow(treatments))){
      
      nut <- treatments$nutrient_addition[i]
      flow <- treatments$flow_rate[i]
      
      pair <- centroids %>%
        filter(
          nutrient_addition == nut,
          flow_rate == flow
        )
      
      if(nrow(pair) != 2) next
      
      d <- sqrt(
        (pair$PC1[1] - pair$PC1[2])^2 +
          (pair$PC2[1] - pair$PC2[2])^2
      )
      
      centroid_distances[[length(centroid_distances)+1]] <-
        data.frame(
          time = as.character(tp),
          nutrient_addition = as.character(nut),
          flow_rate = flow,
          distance = d
        )
      
    }
    
  }
  
  centroid_distances <- bind_rows(
    centroid_distances
  )
  
  #==================================================
  # CENTROID MODEL
  #==================================================
  
  centroid_model <- NULL
  
  if(nrow(centroid_distances) > 5){
    
    centroid_distances$time <- factor(
      centroid_distances$time
    )
    
    centroid_model <- lm(
      distance ~
        nutrient_addition +
        flow_rate +
        time,
      data = centroid_distances
    )
    
  }
  
  #==================================================
  # SAVE STATS TO ONE FILE
  #==================================================
  
  write_section <- function(title, content, con) {
    writeLines("\n========================================", con)
    writeLines(title, con)
    writeLines("----------------------------------------", con)
    writeLines("", con)
    
    if (!is.null(content)) {
      capture.output(print(content), file = con)
    } else {
      writeLines("NULL", con)
    }
  }
  
  con <- file(stats_file, open = "wt")
  
  # Header
  writeLines("========================================", con)
  writeLines(paste("Analysis:", centroid_name), con)
  writeLines("========================================", con)
  
  # PERMANOVA
  write_section("GLOBAL PERMANOVA", perm, con)
  
  # BETADISPER (combined cleanly)
  write_section("BETADISPER ANOVA", disp_anova, con)
  write_section("BETADISPER PERMUTEST", disp_perm, con)
  
  # CENTROID MODEL
  write_section("CENTROID DISTANCE MODEL", centroid_model, con)
  
  # CENTROID RAW DATA
  write_section("CENTROID DISTANCES", centroid_distances, con)
  
  close(con)
  
  #==================================================
  # 5) PLOT CENTROIDS OVER TIME
  #==================================================
  
  #--------------------------------------------------
  # 1. Recode time to numeric days
  #--------------------------------------------------
  
  centroid_distances <- centroid_distances %>%
    mutate(
      day = case_when(
        time == "T1" ~ 0,
        time == "T2" ~ 28,
        time == "T3" ~ 49,
        time == "T4" ~ 73,
        time == "T5" ~ 91,
        TRUE ~ NA_real_
      )
    )
  
  #--------------------------------------------------
  # 2. Ensure correct factor structure + colours
  #--------------------------------------------------
  
  centroid_distances$nutrient_addition <- as.factor(centroid_distances$nutrient_addition)
  centroid_distances$flow_rate <- as.factor(centroid_distances$flow_rate)
  
  nutrient_cols <- c(
    "0" = "#0072B2",
    "25"  = "#56B4E9",
    "50"  = "#E69F00",
    "100" = "#F0E442",
    "200" = "#009E73"
  )
  
  #--------------------------------------------------
  # 3. Plot
  #--------------------------------------------------
  
  p <- ggplot(centroid_distances,
              aes(x = day,
                  y = distance,
                  color = nutrient_addition,
                  group = nutrient_addition)) +
    
    geom_line(na.rm = TRUE) +
    geom_point(size = 2, na.rm = TRUE) +
    
    facet_wrap(~ flow_rate, nrow = 1, scales = "fixed") +
    
    scale_color_manual(
      values = nutrient_cols,
      breaks = c("0", "25", "50", "100", "200"),
      limits = c("0", "25", "50", "100", "200"),
      na.value = NA
    ) +
    
    labs(
      x = "Day",
      y = "Centroid distances between upstream and downstream mesocosms",
      color = "Nutrient Addition"
    ) +
    
    theme_classic() +
    
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  ggsave(
    filename = output_centroid_plot,
    plot = p,
    width = 14,
    height = 4,
    units = "in",
    dpi = 600
  )
  
  
  
  
  
  
  
  
  
  
  #==================================================
  #Sample Paths
  #==================================================
  
  
  #23S
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/23S/23S_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/23S_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_total/23S_BD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_total/23S_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_total/23S_BD_permanova.csv"
  centroid_name = "23S_total_centorids"

#23S Green Algae
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/23S_green_algae/23S_green_algae_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/23S_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_green_algae/23S_green_algaeBD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_green_algae/23S_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_green_algae/23S_green_algaeBD_permanova.csv"
  centroid_name = "23S_green_algae_centorids"
  
  #23S Cyanobacteria
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/23S_cyanobacteria/23S_cyanobacteria_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/23S_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_cyanobacteria/23S_cyanobacteriaBD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_cyanobacteria/23S_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/23S_cyanobacteria/23S_cyanobacteriaBD_permanova.csv"
  centroid_name = "23S_cyanobacteria_centorids"

  #COI Total
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/COI/COI_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/COI_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/COI_total/COI_BD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/COI_total/COI_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/COI_total/COI_BD_permanova.csv"
  centroid_name = "COI_centorids"

  #COI Zooplankton
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/COI_zooplankton/COI_zooplankton_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/COI_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/COI_zooplankton/COI_BD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/COI_zooplankton/COI_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/COI_zooplankton/COI_BD_permanova.csv"
  centroid_name = "COI_zooplankton_centorids"

  #18S Total
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/18S/18S_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/18S_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_total/18S_BD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_total/18S_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_total/18S_BD_permanova.csv"
  centroid_name = "18S_total_centorids"

  #18S Phytoplankton
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/18S_phytoplankton/18S_phytoplankton_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/18S_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_phytoplankton/18S_BD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_phytoplankton/18S_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_phytoplankton/18S_BD_permanova.csv"
  centroid_name = "18S_phytoplankton_centorids"

  #18S micro_heterotophs
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/18S_micro_heterotophs/18S_micro_heterotrophs_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/18S_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_total/18S_BD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_micro_heterotophs/18S_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_micro_heterotophs/18S_BD_permanova.csv"
  centroid_name = "18S_micro_heterotophs_centorids"

  #18S zooplankton
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/18S_zooplankton/18S_zooplankton_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/18S_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_zooplankton/18S_BD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_zooplankton/18S_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_zooplankton/18S_BD_permanova.csv"
  centroid_name = "18S_zooplankton_centorids"

  #18S fungi
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/18S_fungi/18S_fungi_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/18S_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_fungi/18S_BD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_fungi/18S_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_fungi/18S_BD_permanova.csv"
  centroid_name = "18S_fungi_centorids"

  #18S green_algae
  taxa_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_lists/18S_green_algae/18S_green_algae_combined_taxa_aggregated.tsv"
  metadata_file = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/18S_LEAP22_metadata.tsv"
  output_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_green_algae/18S_BD_plot.png"
  output_centroid_plot = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_green_algae/18S_centroid_plot.png"
  output_permanova = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/beta_diversity/18S_green_algae/18S_BD_permanova.csv"
  centroid_name = "18S_green_algae_centorids"

