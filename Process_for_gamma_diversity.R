############## Combine Upstream and Downstream Ponds for Gamma Diversity ###################

#==================================================
# MAIN PIPELINE FUNCTION
#==================================================

run_gamma_diversity_pipeline <- function(taxa_name, metadata,
                                         base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2",
                                         save_figs = TRUE) {
  
  #-----------------------------
  # FILE PATHS
  #-----------------------------
  
  taxa_path <- file.path(base_path, "taxa_lists", "cleaned",
                         paste0(taxa_name, "_genus_relative_abundance.tsv"))
  
  meta_path <- file.path(base_path, "Metadata", metadata)
  
  output_dir  <- file.path(base_path, "taxa_lists/gamma_diversity")
  
  #-----------------------------
  # LOAD DATA
  #-----------------------------
  
  taxa <- fread(taxa_path)
  meta <- fread(meta_path)
  
  taxa[is.na(taxa)] <- 0
  
  #-----------------------------
  # COMBINE UPSTREAM + DOWNSTREAM
  #-----------------------------
  
  ## Samples present in taxa table
  taxa_samples <- colnames(taxa)[-1]
  
  ## Keep only metadata for samples present in taxa
  meta <- meta %>%
    filter(sample_name %in% taxa_samples)
  
  processed_samples <- character()
  removed_samples   <- character()
  
  gamma_taxa <- data.frame(Genus = taxa[[1]])
  gamma_meta <- list()
  
  #==================================================
  # HANDLE 0% CONNECTIVITY SAMPLES
  #==================================================
  
  # Remove T1
  meta <- meta %>%
    filter(time != "T1")
  
  # Identify 0% connectivity samples
  zero_conn <- meta %>%
    filter(flow_rate == 0)
  
  # Process separately for each time point
  for(current_time in unique(zero_conn$time)) {
    
    zero_time <- zero_conn %>%
      filter(time == current_time)
    
    #-----------------------------------
    # Separate nutrient treatments
    #-----------------------------------
    
    zero_nutrient <- zero_time %>%
      filter(nutrient_addition == 0)
    
    nutrient_50 <- zero_time %>%
      filter(nutrient_addition == 50)
    
    nutrient_200 <- zero_time %>%
      filter(nutrient_addition == 200)
    
    
    #==================================================
    # 1. PAIR EVERY 0N POND WITH EVERY OTHER 0N POND
    #==================================================
    
    if(nrow(zero_nutrient) >= 2) {
      
      zero_pairs <- combn(
        zero_nutrient$sample_name,
        2,
        simplify = FALSE
      )
      
      for(pair in zero_pairs) {
        
        sample1 <- pair[1]
        sample2 <- pair[2]
        
        abund <- taxa[[sample1]] + taxa[[sample2]]
        
        if(sum(abund) > 0)
          abund <- abund / sum(abund)
        
        tanks <- sort(c(
          zero_nutrient$tank[
            zero_nutrient$sample_name == sample1
          ],
          zero_nutrient$tank[
            zero_nutrient$sample_name == sample2
          ]
        ))
        
        new_name <- paste0(
          tanks[1], "_",
          tanks[2], "_",
          current_time, "_0N"
        )
        
        gamma_taxa[[new_name]] <- abund
        
        # Both ponds are 0N, so either metadata row is fine
        meta_row <- zero_nutrient[
          zero_nutrient$sample_name == sample1,
        ]
        
        meta_row$sample_name <- new_name
        
        gamma_meta[[length(gamma_meta) + 1]] <- meta_row
      }
    }
    
    
    #==================================================
    # 2. RANDOMLY ASSIGN 50N PONDS TO UNIQUE 0N PONDS
    #==================================================
    
    if(nrow(nutrient_50) > 0 && nrow(zero_nutrient) > 0) {
      
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
      
      # Randomize the 50N ponds
      nutrient_50_samples <- sample(
        nutrient_50$sample_name,
        n_pairs,
        replace = FALSE
      )
      
      for(k in seq_len(n_pairs)) {
        
        sample1 <- zero_50[k]
        sample2 <- nutrient_50_samples[k]
        
        abund <- taxa[[sample1]] + taxa[[sample2]]
        
        if(sum(abund) > 0)
          abund <- abund / sum(abund)
        
        tanks <- sort(c(
          zero_nutrient$tank[
            zero_nutrient$sample_name == sample1
          ],
          nutrient_50$tank[
            nutrient_50$sample_name == sample2
          ]
        ))
        
        new_name <- paste0(
          tanks[1], "_",
          tanks[2], "_",
          current_time, "_50N"
        )
        
        gamma_taxa[[new_name]] <- abund
        
        # IMPORTANT:
        # Use metadata from the 50N pond because it
        # represents the highest nutrient treatment
        meta_row <- nutrient_50[
          nutrient_50$sample_name == sample2,
        ]
        
        meta_row$sample_name <- new_name
        
        gamma_meta[[length(gamma_meta) + 1]] <- meta_row
      }
    }
    
    
    #==================================================
    # 3. RANDOMLY ASSIGN 200N PONDS TO UNIQUE 0N PONDS
    #==================================================
    
    if(nrow(nutrient_200) > 0 && nrow(zero_nutrient) > 0) {
      
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
      
      # Randomize the 200N ponds
      nutrient_200_samples <- sample(
        nutrient_200$sample_name,
        n_pairs,
        replace = FALSE
      )
      
      for(k in seq_len(n_pairs)) {
        
        sample1 <- zero_200[k]
        sample2 <- nutrient_200_samples[k]
        
        abund <- taxa[[sample1]] + taxa[[sample2]]
        
        if(sum(abund) > 0)
          abund <- abund / sum(abund)
        
        tanks <- sort(c(
          zero_nutrient$tank[
            zero_nutrient$sample_name == sample1
          ],
          nutrient_200$tank[
            nutrient_200$sample_name == sample2
          ]
        ))
        
        new_name <- paste0(
          tanks[1], "_",
          tanks[2], "_",
          current_time, "_200N"
        )
        
        gamma_taxa[[new_name]] <- abund
        
        # Use metadata from the 200N pond because it
        # represents the highest nutrient treatment
        meta_row <- nutrient_200[
          nutrient_200$sample_name == sample2,
        ]
        
        meta_row$sample_name <- new_name
        
        gamma_meta[[length(gamma_meta) + 1]] <- meta_row
      }
    }
    
    
    #==================================================
    # MARK ALL 0% CONNECTIVITY SAMPLES AS PROCESSED
    #==================================================
    
    processed_samples <- c(
      processed_samples,
      zero_time$sample_name
    )
  }
  
  
  #==================================================
  # PROCESS ALL >0% CONNECTIVITY SAMPLES
  #==================================================
  
  for(i in seq_len(nrow(meta))){
    
    current <- meta[i, ]
    
    print(
      current[,c(
        "sample_name",
        "tank",
        "connected_tank",
        "flow_rate",
        "time"
      )]
    )
    
    sample_name <- current$sample_name
    
    #-----------------------------------
    # Skip if already processed
    #-----------------------------------
    
    if(sample_name %in% processed_samples)
      next
    
    
    #-----------------------------------
    # Find connected sample
    #-----------------------------------
    
    partner <- meta %>%
      filter(
        tank == current$connected_tank,
        time == current$time
      )
    
    
    #-----------------------------------
    # No connected sample
    #-----------------------------------
    
    if(nrow(partner) == 0){
      
      if(current$flow_rate == 0){
        
        new_name <- paste0(
          current$tank,
          "_",
          current$time
        )
        
        gamma_taxa[[new_name]] <- taxa[[sample_name]]
        
        meta_row <- current
        meta_row$sample_name <- new_name
        
        gamma_meta[[length(gamma_meta) + 1]] <- meta_row
        
      }else{
        
        removed_samples <- c(
          removed_samples,
          sample_name
        )
        
      }
      
      processed_samples <- c(
        processed_samples,
        sample_name
      )
      
      next
    }
    
    
    #-----------------------------------
    # Connected sample found
    #-----------------------------------
    
    partner_name <- partner$sample_name
    
    if(partner_name %in% processed_samples)
      next
    
    abund <- taxa[[sample_name]] + taxa[[partner_name]]
    
    if(sum(abund) > 0)
      abund <- abund / sum(abund)
    
    tanks <- sort(c(
      current$tank,
      partner$tank
    ))
    
    new_name <- paste0(
      tanks[1], "_",
      tanks[2], "_",
      current$time
    )
    
    gamma_taxa[[new_name]] <- abund
    
    meta_row <- current
    meta_row$sample_name <- new_name
    
    gamma_meta[[length(gamma_meta) + 1]] <- meta_row
    
    processed_samples <- c(
      processed_samples,
      sample_name,
      partner_name
    )
  }
  
  
  #==================================================
  # COMBINE METADATA
  #==================================================
  
  gamma_meta <- bind_rows(gamma_meta)
  
  print(removed_samples)
  length(removed_samples)
  
  if(length(removed_samples) > 0){
    
    message(
      "Removed unpaired samples:\n",
      paste(removed_samples, collapse = "\n")
    )
    
  }
  
  #-----------------------------
  # STANDARDIZE RELATIVE ABUNDANCE TO PERCENTAGES
  #-----------------------------
  
  for(j in 2:ncol(gamma_taxa)){
    
    total <- sum(gamma_taxa[[j]], na.rm = TRUE)
    
    if(total > 0){
      gamma_taxa[[j]] <- (gamma_taxa[[j]] / total) * 100
    }
    
  }
  
  #-----------------------------
  # SAVE COMBINED DATA
  #-----------------------------
  
  print(output_dir)
  print(dim(gamma_taxa))
  print(dim(gamma_meta))
  
  fwrite(
    gamma_taxa,
    file.path(output_dir,
              paste0(taxa_name, "_gamma_relative_abundance.tsv")),
    sep = "\t"
  )
  
  fwrite(
    gamma_meta,
    file.path(output_dir,
              paste0(taxa_name, "_gamma_metadata.tsv")),
    sep = "\t"
  )
  
}


#==================================================
# RUN PIPELINE FUNCTION
#==================================================

result_23S <- run_gamma_diversity_pipeline(
  "23S_total", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_23S_green_algae <- run_gamma_diversity_pipeline(
  "23S_green_algae", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_23S_cyanobacteria <- run_gamma_diversity_pipeline(
  "23S_cyanobacteria", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)



result_18S <- run_gamma_diversity_pipeline(
  "18S_total", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_18S_zooplankton <- run_gamma_diversity_pipeline(
  "18S_zooplankton", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_18S_micro_heterotrophs <- run_gamma_diversity_pipeline(
  "18S_micro_heterotrophs", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_18S_phytoplankton <- run_gamma_diversity_pipeline(
  "18S_phytoplankton", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_18S_green_algae <- run_gamma_diversity_pipeline(
  "18S_green_algae", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_18S_fungi <- run_gamma_diversity_pipeline(
  "18S_fungi", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)


result_COI <- run_gamma_diversity_pipeline(
  "COI_total", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_COI_zooplankton <- run_gamma_diversity_pipeline(
  "COI_zooplankton", "23S_LEAP22_metadata_NR.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)
