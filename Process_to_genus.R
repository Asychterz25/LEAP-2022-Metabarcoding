############## Combine Technical Replicates and Process Samples to Genus Level ###################

#==================================================
# LIBRARIES
#==================================================

library(data.table)

#==================================================
# MAIN PIPELINE FUNCTION
#==================================================

run_process_pipeline <- function(taxa_name,
                                 metadata,
                                 base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2",
                                 save_output = TRUE){
  
  #-----------------------------
  # FILE PATHS
  #-----------------------------
  
  taxa_path <- file.path(
    base_path,
    "taxa_lists",
    taxa_name,
    paste0(taxa_name, "_combined_taxa_aggregated.tsv")
  )
  
  meta_path <- file.path(
    base_path,
    "Metadata",
    metadata
  )
  
  output_path <- file.path(
    base_path,
    "taxa_lists",
    "cleaned"
  )
  
  #-----------------------------
  # LOAD DATA
  #-----------------------------
  
  taxa <- fread(taxa_path)
  meta <- fread(meta_path)
  
  #==================================================
  # PROCESS TAXA
  #==================================================
  
  process_taxa <- function(taxa_df){
    
    taxa_df <- as.data.table(taxa_df)
    
    #----------------------------------
    # Remove unwanted columns
    #----------------------------------
    
    taxa_df[, c("False_Positives",
                "Non-Target_Taxa") := NULL]
    
    taxa_df <- taxa_df[
      is.na(Reason) | trimws(Reason) == ""
    ]
    
    taxa_df[, Reason := NULL]
    
    #----------------------------------
    # Keep only genus + sample columns
    #----------------------------------
    
    tax_cols <- c(
      "domain","supergroup","division","subdivision",
      "class","order","family","genus","species",
      "kingdom","phylum","marker"
    )
    
    sample_cols <- setdiff(names(taxa_df), tax_cols)
    
    taxa_df <- taxa_df[
      ,
      c("genus", sample_cols),
      with = FALSE
    ]
    
    #----------------------------------
    # Remove unknown genera
    #----------------------------------
    
    taxa_df <- taxa_df[
      genus != "Unknown" &
        !is.na(genus)
    ]
    
    #----------------------------------
    # Sum duplicate genera
    #----------------------------------
    
    taxa_df <- taxa_df[
      ,
      lapply(.SD, sum),
      by = genus
    ]
    
    #----------------------------------
    # Convert to matrix
    #----------------------------------
    
    taxa_mat <- as.matrix(taxa_df[, -1])
    
    rownames(taxa_mat) <- taxa_df$genus
    
    #----------------------------------
    # Combine technical replicates
    #----------------------------------
    
    sample_names <- colnames(taxa_mat)
    
    base_names <- sub("_R[0-9]+$", "", sample_names)
    
    taxa_counts <- sapply(unique(base_names), function(x){
      
      cols <- which(base_names == x)
      
      rowSums(taxa_mat[, cols, drop = FALSE])
      
    })
    
    # Preserve matrix if only one sample
    taxa_counts <- as.matrix(taxa_counts)
    
    #----------------------------------
    # Recalculate relative abundance
    #----------------------------------
    
    taxa_rel <- sweep(
      taxa_counts,
      2,
      colSums(taxa_counts),
      "/"
    ) * 100
    
    return(taxa_rel)
    
  }
  
  taxa_rel <- process_taxa(taxa)
  
  #==================================================
  # SAVE OUTPUT
  #==================================================
  
  if(save_output){
    
    out_file <- file.path(
      output_path,
      paste0(taxa_name, "_genus_relative_abundance.tsv")
    )
    
    out_df <- data.frame(
      Genus = rownames(taxa_rel),
      taxa_rel,
      check.names = FALSE
    )
    
    fwrite(
      out_df,
      out_file,
      sep = "\t"
    )
    
    message("Saved to: ", out_file)
    
  }
  
  return(taxa_rel)
  
}

#==================================================
# RUN PIPELINE FUNCTION
#==================================================

result_23S <- run_process_pipeline(
  "23S_total", "23S_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_23S_green_algae <- run_process_pipeline(
  "23S_green_algae", "23S_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_23S_cyanobacteria <- run_process_pipeline(
  "23S_cyanobacteria", "23S_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)



result_18S <- run_process_pipeline(
  "18S_total", "18S_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_18S_zooplankton <- run_process_pipeline(
  "18S_zooplankton", "18S_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_18S_phytoplankton <- run_process_pipeline(
  "18S_phytoplankton", "18S_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_18S_micro_heterotrophs <- run_process_pipeline(
  "18S_micro_heterotrophs", "18S_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_18S_green_algae <- run_process_pipeline(
  "18S_green_algae", "18S_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_18S_fungi <- run_process_pipeline(
  "18S_fungi", "18S_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)


result_COI <- run_process_pipeline(
  "COI_total", "COI_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)

result_COI_zooplankton <- run_process_pipeline(
  "COI_zooplankton", "COI_LEAP22_metadata.tsv",
  base_path = "/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2"
)
