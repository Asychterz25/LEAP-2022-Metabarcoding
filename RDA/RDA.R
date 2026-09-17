#==================================================
# LIBRARIES
#==================================================

library(vegan)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)

#==================================================
# FILE PATHS
#==================================================

taxa_path <- file.path("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/taxa_Lists/cleaned/COI_zooplankton_genus_relative_abundance.tsv")

meta_path <- file.path("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/Metadata/23S_LEAP22_metadata_NR.tsv")

fig_path <- file.path("/Users/alekseisychterz/Desktop/Work/PhD/Chapter_2/figures/RDA/COI_zooplankton")

#==================================================
# LOAD DATA
#==================================================

taxa <- fread(taxa_path)
meta <- fread(meta_path)

meta <- meta %>%
  filter(SampleType == "sample")

# Remove T1
meta <- meta %>%
  filter(time != "T1")

#==================================================
# COMMUNITY MATRIX
#==================================================

# Genus = first column
# All remaining columns = samples
comm <- taxa %>%
  select(-Genus) %>%
  as.matrix()

# Transpose: samples = rows, genera = columns
comm <- t(comm)

# Sample names
rownames(comm) <- names(taxa)[-1]

# Restore genus names
colnames(comm) <- taxa$Genus

# Remove T1 samples
comm <- comm[!grepl("_T1$", rownames(comm)), , drop = FALSE]


#==================================================
# REMOVE SAMPLES WITH ZERO ABUNDANCE
#==================================================

# Identify samples where every genus has abundance = 0
zero_samples <- rownames(comm)[rowSums(comm) == 0]

cat("Samples with zero total abundance:", length(zero_samples), "\n")

if (length(zero_samples) > 0) {
  cat("Removing:\n")
  print(zero_samples)
  
  comm <- comm[rowSums(comm) > 0, , drop = FALSE]
}


#==================================================
# HELLINGER TRANSFORMATION
#==================================================

comm <- decostand(
  comm,
  method = "hellinger"
)

#==================================================
# METADATA
#==================================================

meta <- meta %>%
  filter(
    SampleType == "sample",
    time != "T1"
  )

# Clean sample names
meta$sample <- trimws(as.character(meta$sample))
rownames(comm) <- trimws(rownames(comm))


#==================================================
# DIAGNOSE MATCHING
#==================================================

cat("Community samples:", nrow(comm), "\n")
cat("Metadata samples:", nrow(meta), "\n")

cat("\nCommunity but not metadata:\n")
print(setdiff(rownames(comm), meta$sample))

cat("\nMetadata but not community:\n")
print(setdiff(meta$sample, rownames(comm)))


#==================================================
# MATCH COMMUNITY MATRIX TO METADATA
#==================================================

# Keep ONLY samples that have metadata
comm <- comm[rownames(comm) %in% meta$sample, , drop = FALSE]

# Reorder metadata to exactly match community matrix
meta <- meta[match(rownames(comm), meta$sample), ]

# Verify exact correspondence
stopifnot(identical(as.character(meta$sample), rownames(comm)))

# Check
cat("Community samples:", nrow(comm), "\n")
cat("Metadata samples:", nrow(meta), "\n")

#==================================================
# METADATA
#==================================================

meta$flow_rate <- as.numeric(meta$flow_rate)
meta$nutrient_addition <- as.numeric(meta$nutrient_addition)

meta$upstream_downstream <- factor(
  meta$upstream_downstream,
  levels = c("Up","Down")
)

meta$time <- factor(
  meta$time,
  levels = c("T2","T3","T4","T5")
)

# Remove samples with missing RDA variables

complete_samples <- complete.cases(
  meta[, c(
    "flow_rate",
    "nutrient_addition",
    "upstream_downstream",
    "time"
  )]
)

meta <- meta[complete_samples, ]

comm <- comm[complete_samples, ]

#==================================================
# PARTIAL RDA
#==================================================

rda_mod <- rda(
  comm ~
    nutrient_addition +
    flow_rate +
    upstream_downstream +
    Condition(time),
  data = meta
) 

#==================================================
# RESULTS TABLES
#==================================================

#---------------------------
# Overall model statistics
#---------------------------

sum_rda <- summary(rda_mod)

overall_table <- data.frame(
  Metric = c(
    "Total variance",
    "Conditioned variance (Time)",
    "Conditioned (%)",
    "Constrained variance",
    "Constrained (%)",
    "Unconstrained variance",
    "Unconstrained (%)",
    "Adjusted R2"
  ),
  Value = c(
    sum_rda$tot.chi,
    sum_rda$pCCA$tot.chi,
    sum_rda$pCCA$tot.chi / sum_rda$tot.chi * 100,
    sum_rda$CCA$tot.chi,
    sum_rda$CCA$tot.chi / sum_rda$tot.chi * 100,
    sum_rda$CA$tot.chi,
    sum_rda$CA$tot.chi / sum_rda$tot.chi * 100,
    RsquareAdj(rda_mod)$adj.r.squared
  )
)

print(overall_table)

#---------------------------
# Overall model ANOVA
#---------------------------

overall_model <- as.data.frame(anova(rda_mod))
overall_model$Term <- rownames(overall_model)
rownames(overall_model) <- NULL

print(overall_model)

#---------------------------
# Individual predictors
#---------------------------

term_table <- as.data.frame(anova(rda_mod, by = "term"))
term_table$Predictor <- rownames(term_table)
rownames(term_table) <- NULL

print(term_table)

#---------------------------
# Canonical axes
#---------------------------

axis_table <- as.data.frame(anova(rda_mod, by = "axis"))
axis_table$Axis <- rownames(axis_table)
rownames(axis_table) <- NULL

print(axis_table)

#---------------------------
# Variance explained by axes
#---------------------------

axis_var <- data.frame(
  Axis = c("RDA1","RDA2","RDA3"),
  Eigenvalue = rda_mod$CCA$eig,
  Percent_Constrained =
    100 * rda_mod$CCA$eig / sum(rda_mod$CCA$eig),
  Cumulative =
    cumsum(100 * rda_mod$CCA$eig / sum(rda_mod$CCA$eig))
)

print(axis_var)

write.csv(overall_table,
          file.path(fig_path, "RDA_overall_summary.csv"),
          row.names = FALSE)

write.csv(overall_model,
          file.path(fig_path, "RDA_overall_model.csv"),
          row.names = FALSE)

write.csv(term_table,
          file.path(fig_path, "RDA_predictors.csv"),
          row.names = FALSE)

write.csv(axis_table,
          file.path(fig_path, "RDA_axes.csv"),
          row.names = FALSE)

write.csv(axis_var,
          file.path(fig_path, "RDA_axis_variance.csv"),
          row.names = FALSE)

#==================================================
# EXTRACT SCORES
#==================================================

site_scores <- scores(rda_mod,
                      display = "sites",
                      scaling = 2)

species_scores <- scores(rda_mod,
                         display = "species",
                         scaling = 2)

biplot_scores <- scores(rda_mod,
                        display = "bp",
                        scaling = 2)

site_scores <- as.data.frame(site_scores)
site_scores <- bind_cols(site_scores, meta)

species_scores <- as.data.frame(species_scores)
species_scores$Genus <- rownames(species_scores)

biplot_scores <- as.data.frame(biplot_scores)
biplot_scores$Variable <- rownames(biplot_scores)

#Rename Environmental Variables
biplot_scores$Variable <- dplyr::recode(
  biplot_scores$Variable,
  nutrient_addition = "Nutrient Enrichment",
  flow_rate = "Connectivity",
  upstream_downstreamDown = "Node Position (Down)"
)

#==================================================
# TOP 10 GENERA
#==================================================

species_scores$Length <-
  sqrt(species_scores$RDA1^2 +
         species_scores$RDA2^2)

top_species <-
  species_scores %>%
  arrange(desc(Length)) %>%
  slice(1:10)

#==================================================
# COLOURS
#==================================================

nutrient_cols <- c(
  "0"   = "#440154FF",
  "25"  = "#39568CFF",
  "50"  = "#20A387FF",
  "100" = "#73D055FF",
  "200" = "#FDE725FF"
)

alpha_vals <- c(
  "0"  = 0.20,
  "10" = 0.40,
  "20" = 0.60,
  "30" = 0.80,
  "40" = 1.00
)

#==================================================
# RDA PLOT
#==================================================

p <-
  
  ggplot(site_scores,
         aes(RDA1, RDA2)) +
  
  geom_hline(yintercept = 0,
             colour = "grey80") +
  
  geom_vline(xintercept = 0,
             colour = "grey80") +
  
  geom_segment(
    data = biplot_scores,
    aes(
      x = 0,
      y = 0,
      xend = RDA1 * 0.75,
      yend = RDA2 * 0.75
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.25,"cm")),
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(
      colour = factor(nutrient_addition),
      alpha = factor(flow_rate),
      shape = upstream_downstream
    ),
    size = 3.8
  ) +
  
  geom_text_repel(
    data = biplot_scores,
    aes(
      x = RDA1 * 0.8,
      y = RDA2 * 0.8,
      label = Variable
    ),
    bg.color = "white",
    bg.r = 0.15,
    inherit.aes = FALSE,
    size = 4.5
  ) +
  
  geom_text_repel(
    data = top_species,
    aes(
      x = RDA1 * 0.8,
      y = RDA2 * 0.8,
      label = Genus
    ),
    inherit.aes = FALSE,
    colour = "grey20",
    size = 4,
    bg.color = "white",
    bg.r = 0.15,
    segment.colour = "grey60"
  ) +
  
  scale_colour_manual(values = nutrient_cols,
                      name = "Nutrient Enrichment") +
  
  scale_alpha_manual(values = alpha_vals,
                     name = "Connectivity") +
  
  scale_shape_manual(
    values = c(
      Up = 16,
      Down = 18
    ),
    name = "Node"
  ) +
  
  labs(
    x = "RDA1",
    y = "RDA2"
  ) +
  
  theme_classic(base_size = 16)

p

#==================================================
# SAVE
#==================================================

ggsave(
  file.path(fig_path, "COI_zooplankton_RDA.png"),
  p,
  width = 9,
  height = 7,
  dpi = 600
)
