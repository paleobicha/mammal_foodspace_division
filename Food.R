# =============================================================
# Functional Diversity in Terrestrial Mammals across Biomes
# Full Analysis Script - Clean and Reproducible Version
# Author: Gamboa, S., Galvan, S., Sobral, M., Hernández Fernández, M. & Varela, S.
# Description:
#   This script analyzes functional diversity of terrestrial mammals
#   across global biomes using trait-based methods (diet, specialization, 
#   functional dispersion, redundancy, uniqueness, etc.).
#   Includes rarefaction, PCoA, FDis, redundancy, uniqueness,
#   and relation to productivity and specialization gradients.
# =============================================================

# -----------------------------
# 1. Load Required Libraries
# -----------------------------
library(ade4)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggConvexHull)
library(GeoRange)
library(tidyverse)
library(purrr)
library(FD)
library(xlsx)
library(reshape2)
library(scales)
library(nlme)
library(ggeffects)

# -----------------------------
# 2. Load and Preprocess Data
# -----------------------------
setwd("~/Dropbox/Mac/Documents/Sara Gamboa/Papers/mammalbiofun/script_material/")

# Load mammal and biome data
mamm_data <- read.csv("./mamm_data_terr.txt")
biome_data <- read.csv("./biome_mammals.csv", header = TRUE, sep = ";")

# Harmonize binomial names
biome_data$binomial <- sub(" ", "_", biome_data$binomial)
mamm_data$Scientific <- sub(" ", "_", mamm_data$Scientific)

# Merge and recode biomes
mamm_data <- merge(mamm_data, biome_data, by.x = "Scientific", by.y = "binomial", all.y = TRUE)
biome_labels <- c(
  "1" = "Evergreen Equatorial Rainforest", "2" = "Tropical Deciduous Woodland",
  "3" = "Subtropical Desert", "4" = "Sclerophyllous Woodland and Shrubland",
  "5" = "Temperate Evergreen Forest", "6" = "Broad-leaf Deciduous Forest",
  "7" = "Steppe", "8" = "Taiga", "9" = "Tundra", "23" = "Savanna"
)
mamm_data$BIOME <- recode(as.character(mamm_data$BIOME), !!!biome_labels)
mamm_data <- subset(mamm_data, !(BIOME %in% c("98", "99")))

# =============================================================
# 3. Define Dietary Specialization Categories
# =============================================================
diet_cols <- c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk",
               "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO")

# Remove incomplete cases
mamm_data <- mamm_data[complete.cases(mamm_data[, diet_cols]), ]

# Define generalists and specialists
mamm_data$D_General <- apply(mamm_data[, diet_cols], 1, function(x) ifelse(any(x == 100), 0, 1))
mamm_data$D_Special <- ifelse(mamm_data$D_General == 1, 0, 1)
mamm_data$Biome_specialist <- as.logical(mamm_data$Specialist)
mamm_data$D_Special <- as.logical(mamm_data$D_Special)
mamm_data$SS <- mamm_data$Biome_specialist & mamm_data$D_Special

# Calculate Biome Specialization Index (BSI)
mamm_data$BSI <- ave(mamm_data$Scientific, mamm_data$Scientific, FUN = length)

# =============================================================
# 4. Principal Coordinates Analysis (PCoA)
# =============================================================
diet_data <- prep.fuzzy(mamm_data[, diet_cols], 10, labels = "diet")
data_ktab <- ktab.list.df(list(diet_data))
mat_dissim <- dist.ktab(data_ktab, type = "F", option = "scaledBYrange")
pcoa_res <- pcoa(mat_dissim)
pcoa_axes <- as.data.frame(pcoa_res$vectors)
rownames(pcoa_axes) <- make.unique(mamm_data$Scientific)
mamm_data$Scientific_unique <- make.unique(as.character(mamm_data$Scientific))

# Save output
write.xlsx(pcoa_axes, "./output2/pcoa_mamm_terr.xlsx", row.names = TRUE)
write.xlsx(pcoa_res$values, "./output2/pcoa_values.xlsx", row.names = TRUE)

# Correlation of diet variables with PCoA axes
scores <- pcoa_axes[, 1:2]
diet_vars <- mamm_data[, diet_cols]
loading_table <- data.frame(
  Diet = colnames(diet_vars),
  Axis1 = round(apply(diet_vars, 2, function(x) cor(x, scores[,1])), 3),
  Axis2 = round(apply(diet_vars, 2, function(x) cor(x, scores[,2])), 3)
)
print(loading_table)

# =============================================================
# 5. Functional Space Occupation: Convex Hull Area per Biome
# =============================================================
# Calculate total functional space area using all species
total_area <- CHullArea(pcoa_axes$Axis.1, pcoa_axes$Axis.2)
pcoa_axes$Scientific <- rownames(pcoa_axes)

# Merge PCoA coordinates with biome information
joined_data <- mamm_data %>%
  mutate(Scientific_unique = make.unique(as.character(Scientific))) %>%
  select(Scientific_unique, BIOME) %>%
  left_join(pcoa_axes, by = c("Scientific_unique" = "Scientific"))

# Calculate functional space occupation per biome (as % of total area)
perc_hull <- joined_data %>%
  group_by(BIOME) %>%
  summarise(area = CHullArea(Axis.1, Axis.2), .groups = 'drop') %>%
  mutate(percent = round((area / total_area) * 100, 2))

# Preview results
print(perc_hull)

# -------------------------------------------------------------
# 6. Plot: Functional Space Occupied per Biome
# -------------------------------------------------------------
# Define consistent color palette
colores <- c("#2E8A3C", "#813BB0", "#1C542D", "#A43749", "#CB9162", "#EEB003",
             "#2D866A", "#A13678", "#84B43C", "#308091")

# Plot occupation percentage using barplot
perc_hull$BIOME <- factor(perc_hull$BIOME)

ggplot(perc_hull, aes(x = BIOME, y = percent)) +
  geom_col(color = colores, fill = colores, alpha = 0.6, width = 0.8) +
  scale_y_continuous(name = "% Functional Space") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, hjust = 0.5)
  ) +
  labs(title = "Percentage of Functional Space Occupied per Biome")
# -------------------------------------------------------------
# 7. PCoA Plots per Biome by Specialization Group
# -------------------------------------------------------------
# Prepare colors
biomes <- unique(mamm_data$BIOME)
biome_colors <- setNames(colores, biomes)

# Helper function to plot convex hull per biome/group
plot_biome_group <- function(biom, condition, title_suffix, color_fill) {
  biom_sel <- mamm_data$BIOME == biom
  group_sel <- condition
  biom_data_subset <- pcoa_axes[biom_sel & group_sel, ]
  
  if (nrow(biom_data_subset) > 2) {
    file_name <- paste0("./plots/PCoA_groups/", title_suffix, "_", gsub(" ", "_", biom), ".pdf")
    pdf(file = file_name, height = 8, width = 10)
    
    eigvals <- round(pcoa_res$values$Eigenvalues / sum(pcoa_res$values$Eigenvalues) * 100, 2)
    
    p <- ggplot(pcoa_axes, aes(x = Axis.1, y = Axis.2)) +
      geom_point(color = "grey90", size = 1.5) +
      geom_point(data = biom_data_subset, color = color_fill, size = 2) +
      geom_convexhull(data = biom_data_subset, alpha = 0.3, fill = color_fill) +
      ggtitle(paste0(biom, " - ", title_suffix)) +
      xlab(paste0("Axis 1 (", eigvals[1], "%)")) +
      ylab(paste0("Axis 2 (", eigvals[2], "%)")) +
      coord_equal() +
      theme_minimal(base_size = 14)
    
    print(p)
    dev.off()
  }
}

# Create output directory
if (!dir.exists("./plots/PCoA_groups")) dir.create("./plots/PCoA_groups")

# Loop by biome and group
for (biom in biomes) {
  plot_biome_group(biom, mamm_data$BSI == 1, "Specialists", "firebrick")
  plot_biome_group(biom, mamm_data$BSI >= 2 & mamm_data$BSI <= 4, "Moderate_Generalists", "dodgerblue3")
  plot_biome_group(biom, mamm_data$BSI >= 5, "Extreme_Generalists", "darkgreen")
}

# =============================================================
# 8. Rarefaction Analysis: % of Functional Space per Biome
# =============================================================
# Parameters
set.seed(42)
sample_size <- 100
iterations <- 100

# Define biome order and colors for plotting
biome_order <- c(
  "Evergreen Equatorial Rainforest", "Tropical Deciduous Woodland", "Savanna",
  "Subtropical Desert", "Sclerophyllous Woodland and Shrubland",
  "Temperate Evergreen Forest", "Broad-leaf Deciduous Forest",
  "Steppe", "Taiga", "Tundra"
)

biome_colors <- c(
  "Evergreen Equatorial Rainforest" = "#813BB0",
  "Tropical Deciduous Woodland" = "#84B43C",
  "Savanna" = "#1C542D",
  "Subtropical Desert" = "#EEB003",
  "Sclerophyllous Woodland and Shrubland" = "#A43749",
  "Temperate Evergreen Forest" = "#A13678",
  "Broad-leaf Deciduous Forest" = "#2E8A3C",
  "Steppe" = "#CB9162",
  "Taiga" = "#2D866A",
  "Tundra" = "#308091"
)

# Darker color for outlines
darken_color <- function(color, factor = 1.4) {
  rgb_vals <- col2rgb(color) / 255
  rgb_dark <- pmax(0, rgb_vals / factor)
  rgb(rgb_dark[1], rgb_dark[2], rgb_dark[3])
}
biome_lines <- sapply(biome_colors, darken_color)

# Prepare data: combine traits and PCoA
pcoa_data <- pcoa_axes
pcoa_data$Scientific <- rownames(pcoa_data)

joined <- left_join(mamm_data, pcoa_data, by = c("Scientific_unique" = "Scientific"))

# Matrix to store rarefied functional space areas
space_matrix <- matrix(NA, nrow = iterations, ncol = length(biome_order))
colnames(space_matrix) <- biome_order

# Loop for rarefaction: subsample species and compute convex hull area
for (i in 1:iterations) {
  for (j in seq_along(biome_order)) {
    biome_name <- biome_order[j]
    subset_df <- joined %>% filter(BIOME == biome_name)
    
    if (nrow(subset_df) >= 3) {
      sampled <- if (nrow(subset_df) >= sample_size) {
        subset_df[sample(nrow(subset_df), sample_size), ]
      } else {
        subset_df
      }
      
      hull <- chull(sampled$Axis.1, sampled$Axis.2)
      area <- CHullArea(sampled$Axis.1[hull], sampled$Axis.2[hull])
      space_matrix[i, j] <- (area / total_area) * 100
    }
  }
}

# Convert to long format for plotting
space_long <- melt(space_matrix, variable.name = "Biome", value.name = "Area")
names(space_long) <- c("ID", "Biome", "Area")
space_long$Biome <- factor(space_long$Biome, levels = biome_order)

# =============================================================
# 9. Visualization: Rarefied Functional Space and Empirical Comparison
# =============================================================

# Calculate empirical convex hull areas for specialists in each biome
specialist_empirical <- joined %>%
  filter(BSI == 1, BIOME %in% biome_order) %>%
  group_by(BIOME) %>%
  summarise(
    Area = if (n() > 2) (CHullArea(Axis.1, Axis.2) / total_area) * 100 else NA_real_,
    .groups = 'drop'
  ) %>%
  filter(!is.na(Area)) %>%
  rename(Biome = BIOME) %>%
  mutate(Biome = factor(Biome, levels = biome_order))

# Boxplot comparing rarefied areas with specialist empirical values
ggplot(space_long, aes(x = Biome, y = Area, fill = Biome, color = Biome)) +
  geom_boxplot(
    width = 0.6,
    size = 0.8,
    outlier.shape = 21,
    outlier.size = 1.2,
    outlier.colour = "grey30"
  ) +
  geom_point(
    data = specialist_empirical,
    aes(x = Biome, y = Area),
    inherit.aes = FALSE,
    shape = 8, size = 3.5,
    color = "black"
  ) +
  scale_fill_manual(values = biome_colors) +
  scale_color_manual(values = biome_lines) +
  ylim(0, 100) +
  labs(
    title = "Rarefied Functional Space per Biome\n(★ = Specialist Species Only)",
    y = "% of Total Functional Space",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# -----------------------------------------------
# 10. Functional Redundancy by Biome and Specialization
# -----------------------------------------------
# Round PCoA axes for discrete trait entity definition
traits_combined <- traits_combined %>%
  mutate(
    Axis1_round = round(Axis.1, 2),
    Axis2_round = round(Axis.2, 2),
    FuncEntity = paste0(Axis1_round, "_", Axis2_round)
  )

# Initialize results table
redundancy_results <- data.frame()

# Calculate redundancy by biome × specialization
for (b in unique(traits_combined$BIOME)) {
  for (s in unique(traits_combined$Specialization)) {
    subset_data <- traits_combined %>%
      filter(BIOME == b, Specialization == s)
    
    if (nrow(subset_data) >= 1) {
      n_species <- nrow(subset_data)
      n_entities <- length(unique(subset_data$FuncEntity))
      redundancy <- n_species - n_entities
      
      redundancy_results <- rbind(redundancy_results, data.frame(
        BIOME = b,
        Specialization = s,
        N_species = n_species,
        N_entities = n_entities,
        Redundancy = redundancy
      ))
    }
  }
}

# Redundancy for all species per biome (not split by group)
redundancy_all <- traits_combined %>%
  group_by(BIOME) %>%
  summarise(
    N_species = n(),
    N_entities = n_distinct(FuncEntity),
    .groups = "drop"
  ) %>%
  mutate(
    Specialization = "All",
    Redundancy = N_species - N_entities
  )

# Combine all results
redundancy_final <- bind_rows(redundancy_results, redundancy_all) %>%
  arrange(BIOME, Specialization)

# Plot: Functional Redundancy by Biome and Group
group_colors <- c(
  "Specialist" = "#FF5959",
  "Moderate Generalist" = "#FFAD5A",
  "Extreme Generalist" = "#4F9DA6",
  "All" = "gray60"
)

ggplot(redundancy_final, aes(x = Specialization, y = Redundancy, fill = Specialization)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ BIOME, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Effective Functional Redundancy by Biome and Specialization",
    y = "Redundancy (N species - N entities)",
    x = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
# =============================================================
# 11. Functional Dispersion (FDis) by Biome and Specialization
# =============================================================

# Prepare functional space matrix
traits_fspace <- pcoa_axes[, c("Axis.1", "Axis.2")]
rownames(traits_fspace) <- rownames(pcoa_axes)

# Metadata with BIOME, BSI and Specialization category
trait_metadata <- mamm_data %>%
  mutate(Scientific_unique = make.unique(as.character(Scientific))) %>%
  select(Scientific_unique, BIOME, BSI) %>%
  mutate(Specialization = case_when(
    BSI == 1 ~ "Specialist",
    BSI >= 2 & BSI <= 4 ~ "Moderate Generalist",
    BSI >= 5 ~ "Extreme Generalist"
  ))

# Combine coordinates and metadata
traits_combined <- trait_metadata %>%
  left_join(traits_fspace %>% rownames_to_column("Scientific_unique"), by = "Scientific_unique") %>%
  filter(!is.na(BIOME), !is.na(Specialization))

# Initialize FDis results table
fdis_results <- data.frame()

# Compute FDis for each biome × specialization group
for (b in unique(traits_combined$BIOME)) {
  for (s in unique(traits_combined$Specialization)) {
    subset_data <- traits_combined %>% filter(BIOME == b, Specialization == s)
    
    if (nrow(subset_data) > 2) {
      coords <- subset_data[, c("Axis.1", "Axis.2")]
      rownames(coords) <- subset_data$Scientific_unique
      
      abund <- matrix(1, nrow = 1, ncol = nrow(coords))
      colnames(abund) <- rownames(coords)
      
      fdis_val <- FD::dbFD(x = coords, a = abund, messages = FALSE)$FDis
      
      fdis_results <- rbind(fdis_results, data.frame(
        BIOME = b,
        Specialization = s,
        FDis = round(fdis_val, 4),
        N_species = nrow(subset_data)
      ))
    }
  }
}

# Compute FDis per biome using all species (not grouped)
fdis_total <- traits_combined %>%
  group_by(BIOME) %>%
  filter(n() > 2) %>%
  group_split() %>%
  purrr::map_dfr(function(df) {
    coords <- df[, c("Axis.1", "Axis.2")]
    rownames(coords) <- df$Scientific_unique
    abund <- matrix(1, nrow = 1, ncol = nrow(coords))
    colnames(abund) <- rownames(coords)
    
    tibble(
      BIOME = unique(df$BIOME),
      Specialization = "All Species",
      FDis = round(FD::dbFD(x = coords, a = abund, messages = FALSE)$FDis, 4),
      N_species = nrow(coords)
    )
  })

# Combine grouped and total results
fdis_results_all <- bind_rows(fdis_results, fdis_total)

# Save results
write.xlsx(fdis_results_all, "./output2/FDis_by_biome_group.xlsx", row.names = FALSE)

# Visualization: FDis by Specialization
group_colors <- c(
  "Specialist" = "#FF5959",
  "Moderate Generalist" = "#FFAD5A",
  "Extreme Generalist" = "#4F9DA6"
)

# Plot: Boxplot of FDis by specialization group
ggplot(fdis_results, aes(x = Specialization, y = FDis, fill = Specialization)) +
  geom_boxplot(alpha = 0.3, width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = Specialization), width = 0.15, size = 2.5, shape = 21, stroke = 0.7) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Functional Dispersion (FDis) by Specialization Group",
    x = NULL,
    y = "FDis"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"
  )

# Statistical Test: Kruskal-Wallis and Wilcoxon

kruskal_test <- kruskal.test(FDis ~ Specialization, data = fdis_results)
print(kruskal_test)

pairwise_test <- pairwise.wilcox.test(fdis_results$FDis, fdis_results$Specialization, p.adjust.method = "BH")
print(pairwise_test)

# =============================================================
# 12. Functional Uniqueness by Biome and Specialization
# =============================================================
# Compute pairwise distances using the two first PCoA axes
dist_matrix <- dist(pcoa_axes[, c("Axis.1", "Axis.2")])
dist_matrix <- as.matrix(dist_matrix)

# Compute uniqueness score: mean distance to the 5 nearest neighbors
k <- 5
uniqueness_scores <- apply(dist_matrix, 1, function(x) {
  mean(sort(x[x > 0])[1:k], na.rm = TRUE)  # exclude self-distance
})

# Add uniqueness to PCoA data
pcoa_axes$Uniqueness <- uniqueness_scores
pcoa_axes$Scientific_unique <- rownames(pcoa_axes)

# Merge with mammal metadata and define specialization groups
uniqueness_data <- mamm_data %>%
  mutate(Scientific_unique = make.unique(as.character(Scientific))) %>%
  left_join(pcoa_axes[, c("Scientific_unique", "Uniqueness")], by = "Scientific_unique") %>%
  mutate(Specialization = case_when(
    BSI == 1 ~ "Specialist",
    BSI >= 2 & BSI <= 4 ~ "Moderate Generalist",
    BSI >= 5 ~ "Extreme Generalist"
  )) %>%
  filter(!is.na(Uniqueness))

# Summarize Uniqueness by Biome and Specialization
# Mean and SD by biome
uniqueness_by_biome <- uniqueness_data %>%
  group_by(BIOME) %>%
  summarise(
    mean_uniqueness = round(mean(Uniqueness, na.rm = TRUE), 4),
    sd_uniqueness = round(sd(Uniqueness, na.rm = TRUE), 4),
    n_species = n()
  )

# Mean and SD by specialization
uniqueness_by_group <- uniqueness_data %>%
  group_by(Specialization) %>%
  summarise(
    mean_uniqueness = round(mean(Uniqueness, na.rm = TRUE), 4),
    sd_uniqueness = round(sd(Uniqueness, na.rm = TRUE), 4),
    n_species = n()
  )

# Top 10 most unique species
top_unique_species <- uniqueness_data %>%
  arrange(desc(Uniqueness)) %>%
  select(Scientific, BIOME, Specialization, Uniqueness) %>%
  slice_head(n = 10)

# Statistical Tests: Kruskal-Wallis and Pairwise Wilcoxon
# Overall test across groups
kruskal.test(Uniqueness ~ Specialization, data = uniqueness_data)

# Pairwise comparisons
pairwise.wilcox.test(
  uniqueness_data$Uniqueness,
  uniqueness_data$Specialization,
  p.adjust.method = "BH"
)

# Visualization: Uniqueness by Specialization and Biome
group_colors <- c(
  "Specialist" = "#FF5959",
  "Moderate Generalist" = "#FFAD5A",
  "Extreme Generalist" = "#4F9DA6"
)

biome_order <- c(
  "Evergreen Equatorial Rainforest", "Tropical Deciduous Woodland", "Savanna",
  "Broad-leaf Deciduous Forest", "Temperate Evergreen Forest",
  "Sclerophyllous Woodland and Shrubland", "Subtropical Desert",
  "Steppe", "Taiga", "Tundra"
)

# Prepare plotting data
plot_data <- uniqueness_data %>%
  filter(BIOME %in% biome_order) %>%
  mutate(
    Specialization = factor(Specialization, levels = c("Specialist", "Moderate Generalist", "Extreme Generalist")),
    BIOME = factor(BIOME, levels = biome_order)
  )

# Plot: Uniqueness by Specialization across biomes
ggplot(plot_data, aes(x = Specialization, y = Uniqueness, fill = Specialization)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Specialization), width = 0.15, size = 2, alpha = 0.6, show.legend = FALSE) +
  facet_wrap(~ BIOME, ncol = 2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  stat_compare_means(
    method = "wilcox.test",
    label = "letters",
    comparisons = list(
      c("Specialist", "Moderate Generalist"),
      c("Specialist", "Extreme Generalist"),
      c("Moderate Generalist", "Extreme Generalist")
    ),
    label.y = 0.0025,
    p.adjust.method = "BH"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Functional Uniqueness by Specialization Group Across Biomes",
    x = NULL,
    y = "Functional Uniqueness"
  ) +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none"
  )

# =============================================================
# 13. Relationship Between Geographic Breadth (BSI)
#     and Trophic Breadth (Number of Diet Categories)
# =============================================================
# Count number of diet categories used per species
mamm_data$diet_breadth <- apply(mamm_data[, diet_cols], 1, function(x) sum(x > 0))

# One row per species for modeling
diet_breadth_data <- mamm_data %>%
  group_by(Scientific_unique) %>%
  summarise(
    diet_breadth = first(diet_breadth),
    BSI = first(BSI),
    family = first(family),
    order = first(order)
  )

#Linear Model
lm_simple <- lm(diet_breadth ~ BSI, data = diet_breadth_data)
summary(lm_simple)

#Linear Mixed Model (family as random effect)
# Convert BSI to numeric
diet_breadth_data$BSI <- as.numeric(as.character(diet_breadth_data$BSI))

# Fit LMM with random intercept by family
lme_model <- lme(diet_breadth ~ BSI, random = ~1 | family, data = diet_breadth_data)
summary(lme_model)

#Visualization: Relationship Between BSI and Diet Breadth
ggplot(diet_breadth_data, aes(x = BSI, y = diet_breadth)) +
  geom_jitter(width = 0.3, height = 0.2, alpha = 0.5, color = "grey40") +
  geom_smooth(method = "lm", color = "#2E8A3C", fill = "#2E8A3C", alpha = 0.2, size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Biome Specialization Index (BSI)",
    y = "Trophic Breadth (Number of Diet Categories)",
    title = "Relationship Between Geographic and Trophic Breadth"
  )

# Mixed Model Predictions: BSI as factor
# Create factor version of BSI
diet_breadth_data$BSI_factor <- as.factor(diet_breadth_data$BSI)

# Model with BSI as factor
lme_bsi_cat <- lme(diet_breadth ~ BSI_factor, random = ~1 | family, data = diet_breadth_data, method = "REML")

# Predictions with confidence intervals
pred_bsi_cat <- ggpredict(lme_bsi_cat, terms = "BSI_factor")

# Plot predictions
ggplot(pred_bsi_cat, aes(x = x, y = predicted)) +
  geom_point(size = 3, color = "#2E8A3C") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, color = "#2E8A3C") +
  geom_line(group = 1, color = "#2E8A3C", linetype = "dashed") +
  labs(
    title = "Predicted Trophic Breadth by Biome Specialization Index",
    x = "Biome Specialization Index (BSI)",
    y = "Predicted Trophic Breadth"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

# Model summary table (rounded)
round(summary(lme_bsi_cat)$tTable, 3)

# =============================================================
#  14. Productivity Gradient & Functional Space Occupation
# =============================================================
# Assign ordinal productivity scores to biomes (1 = least productive, 10 = most)
productivity_rank <- tibble(
  BIOME = c(
    "Evergreen Equatorial Rainforest", "Tropical Deciduous Woodland",
    "Temperate Evergreen Forest", "Broad-leaf Deciduous Forest",
    "Savanna", "Taiga", "Sclerophyllous Woodland and Shrubland",
    "Steppe", "Tundra", "Subtropical Desert"
  ),
  Productivity = c(10, 9, 8, 7, 6, 5, 4, 3, 2, 1)
)

# Merge with percentage of convex hull space (from previous section)
space_productivity <- left_join(perc_hull, productivity_rank, by = "BIOME")

# Plot: functional space vs productivity
ggplot(space_productivity, aes(x = Productivity, y = percent)) +
  geom_point(size = 3.5, color = "#2E8A3C") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_x_continuous(breaks=seq(1:10))+
  labs(
    x = "Biome Productivity (Ordinal Scale)",
    y = "% Functional Space Occupied",
    title = "Rarefied Functional Space vs Productivity"
  ) +
  theme_minimal(base_size = 14)

# Statistical test: Spearman correlation
cor_test_result <- cor.test(space_productivity$Productivity, space_productivity$percent, method = "spearman")
print(cor_test_result)

# Rarefaction-Corrected Relationship
#Calculate average rarefied space across iterations
rarefied_means <- tibble(
  BIOME = colnames(space_matrix),
  Rarefied_Percent = colMeans(space_matrix, na.rm = TRUE)
)

#Merge with productivity rank
rarefied_productivity <- left_join(rarefied_means, productivity_rank, by = "BIOME")

#Plot rarefied functional space vs productivity
ggplot(rarefied_productivity, aes(x = Productivity, y = Rarefied_Percent)) +
  geom_point(size = 3.5, color = "#2E8A3C") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  labs(
    title = "Rarefied Functional Space vs Productivity",
    x = "Biome Productivity (Ordinal Scale)",
    y = "% Functional Space (Rarefied)"
  ) +
  theme_minimal(base_size = 14)

#Linear model and correlation test
lm_rare <- lm(Rarefied_Percent ~ Productivity, data = rarefied_productivity)
summary(lm_rare)

cor.test(rarefied_productivity$Productivity, rarefied_productivity$Rarefied_Percent, method = "spearman")

# =============================================================
# 15. Diet Composition by Biome and Specialization Group
# =============================================================

# Select columns related to diet categories
diet_cols <- c(
  "Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk",
  "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO"
)

# Classify species into specialization groups based on Biome Specialization Index (BSI)
mamm_data <- mamm_data %>%
  mutate(Specialization = case_when(
    BSI == 1 ~ "Specialist",
    BSI >= 2 & BSI <= 4 ~ "Moderate Generalist",
    BSI >= 5 ~ "Extreme Generalist"
  ))

# Compute the mean dietary composition by biome and specialization
diet_summary <- mamm_data %>%
  filter(!is.na(Specialization), !is.na(BIOME)) %>%
  group_by(BIOME, Specialization) %>%
  summarise(across(all_of(diet_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Export the summary table if needed
write.xlsx(diet_summary, "./output2/Dietary_Profile_by_Biome_and_Specialization.xlsx", row.names = FALSE)

# Define custom order for plotting
diet_order <- c(
  "Diet.PlantO", "Diet.Seed", "Diet.Fruit", "Diet.Nect",
  "Diet.Inv", "Diet.Vfish", "Diet.Vect", "Diet.Vunk",
  "Diet.Vend", "Diet.Scav"
)

biome_order <- rev(c(
  "Evergreen Equatorial Rainforest",
  "Tropical Deciduous Woodland",
  "Savanna",
  "Subtropical Desert",
  "Sclerophyllous Woodland and Shrubland",
  "Temperate Evergreen Forest",
  "Broad-leaf Deciduous Forest",
  "Steppe",
  "Taiga",
  "Tundra"
))

specialization_order <- c("Specialist", "Moderate Generalist", "Extreme Generalist")

# Reshape the data to long format for plotting
diet_long <- diet_summary %>%
  pivot_longer(cols = all_of(diet_cols), names_to = "DietCategory", values_to = "MeanPercent") %>%
  mutate(
    BIOME = factor(BIOME, levels = biome_order),
    Specialization = factor(Specialization, levels = specialization_order),
    DietCategory = factor(DietCategory, levels = diet_order)
  )

# Define a custom color palette
custom_palette <- c(
  "#14130f", "#005682", "#2a7eaa", "#6eb6c8",
  "#eae6b9", "#f4a261", "#f78513"
)

# Heatmap: Mean diet composition by biome and specialization
ggplot(diet_long, aes(x = DietCategory, y = BIOME, fill = MeanPercent)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~Specialization, ncol = 1) +
  scale_fill_gradientn(colours = custom_palette, name = "Mean % Diet") +
  labs(
    title = "Average Diet Composition by Biome and Specialization",
    x = "Diet Category",
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 16, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20)
  )

# =============================================================
# 16. PERMANOVA: Effect of Biome and Specialization on Diet Profile
# =============================================================

# Prepare diet matrix and grouping metadata
diet_individual <- mamm_data %>%
  filter(!is.na(BIOME), !is.na(Specialization)) %>%
  select(all_of(diet_cols), BIOME, Specialization)

diet_matrix <- as.matrix(diet_individual[, diet_cols])
group_data <- diet_individual[, c("BIOME", "Specialization")]

# Perform PERMANOVA to test differences in dietary profiles
permanova_result <- adonis2(
  diet_matrix ~ BIOME * Specialization,
  data = group_data,
  method = "bray",
  permutations = 999
)

# Print results
print(permanova_result)

# =============================================================
# 17. Functional Space Overlap: Nestedness & Turnover Analysis
# =============================================================

# Required libraries
library(dplyr)
library(tibble)
library(sf)

# Function to calculate nestedness and turnover between two groups
calculate_hull_overlap <- function(data, biome_name, group1_filter, group2_filter) {
  # Filter data for the given biome
  subset_data <- data %>% filter(BIOME == biome_name)
  
  # Apply dynamic filters to define groups (interpreted as logical expressions)
  group1 <- subset_data %>% filter(eval(parse(text = group1_filter)))
  group2 <- subset_data %>% filter(eval(parse(text = group2_filter)))
  
  # At least 3 points are required to form a convex hull
  if (nrow(group1) < 3 | nrow(group2) < 3) {
    return(tibble(
      Biome = biome_name,
      Intersection_Area = NA,
      Group1_Area = NA,
      Group2_Area = NA,
      Nestedness = NA,
      Turnover = NA
    ))
  }
  
  # Extract convex hull coordinates for both groups
  coords1 <- group1[chull(group1$Axis.1, group1$Axis.2), c("Axis.1", "Axis.2")]
  coords2 <- group2[chull(group2$Axis.1, group2$Axis.2), c("Axis.1", "Axis.2")]
  
  # Ensure polygons are closed (first point repeated at the end)
  coords1 <- rbind(coords1, coords1[1, ])
  coords2 <- rbind(coords2, coords2[1, ])
  
  # Create spatial polygons (in functional trait space, no CRS)
  poly1 <- sf::st_sfc(sf::st_polygon(list(as.matrix(coords1))), crs = sf::NA_crs_)
  poly2 <- sf::st_sfc(sf::st_polygon(list(as.matrix(coords2))), crs = sf::NA_crs_)
  
  # Calculate intersection polygon
  intersection <- suppressWarnings(sf::st_intersection(poly1, poly2))
  
  # Compute polygon areas
  a1 <- as.numeric(sf::st_area(poly1))
  a2 <- as.numeric(sf::st_area(poly2))
  ai <- if (length(intersection) > 0) as.numeric(sf::st_area(intersection)) else 0
  
  # Nestedness = proportion of group 1 contained in group 2
  nestedness <- ai / a1
  turnover <- 1 - (ai / (a1 + a2 - ai))
  
  return(tibble(
    Biome = biome_name,
    Intersection_Area = ai,
    Group1_Area = a1,
    Group2_Area = a2,
    Nestedness = nestedness,
    Turnover = turnover
  ))
}

# Prepare dataset for overlap analysis
# Combine biome, BSI and PCoA axes with species identity
traits_joined <- mamm_data %>%
  mutate(Scientific_unique = make.unique(as.character(Scientific))) %>%
  select(Scientific_unique, BIOME, BSI) %>%
  mutate(Specialization = case_when(
    BSI == 1 ~ "Specialist",
    BSI >= 2 & BSI <= 4 ~ "Moderate_Generalist",
    BSI >= 5 ~ "Extreme_Generalist"
  )) %>%
  left_join(
    pcoa_axes %>% rownames_to_column("ID_temp") %>%
      mutate(Scientific_unique = make.unique(ID_temp)),
    by = "Scientific_unique"
  ) %>%
  filter(!is.na(BIOME), !is.na(Specialization))

# List of unique biomes
biomes <- unique(traits_joined$BIOME)

# Run overlap function: Specialists vs Moderate Generalists
results_nestedness_mod <- purrr::map_dfr(
  biomes,
  ~calculate_hull_overlap(traits_joined, .x, "BSI == 1", "BSI >= 2 & BSI <= 4")
)

# Run overlap function: Specialists vs Extreme Generalists
results_nestedness_ext <- purrr::map_dfr(
  biomes,
  ~calculate_hull_overlap(traits_joined, .x, "BSI == 1", "BSI >= 5")
)

# Summarize Functional Uniqueness by Specialization
uniqueness_summary <- uniqueness_data %>%
  group_by(Specialization) %>%
  summarise(
    Mean_Uniqueness = round(mean(Uniqueness, na.rm = TRUE), 4),
    SD_Uniqueness = round(sd(Uniqueness, na.rm = TRUE), 4),
    N_Species = n()
  )

# Print summary table
print(uniqueness_summary)




