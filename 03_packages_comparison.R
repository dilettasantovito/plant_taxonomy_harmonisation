## 03 PACKAGES COMPARISON ##

#### Required packages ####

# Data manipulation and visualisation
library(tidyverse)

# Data analysis
library(vegan) 

# Auxiliary packages
library(styler) # pretty script
library(ggrepel) # nice labels
library(viridis) # nice palettes
#library(ggplotify) # convert plots in ggplot2 objects

#### Data transformations ####

# Loading data

accepted_species <- read.csv("./data/output/accepted_species.csv", header = TRUE)

#### Summary table ####

summary_table <- accepted_species |> 
  summarise(across(everything(), \(x) sum(!is.na(x)))) |> 
  pivot_longer(everything(), names_to = "Package", values_to = "Accepted_names") |> 
  arrange(desc(Accepted_names)) |> 
  mutate(Accepted_names_perc = paste(round(Accepted_names/max(Accepted_names)*100, 2), " %"))

# Exporting the summary table

# write.csv(summary_table, "./data/output/summary_table.csv", row.names = FALSE)

#### Data exploration ####

# Add a column containing the number of single values per row
# (excluding missing values)

accepted_species_values <- accepted_species |> 
  mutate(
    Number_of_values =
      apply(accepted_species, 1, \(x) length(unique(na.omit(x))))
  ) |> 
  arrange(desc(Number_of_values))
# Consider a way to exclude Input_species from this count

# Add a column containing the number of NAs per row

matched_species_NAs <- accepted_species_values |> 
  (\(df) df |> mutate(Na_values = rowSums(is.na(df))))() |>
  arrange(desc(Na_values))

#### Building a "community" matrix ####

# In order to measure the concordance between the different packages, we'll start
# building a "community" matrix that we will subsequently use for the calculation
# of dissimilarities

matches_matrix <- accepted_species |> 
  pivot_longer(!(Input_species)) |> 
  pivot_wider(
    names_from = Input_species,
    values_from = value
  ) |> 
  column_to_rownames(var = "name") |> 
  mutate_all(~ ifelse(is.na(.), 0, 1))

matches_matrix_t <- matches_matrix |> 
  t() |> 
  as.data.frame()

#### NMDS ####

# Computing the distance using Jaccard

matches_MDS1 <- metaMDS(matches_matrix, distance = "jaccard", k = 2, trymax = 1000)

# Extract NMDS coordinates from the points element

nmds_coords <- as.data.frame(matches_MDS1$points)

# Adding point identifiers 

nmds_coords$ID <- rownames(matches_matrix)

# Plotting with ggplot2

ggplot(nmds_coords, aes(x = MDS1, y = MDS2, label = ID)) +
  geom_point() +
  geom_text(vjust = -0.5, hjust = 0.5) +
  labs(x = "NMDS1", y = "NMDS2",  font = "serif", size = 4) +
  theme_minimal()

ggplot(nmds_coords, aes(x = MDS1, y = MDS2, label = ID)) +
  geom_point(size = 2, shape = 21, alpha = 0.3, fill = "black") +
  geom_text_repel(
    size = 4,
    # family = "serif",
    max.overlaps = Inf,      # ensure all labels are shown
    box.padding = 0.5,       # space around text
    point.padding = 0.3,     # space between text and point
    segment.size = 0.2       # line connecting point to label
  ) +
  labs(x = "NMDS1", y = "NMDS2", font = "serif", size = 4) +
  theme_minimal()

# There are some outliers
# Removing them and repeating the analysis

matches_matrix_noout <- accepted_species |> 
  select(-c(taxadb_ITIS, taxonomyCleanr_ITIS, taxizedb_ITIS, taxizedb_GBIF, taxonomyCleanr_TROPICOS)) |>   
  pivot_longer(!(Input_species)) |> 
  pivot_wider(
    names_from = Input_species,
    values_from = value
  ) |> 
  column_to_rownames(var = "name") |> 
  mutate_all(~ ifelse(is.na(.), 0, 1))

set.seed(42)

matches_MDS1_noout <- metaMDS(matches_matrix_noout, distance = "jaccard", k = 2, trymax = 1000)

nmds_coords_noout <- as.data.frame(matches_MDS1_noout$points)

nmds_coords_noout$ID <- rownames(matches_matrix_noout)

# mypal <- viridis_pal(option = "turbo", begin = 0, end = 0.8)(nrow(matches_matrix_noout))

nmds_plot <- ggplot(nmds_coords_noout, aes(x = MDS1, y = MDS2, label = ID)) +
  geom_point(size = 2, shape = 21, alpha = 0.3, fill = "black") +
  #geom_jitter(size = 4, shape = 21, stroke = 1, width = 0.001, height = 0.001) +
  geom_text_repel() +
  labs(x = "NMDS1", y = "NMDS2", font = "serif", size = 4) +
  #scale_color_viridis_d(option = "turbo", begin = 0, end = 0.8) +
  theme_minimal()

nmds_plot
# ggsave("./data/output/nmds_plot.png", dpi = 320)

#### Intersection plot ####

library(UpSetR) # intersection plots

# Creating an intersection plot

mypal <- mako(ncol(accepted_species) - 1, end = 0.8)

intersect_plot <- upset(matches_matrix_t, # data matrix
                      nsets = ncol(matches_matrix_t), # number of your packages
                      sets =  colnames(matches_matrix_t), # names of your packages
                      nintersects = 20,
                      order.by = "freq", # order of intersection bars
                      text.scale = 1.2,
                      sets.bar.color = mypal
                      ) 

# intersect_plot <- as.ggplot(intersect_plot)

intersect_plot
