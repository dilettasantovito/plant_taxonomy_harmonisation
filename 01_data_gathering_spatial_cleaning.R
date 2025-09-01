## 01. DATA GATHERING, SPATIAL CLEANING AND VISUALISATION ##

#### Required packages ####

# Data manipulation and visualisation
library(tidyverse) # tidy data
library(viridis) # nice palettes
library(ggpubr) # arrange plots

# Biodiversity database access
library(BIEN) # BIEN
library(rgbif) # GBIF

# Dealing with dates
library(lubridate) # date manipulation

# Spatial cleaning
library(CoordinateCleaner) # clean coordinates

# Spatial data handling
library(sf) # simple features
library(terra) # spatial data analysis
library(mapdata) # geospatial data

# Auxiliary packages
library(beepr) # sound notifications
library(styler) # pretty scripts

#### DATA GATHERING ####

# Download BIEN data

# bien_it <- BIEN_occurrence_country("Italy")
# beep(5)
#
# write.csv(bien_it, file = "./data/occurrences/bien_it.csv", row.names = FALSE)

# Download GBIF data (after setting up a GBIF account on https://www.gbif.org)

#  gbif_it <- occ_download(
#    pred_in("basisOfRecord", c("OBSERVATION", "HUMAN_OBSERVATION")),
#    pred("country", "IT"),
#    pred("taxonKey", 7707728),
#    pred("hasCoordinate", TRUE),
#    pred("hasGeospatialIssue", FALSE),
#    pred_gte("year", 1950),
#    pred_gte("coordinateUncertaintyInMeters", 0),
#    pred_lte("coordinateUncertaintyInMeters", 5000),
#    format = "SIMPLE_CSV",
#    user = "your_username", pwd = "your_password", email = "your_email@example.com"
# )
#
# gbif_status <- occ_download_wait(gbif_it[1])
#
# gbif_it <- occ_download_get(gbif_it[1]) |>
#   occ_download_import()
# beep(5)
#
# write.csv(gbif_it, file = "./data/occurrences/gbif_it.csv", row.names = FALSE)

# Data downloaded in January 2025

# gbif_it <- read.csv("./data/occurrences/gbif_it.csv", header = TRUE)
# bien_it <- read.csv("./data/occurrences/bien_it.csv", header = TRUE)

# BIEN data: https://doi.org/10.5281/zenodo.16317606
# GBIF data: https://doi.org/10.15468/dl.efd2p7

# Modify the BIEN dataset to make it more similar to the GBIF dataset and remove
# records with no species name, longitude and latitude

bien_it <- bien_it |>
  filter(date_collected >= "1950-01-01") |>
  rename(
    species = scrubbed_species_binomial,
    decimalLatitude = latitude, decimalLongitude = longitude,
    eventDate = date_collected
  ) |>
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude) & !is.na(species))

# Remove GBIF records with missing species names

gbif_it <- gbif_it |>
  filter(!species == "")

# write.csv(bien_it, file = "./data/occurrences/bien_it_clean.csv", row.names = FALSE)

#### DATA VISUALISATION ####

# Download the map of Italy

italy <- map_data("italy")

# Convert the occurrences in spatial objects

gbif_it_sp <- gbif_it |>
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "epsg:3035")

bien_it_sp <- bien_it |>
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "epsg:3035")

gbif_coords <- st_coordinates(gbif_it_sp$geometry)
bien_coords <- st_coordinates(bien_it_sp$geometry)

# Create hexbin plots to get counts

bien_hex <- ggplot() +
  geom_polygon(
    data = italy, aes(x = long, y = lat, group = group),
    fill = "gray30", linewidth = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(data = as.data.frame(bien_coords), aes(x = X, y = Y), bins = 50)

gbif_hex <- ggplot() +
  geom_polygon(
    data = italy, aes(x = long, y = lat, group = group),
    fill = "gray30", linewidth = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(data = as.data.frame(gbif_coords), aes(x = X, y = Y), bins = 50)

# Extract hexbin counts

bien_counts <- ggplot_build(bien_hex)$data[[2]]$count
gbif_counts <- ggplot_build(gbif_hex)$data[[2]]$count

# Get the maximum count

global_max_count <- max(c(bien_counts, gbif_counts))

# Build the final plots

bien_hex_plot <- ggplot() +
  geom_polygon(
    data = italy, aes(x = long, y = lat, group = group),
    fill = "gray30", size = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(
    data = as.data.frame(bien_coords), aes(x = X, y = Y),
    bins = 50, color = "gray20", size = 0.2
  ) +
  scale_fill_viridis(
    option = "G", limits = c(0, global_max_count),
    direction = -1, transform = "sqrt"
  ) +
  theme_void() +
  labs(
    title = "BIEN",
    x = "decimalLongitude",
    y = "decimalLatitude",
    fill = "No of\noccurrences"
  ) +
  theme(
    legend.text = element_text(size = 8), legend.title = element_text(
      size = 9,
      vjust = 0.8
    ),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

gbif_hex_plot <- ggplot() +
  geom_polygon(
    data = italy, aes(x = long, y = lat, group = group),
    fill = "gray30", size = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(
    data = as.data.frame(gbif_coords), aes(x = X, y = Y),
    bins = 50, color = "black", size = 0.2
  ) +
  scale_fill_viridis(
    option = "G", limits = c(0, global_max_count),
    direction = -1, transform = "sqrt"
  ) +
  theme_void() +
  labs(
    title = "GBIF",
    x = "decimalLongitude",
    y = "decimalLatitude",
    fill = "Number of\noccurrences"
  ) +
  theme(
    legend.text = element_text(size = 8), legend.title = element_text(
      size = 9,
      vjust = 0.8
    ),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Arrange the plots side by side for comparison (Fig. 1)

combined_plot <- ggarrange(bien_hex_plot, gbif_hex_plot,
  ncol = 2, nrow = 1,
  common.legend = TRUE, legend = "right"
)

# Display the combined plot

print(combined_plot)

#### DATA INTEGRATION ####

# Before integrating the two datasets, we will add a column to each dataset containing an identifier,
# so that we will always be able to link the data back to the original dataset

gbif_id <- gbif_it |>
  mutate(Identifier = paste("GBIF", row_number(), sep = ""))

bien_id <- bien_it |>
  mutate(Identifier = paste("BIEN", row_number(), sep = ""))

# In BIEN dataset the collection date only includes the day, month and year of collection,
# while in GBIF it sometimes includes the time of collection as well
# Before merging the two datasets, we will have to standardise the collection time

gbif_id <- gbif_id |>
  mutate(eventDate = as.Date(eventDate))

# Now we can merge the two datasets, maintaining only the columns we strictly need

gbif_subset <- gbif_id |>
  select(species, decimalLatitude, decimalLongitude, eventDate, Identifier)

bien_subset <- bien_id |>
  select(species, decimalLatitude, decimalLongitude, eventDate, Identifier)

occurrences <- rbind(gbif_subset, bien_subset)

# write.csv(occurrences, file = "./data/occurrences/occurrences.csv", row.names = FALSE)

#### SPATIAL CLEANING ####

# We will now use package CoordinateCleaner for flagging problematic occurrences

occurrences_flags <- occurrences %>%
  mutate(
    Capitals = cc_cap(., value = "flagged"),
    Centroids = cc_cen(., value = "flagged", buffer = 100),
    Institutions = cc_inst(., value = "flagged"),
    Urban = cc_urb(., value = "flagged")
  )
# beep(5)

# Now we will plot the whole dataset and the flagged occurrences

occurrences_sp <- occurrences_flags |>
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = "epsg:3035")

flag_sp <- occurrences_sp |>
  filter(Capitals == FALSE | Centroids == FALSE | Institutions == FALSE | Urban == FALSE)

no_flag_sp <- occurrences_sp |>
  filter(Capitals == TRUE, Centroids == TRUE, Institutions == TRUE, Urban == TRUE)

flag_coords <- st_coordinates(flag_sp$geometry)
no_flag_coords <- st_coordinates(no_flag_sp$geometry)

# Create hexbin plots to get counts

flag_hex <- ggplot() +
  geom_polygon(
    data = italy, aes(x = long, y = lat, group = group),
    fill = "gray30", linewidth = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(data = as.data.frame(flag_coords), aes(x = X, y = Y), bins = 50)

no_flag_hex <- ggplot() +
  geom_polygon(
    data = italy, aes(x = long, y = lat, group = group),
    fill = "gray30", linewidth = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(data = as.data.frame(no_flag_coords), aes(x = X, y = Y), bins = 50)

# Extract hexbin counts

flag_counts <- ggplot_build(flag_hex)$data[[2]]$count
no_flag_counts <- ggplot_build(no_flag_hex)$data[[2]]$count

# Get the maximum count

global_max_count_occ <- max(c(flag_counts, no_flag_counts))

# Build the final plots

flag_hex_plot <- ggplot() +
  geom_polygon(
    data = italy, aes(x = long, y = lat, group = group),
    fill = "gray30", size = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(
    data = as.data.frame(flag_coords), aes(x = X, y = Y),
    bins = 50, color = "gray20", size = 0.2
  ) +
  scale_fill_viridis(
    option = "F", limits = c(0, global_max_count_occ),
    direction = -1, transform = "sqrt"
  ) +
  theme_void() +
  labs(
    title = "Flagged occurrences",
    x = "decimalLongitude",
    y = "decimalLatitude",
    fill = "No of flagged \noccurrences"
  ) +
  theme(
    legend.text = element_text(size = 8), legend.title = element_text(
      size = 9,
      vjust = 0.8
    ),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

no_flag_hex_plot <- ggplot() +
  geom_polygon(
    data = italy, aes(x = long, y = lat, group = group),
    fill = "gray30", size = 0.2, alpha = 0.3
  ) +
  stat_bin_hex(
    data = as.data.frame(no_flag_coords), aes(x = X, y = Y),
    bins = 50, color = "black", size = 0.2
  ) +
  scale_fill_viridis(
    option = "G", limits = c(0, global_max_count_occ),
    direction = -1, transform = "sqrt"
  ) +
  theme_void() +
  labs(
    title = "Occurrences after spatial cleaning",
    x = "decimalLongitude",
    y = "decimalLatitude",
    fill = "No of\noccurrences"
  ) +
  theme(
    legend.text = element_text(size = 8), legend.title = element_text(
      size = 9,
      vjust = 0.8
    ),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# Arrange the plots side by side for comparison

combined_plot_occ <- ggarrange(flag_hex_plot, no_flag_hex_plot,
  ncol = 2, nrow = 1,
  legend = "right"
)

# Display the combined plot (Fig. 2)

print(combined_plot_occ)

# ggarrange(bien_hex_plot, gbif_hex_plot,
#   flag_hex_plot, no_flag_hex_plot,
#   ncol = 2, nrow = 2,
#   legend = "bottom"
# )

dev.off()

#### CHECKING FLAGGED RECORDS ####

no_flags <- occurrences_flags |> 
  summarise(
    Capitals = sum(!Capitals, na.rm = TRUE),
    Centroids = sum(!Centroids, na.rm = TRUE),
    Institutions = sum(!Institutions, na.rm = TRUE),
    Urban = sum(!Urban, na.rm = TRUE)
  )

#### CLEANED OBSERVATIONS ####

# Removing duplicate records

occ_clean <- occurrences_flags |>
  filter(Capitals == TRUE, Centroids == TRUE, Institutions == TRUE, Urban == TRUE) |>
  select(!c(Capitals, Centroids, Institutions, Urban)) %>%
  cc_dupl(.)

# Saving the final dataset to .csv

# write.csv(occ_clean, "./data/occurrences/occ_clean.csv", row.names = FALSE)
