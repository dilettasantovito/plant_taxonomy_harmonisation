## 02. TAXONOMIC HARMONISATION ##

#### Required packages ####

# Data manipulation and visualisation
library(tidyverse)

# Auxiliary packages
library(styler) # pretty scripts
library(beepr) # sound notifications
library(tictoc) # timings

#### Data loading ####

# gbif_it <- read.csv("./data/occurrences/gbif_it.csv", header = TRUE)
# bien_it <- read.csv("./data/occurrences/bien_it_clean.csv", header = TRUE)
occ_clean <- read.csv("./data/occurrences/occ_clean.csv", header = TRUE)

# Renaming the column species in order to better distinguish input names and
# output names

occ_clean <- occ_clean |>
  rename(Input_species = species)

# Total number of unique species

unique_species <- occ_clean |>
  distinct(Input_species)

# Saving the data in .csv format

# write.csv(unique_species, file = "./data/occurrences/unique_species.csv", row.names = F)

#### lcvplants ####

# This package lets you harmonise plant names by fuzzy matching and synonymy
# resolution against the Leipzig Plant Catalogue as taxonomic backbone

# Installation
# devtools::install_github("idiv-biodiversity/LCVP") # data package
# devtools::install_github("idiv-biodiversity/lcvplants")

library(lcvplants)

# Matching

# lcvp_search() is the core function of the package. Connects to the LCVP table
# and validates the names of a vector of plant taxa. The algorithm will first try
# to exactly match the binomial names provided. If no match is found, it will try
# to find the closest name given the maximum distance defined in max_distance.

# From now on, I'll take the time required by the different packages to work,
# to have a general idea about the performance of the various functions

tic()

lcvp_species <- unique_species$Input_species |>
  lcvp_search(show_correct = TRUE)

lcvp_time <- toc()

# To summarise the results from a multiple species search and know the number of
# species matched, and how many of them were exactly or fuzzy matched, there's
# the lcvp_summary() function.

lcvp_summary(lcvp_species)

# Species searched: 5902
# Species matched: 5878 (99.59%)
# Species exactly matched: 5796 (98.2%)
# Species fuzzy matched: 82 (1.39%)
# Authors exactly matched: 0 (0%)
# Infracategories exactly matched: 5875 (99.54%)

# Some species were matched multiple times

# lcvp_sps_mult <- attr(lcvp_species, "matched_mult")
# lcvp_sps_mult

# We can have a look at it using the function lcvp_fuzzy_search()

# lcvp_fuzzy_match <- lcvp_sps_mult |>
#  lcvp_fuzzy_search()
# view(lcvp_fuzzy_match)

# Since the output includes authors, only keep the genus and specific epithet

lcvp_species_parsed <- lcvp_species |>
  mutate(
    Output.Author = word(Output.Taxon, start = 3),
    Output.Taxon = word(Output.Taxon, 1, 2)
  )

lcvp_notmatched <- lcvp_species_parsed |>
  filter(is.na(Output.Taxon)) # non matched species

lcvp_notmatched_short <- lcvp_notmatched |> 
  select(Search)

# write.csv(lcvp_notmatched, "./data/output/lcvp_notmatched.csv", row.names = FALSE)

# lcvp_fuzzy <- lcvp_fuzzy_search(unique_species$species) # to know every possible match

# N.B: decide wether to keep the 'unresolved' names

#### rWCVP ####

# rWCVP is a package designed to make some common tasks that use WCVP easier,
# including the standardisation of taxon names against the World Checklist of
# Vascular Plants

# Installation

# devtools::install_github("matildabrown/rWCVP")
library(rWCVP)

# In order to make rWCVP functions work, we have to install the associated data
# package, rWCVPdata

# if (!require(rWCVPdata)) {
#  install.packages("rWCVPdata",
#    repos = c(
#      "https://matildabrown.github.io/drat",
#      "https://cloud.r-project.org"
#    )
#  )
# }

library(rWCVPdata)
data(wcvp_names) # loading the world checklist

# Matching names against the WCVP checklist
# The package offers several functions for name matching:
# - wcvp_match_exact for exact matching
# - wcvp_match_fuzzy for fuzzy matching
# - wcvp_match_names for both exact and fuzzy matching

# Matching names using wcvp_match_names

tic()

wcvp_species <- wcvp_match_names(unique_species,
  wcvp_names,
  name_col = "Input_species"
)

# ✔ Matched 6099 of 6103 names
# ℹ Exact (without author): 6064
# ℹ Fuzzy (edit distance): 13
# ℹ Fuzzy (phonetic): 22
# ! Names with multiple matches: 1152

rWCVP_time <- toc()

wcvp_accepted_species <- wcvp_species |>
  dplyr::filter(
    wcvp_status == "Accepted",
    match_similarity >= 0.8
  )

wcvp_accepted_species_short <- wcvp_accepted_species |>
  select(Input_species, wcvp_name)

#### taxadb ####

# taxadb provides fast, consistent access to taxonomic data. Existing approaches
# to these problems typically relies on APIs. taxadb creates a local database
# of most readily available taxonomic authorities

# Installation

# install.packages("taxadb")

library(taxadb)

# Authorities currently recognized by taxadb (list found in package documentation):

# - itis: Integrated Taxonomic Information System
# - ncbi: National Center for Biotechnology Information
# - col: Catalogue of Life
# - gbif: Global Biodiversity Information Facility
# - ott: OpenTreeTaxonomy
# - iucn: IUCN Red List
# - itis_test

# Retrieving the taxonomic identifiers
# Create a vector of names containing all taxadb providers

taxadb_prov <- c("itis", "col", "gbif") # excluding ott and iucn because
# I receive an error

# ott: Error: no match found for schema: common provider: ott version: 22.12
# iucn: Error: no match found for schema: dwc provider: iucn version: 22.12

# Also excluding NCBI because it doesn't work with get_names()

# Create a local taxonomic database

td_create(provider = taxadb_prov, dbdir = "./data/taxonomy", overwrite = TRUE)

# Retrieve taxonomic IDs and accepted names, separately for every provider

taxadb_times <- vector(mode = "list", length = length(taxadb_prov))
names(taxadb_times) <- taxadb_prov

taxadb_species <- purrr::map(taxadb_prov, \(x) {
  tic()

  result <- unique_species |>
    mutate(id = taxadb::get_ids(Input_species, provider = x)) |>
    mutate(accepted_name = taxadb::get_names(id, provider = x))

  taxadb_times[[x]] <<- toc()

  return(result)
})

# Retrieving matched non matched names

taxadb_match <- purrr::map(taxadb_species, \(x) x |> filter(!is.na(accepted_name)))
taxadb_nonmatch <- purrr::map(taxadb_species, \(x) x |> filter(is.na(accepted_name)))

# Resolving multiple matches
# Since taxadb doesn't have an automated way to resolve multiple matches,
# we'll just keep the accepted names match the input names

taxadb_nonmatch_multi <- purrr::map2(taxadb_nonmatch, taxadb_prov, \(x, y)
filter_name(x$Input_species, provider = y) |>
  mutate(accepted_name = get_names(acceptedNameUsageID, provider = y)) |>
  select(scientificName, accepted_name, acceptedNameUsageID) |>
  rename(Input_species = scientificName, id = acceptedNameUsageID) |>
  distinct(id, .keep_all = TRUE) |>
  mutate(check = Input_species == accepted_name) |>
  filter(check == TRUE))

# Now that we got the accepted names that match the input names, we want to check
# for duplicated values

taxadb_duplicated <- purrr::map(taxadb_nonmatch_multi, \(x)
x |>
  group_by(Input_species) |>
  filter(n() > 1) |>
  summarize(count = n()) |>
  arrange(desc(count)))

# For some reason there are different IDs that give back the same accepted names,
# so we'll just keep the unique values

taxadb_nonmatch_multi <- purrr::map(taxadb_nonmatch_multi, \(x) x |>
  distinct(accepted_name, .keep_all = TRUE) |>
  select(-check))

# Final taxadb matches tables

taxadb_species_final <- purrr::map2(taxadb_match, taxadb_nonmatch_multi, rbind)

# Extracting the dataframes from the list

taxadb_df_names <- c("taxadb_itis", "taxadb_col", "taxadb_gbif")
names(taxadb_species_final) <- taxadb_df_names

taxadb_species_final |>
  purrr::map(\(x) x |>
    select(-id)) |>
  list2env(.GlobalEnv)

#### taxonomycleanR and taxize ####

# The taxonomyCleanr is a user friendly workflow and collection of functions
# to help you: (1) identify and correct misspelled taxa, (2) resolve taxa to an
# authority, (3) get hierarchical rank values, (4) get common names, (5) render
# taxonomic information in the Ecological Metadata Language (EML)
# The taxonomyCleanr package in based on the taxize package developed by
# Chamberlain et al. (2016). That's why we will be using the two together.

# Installation

# remotes::install_github("ropensci/taxize")
# remotes::install_github("EDIorg/taxonomyCleanr")

library(taxonomyCleanr)
library(taxize)

view_taxa_authorities() # taxa authorities recognized by taxonomyCleanr

# id                                       authority resolve_sci_taxa resolve_comm_taxa
# 1   3  Integrated Taxonomic Information System (ITIS)        supported         supported
# 2   9        World Register of Marine Species (WORMS)        supported     not supported
# 3  11 Global Biodiversity Information Facility (GBIF)        supported     not supported
# 4 165            Tropicos - Missouri Botanical Garden        supported     not supported

# In order to use tropicos we need a usage key that can be obtained via the taxize package
# (that's why I put the two packages together)

# use_tropicos() # this function is via the taxize package
# TROPICOS_KEY: "your_tropicos_key"

tropicoskey <- "your_tropicos_key"

Sys.setenv(TROPICOS_KEY = tropicoskey)

# taxonomycleanR

tc_sources <- c(3, 11, 165) # vector containing the codes associated with the
# authorities I need to standardize plant names

tc_times <- vector(mode = "list", length = length(tc_sources))

tc_names <- c("tc_itis", "tc_gbif", "tc_tropicos")

names(tc_times) <- tc_names

tc_species <- purrr::map2(tc_sources, tc_names, \(sources, names) # creating a loop to use this package
                          # core functions for every different chosen authority
                          {
                            tic()

                            my_path <- tempdir()

                            taxa_map <- create_taxa_map(path = my_path, x = unique_species, col = "Input_species")

                            tc_count <- count_taxa(x = unique_species, col = "Input_species", path = my_path)

                            tc_trim <- trim_taxa(path = my_path)

                            tc_resolve <- resolve_sci_taxa(path = my_path, data.sources = sources)

                            tc_times[[names]] <<- toc()

                            return(tc_resolve)
}
)
beep(5)

tc_species_accepted <- purrr::map(tc_species, \(x) x |>
                                    filter(!is.na(authority_id)) |>
                                    rename(Input_species = taxa_raw) |>
                                    select(Input_species, taxa_clean))

names(tc_species_accepted) <- tc_names

list2env(tc_species_accepted, .GlobalEnv)

#### taxize (not included in the analysis) ####

# GBIF
# taxize_gbif_ids <- unique_species |>
#   mutate(id = get_gbifid(unique_species$Input_species, rows = 1)) # done
#
# # Plants of the World
# #taxize_pow_ids <- unique_species |>
# #  mutate(id = get_pow(unique_species$Input_species, rows = 1)) # too many requests
#
# # Tropicos
# taxize_tropicos_ids <- unique_species |>
#   mutate(id = get_tpsid(unique_species$Input_species, rows = 1)) # done
#
# # ITIS
# #taxize_itis_ids <- unique_species |>
# #  mutate(id = get_tsn(unique_species$Input_species, rows = 1))
#
# #! SSL peer certificate or SSH remote key was not OK: [www.itis.gov]
# #SSL certificate problem: unable to get local issuer certificate
#
# #taxize_gbif_species <- id2name(as.vector(taxize_gbif_ids$id), db = "gbif")
# # Error in data.frame(id = x, name = z$canonicalName, rank = tolower(z$rank),  :
# #arguments imply differing number of rows: 1, 0
#
# #taxize_tropicos_species <- id2name(taxize_tropicos_ids$id, db = "tropicos") # not an option
#
# # Let's try the gnr_resolve function
#
# gnr_ds <- gnr_datasources()
# gnr_ds
# gnr_usefulds <- c(1, 165, 196, 197, 198)
# # Catalogue of Life (1), GBIF (11), Tropicos (165),
# # WorldFloraOnline (196), World Checklist of Vascular Plants (197), Leipzig Catalogue
# # of Vascular Plants (198)
#
# # Excluded GBIF (11) because it gives me the following error
# # Error: Internal Server Error (HTTP 500)
#
# tic()
#
# taxize_gnr_species <- map(gnr_usefulds, \(x)
#                           gnr_resolve(unique_species$Input_species,
#                                       data_source_ids = x,
#                                       canonical = TRUE, best_match_only = TRUE
#                           ))
#
# taxize_gnr_time <- toc()
#
# taxize_gnr_names <- c(
#   "taxize_col", "taxize_tropicos", "taxize_wfo", "taxize_wcvp",
#   "taxize_lcvp"
# )
#
# names(taxize_gnr_species) <- taxize_gnr_names
#
# taxize_gnr_species |>
#   #map(\(x) x |>
#   #      mutate(matched_author = word(matched_name2, start = 3),
#   #             matched_species = word(matched_name2, 1, 2))) |>
#   map(\(x) x |>
#         rename(Input_species = user_supplied_name)) |>
#   map(\(x) x |>
#         select(c(matched_name2, Input_species))) |>
#   list2env(.GlobalEnv)
#
# # Filtering species with a match score > 0.8
#
# taxize_matched_names <- c(
#   "taxize_col_matched",
#   "taxize_tropicos_matched",
#   "taxize_wfo_matched",
#   "taxize_wcvp_matched",
#   "taxize_lcvp_matched"
# )
#
# taxize_matched <- taxize_gnr_species |>
#   map(\(x) x |> filter(score >= 0.8)) |>
#   map(\(x) x |> rename(Input_species = user_supplied_name)) |>
#   map(\(x) x |> select(c(matched_name2, Input_species)))
#
# names(taxize_matched) <- taxize_matched_names
#
# taxize_matched |>
#   list2env(.GlobalEnv)
#
# # Filtering species with a match score < 0.8
#
# taxize_notmatched_names <- c(
#   "taxize_col_notmatched",
#   "taxize_tropicos_notmatched",
#   "taxize_wfo_notmatched",
#   "taxize_wcvp_notmatched",
#   "taxize_lcvp_notmatched"
# )
#
# taxize_notmatched <- taxize_gnr_species |>
#   map(\(x) x |> filter(score < 0.8))
#
# names(taxize_notmatched) <- taxize_notmatched_names
#
# taxize_notmatched |>
#   list2env(.GlobalEnv)

#### taxizedb ####

# taxizedb provides tools for working with taxonomic databases on your machine,
# similarly to taxadb

# Taxonomic databases supported:
# - NCBI
# - ITIS
# - COL
# - GBIF
# - Wikidata
# - World Flora Online

# Installation

# remotes::install_github("ropensci/taxizedb")
library(taxizedb)

# Downloading taxonomic databases

# NCBI
# db_download_ncbi(overwrite = TRUE)
# Error in curl::curl_download(db_url, db_path_file, quiet = TRUE) :
#  Timeout was reached: [] Connection timed out after 10000 milliseconds

# db_download_col(overwrite = TRUE)
# Error in curl::curl_download(db_url, db_path_file, quiet = TRUE) :
# HTTP error 403.

db_download_gbif(overwrite = TRUE) # "~/.cache/R/taxizedb/gbif.sqlite"
db_download_itis(overwrite = TRUE) # "~/.cache/R/taxizedb/ITIS.sqlite"
# db_download_wfo(overwrite = TRUE)
# Error in curl::curl_download(db_url, db_path, quiet = TRUE) :
#  Timeout was reached: [] Connection timed out after 10001 milliseconds

# db_download_wikidata(overwrite = TRUE)
# Downloaded but we don't need that, because I can't use it for retrieving
# scientific names

# Downloaded databases: gbif, itis
# Retrieving species ids

taxizedb_dbs <- c("itis", "gbif")
taxizedb_times <- vector(mode = "list", length = length(taxizedb_dbs))
names(taxizedb_times) <- taxizedb_dbs

taxizedb_species <- purrr::map2(
  taxizedb_dbs,
  purrr::map(taxizedb_dbs, \(x) {
    tic()
    result <- name2taxid(unique_species$Input_species, db = x, out_type = "summary")
    taxizedb_times[[x]] <<- toc()
    result
  }),
  \(x, y) {
    y |>
      mutate(newnames = taxid2name(y$id, db = x)) |>
      rename(Input_species = name)
  }
)

taxizedb_df_names <- c("taxizedb_itis", "taxizedb_gbif")
names(taxizedb_species) <- taxizedb_df_names

taxizedb_species |>
  purrr::map(\(x) x |>
    select(-id) |>
    distinct_all()) |>
  list2env(.GlobalEnv)

#### TNRS ####

# The TNRS package provides access to the Taxonomic Name Resolution Service API,
# which is a tool for automated standardization of plant scientific names. The
# TNRS corrects spelling errors and alternative spellings to a standard list of
# names

# Installation
# install.packages("TNRS")
library(TNRS)

# Selected databases

# TNRS_dbs <- c("tropicos", "wfo", "wcvp") # Tropicos not working

TNRS_dbs <- c("wfo", "wcvp")

TNRS_times <- vector(mode = "list", length = length(TNRS_dbs))
names(TNRS_times) <- TNRS_dbs

TNRS_species <- purrr::map(TNRS_dbs, \(x) {
  tic()

  result <- TNRS(unique_species$Input_species, sources = x) |>
    filter(Overall_score >= 0.8)

  TNRS_times[[x]] <<- toc()

  return(result)
})

TNRS_names <- c("TNRS_wfo", "TNRS_wcvp")

names(TNRS_species) <- TNRS_names

TNRS_species |>
  list2env(.GlobalEnv)

TNRS_wcvp_short <- TNRS_wcvp |>
  select(Name_submitted, Accepted_name) |>
  rename(Input_species = Name_submitted)

TNRS_wfo_short <- TNRS_wfo |>
  select(Name_submitted, Accepted_name) |>
  rename(Input_species = Name_submitted)

TNRS_wcvp_notmatched <- TNRS(unique_species$Input_species, sources = "wcvp") |> 
  filter(Accepted_name == "" | Overall_score < 0.8)

TNRS_wcvp_notmatched_short <- TNRS_wcvp_notmatched |> 
  select(Name_submitted)
  
#### WorldFlora ####

# Installation

# install.packages("WorldFlora")
library(WorldFlora)

# First you need to download and unzip the World Flora Online taxonomic backbone from
# www.worldfloraonline.org/downloadData
# Downloading version v.2024.06

# Load WorldFlora classification

wfo_classification <- read.table("./data/taxonomy/wfo_classification.csv",
  header = TRUE, sep = "\t", quote = "", comment.char = "",
  strip.white = TRUE
)

# Match the species list with WorldFlora classification

tic()

wfo_species <- WFO.match(unique_species$Input_species, WFO.data = wfo_classification)

wfo_time <- toc()

beep(5)

# Matching one-to-one

wfo_species_one <- WFO.one(wfo_species, priority = "Accepted") |> 
  mutate(across(everything(), ~ str_remove_all(., '"'))) |>
  mutate(across(everything(), ~ replace(.x, is.na(.x), ""))) |>
  mutate(Output_species = paste(genus, specificEpithet, sep = " ")) |>
  map_dfr(na_if, " ") |>
  filter(Fuzzy.dist <= 4)

wfo_species_short <- wfo_species_one |>
  select(spec.name, Output_species) |>
  rename(Input_species = spec.name)

wfo_unmatched <- wfo_species_short |> 
  filter(is.na(Output_species))

# In order: removed '"' in the df, converted NAs in blank spaces to avoid the
# conversion of NAs in characters, converted blanks spaces in NAs again

#### Matches table ####

accepted_species <- list(
  unique_species, taxadb_col, taxadb_gbif, taxadb_itis, taxizedb_gbif, taxizedb_itis,
  tc_gbif, tc_itis, tc_tropicos, wcvp_accepted_species_short, 
  TNRS_wcvp_short, TNRS_wfo_short,wfo_species_short 
  ) |>
  reduce(left_join, by = "Input_species") |>
  cbind(
    lcvp_species_parsed$Output.Taxon
  ) |>
  mutate_all(na_if, "")

colnames(accepted_species) <- c(
  "Input_species", "taxadb_COL", "taxadb_GBIF", "taxadb_ITIS", "taxizedb_GBIF", "taxizedb_ITIS",
  "taxonomyCleanr_GBIF", "taxonomyCleanr_ITIS", "taxonomyCleanr_TROPICOS", "rWCVP",
  "TNRS_WCVP", "TNRS_WFO", "WorldFlora", "lcvplants"
)

# Save the output in .csv format

# write.csv(accepted_species, "./data/output/accepted_species.csv", row.names = FALSE)