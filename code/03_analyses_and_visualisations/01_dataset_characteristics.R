
# Script for some dataset info for in the manuscript
# For project Urban Bat Diet
# NIOO-KNAW
# By Sander W. J. Buddendorf, 2025-2026

# Prelims ----
source("code/00_source/packages_and_functions.R")
get_packages("writexl")

# Import ----
df_meta <- read_rds("data/02_processed/df_diet_meta.rds")
df <- read_rds("data/02_processed/df_diet.rds") |>
  separate(taxon, into = c("kingdom","phylum","class","order","family","genus","species"), sep = "_")

# Sample collection overview ----
## Sampled sites per sampling week
sampled_sites_week <- df_meta |> 
  mutate(week_nr = strftime(date_sampling, format = "%V"))|>
  count(site, habitat, week_nr) |>
  arrange(week_nr) |>
  pivot_wider(names_from = week_nr, values_from = n, values_fill = 0) |>
  arrange(habitat)
sampled_sites_week
write_xlsx(sampled_sites_week, "output/suppl_materials/Table_S1.xlsx")

## Taxa 
taxa <- df |> distinct(across(kingdom:species))
nrow(taxa) # 747 taxa
taxa |> distinct(species) # 651 species
taxa |> distinct(genus) # 402 genus

## Reads and ASVs per sample
df |> left_join(df_meta, by=c("sample_id")) |>
  mutate(habitat=ifelse(habitat=="Rural", "Natural", "Urban"),
         taxon=paste(order, family, genus, species, sep="_")) |>
  group_by(sample_id, site, habitat, date_sampling) |>
  summarise(total_reads=dplyr::first(total_reads),
            n_asvs=n_distinct(taxon),
            .groups="drop") |>
  arrange(habitat, site, date_sampling) |>
  write_csv("output/suppl_materials/Table_S2.csv")



