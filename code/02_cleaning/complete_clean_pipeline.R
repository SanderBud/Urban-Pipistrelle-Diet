
# Complete cleaning pipeline: From ASVs and their reads to filtered data ready for ecological analysis
# For project Urban Bat Diet
# NIOO-KNAW
# By Sander W. J. Buddendorf, 2025-2026

# Prelims ----
source("code/00_source/packages_and_functions.R")

# Import ----
feature_table <- read.delim("data/01_raw/feature-table_BOLD_NCBI.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) |> # ASV ID with reads per sample
  select(-starts_with("LV")) # Remove samples that were not part of this experiment (From colleague)
taxonomy <- read.delim("data/01_raw/taxonomy_BOLD_NCBI.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) # Taxonomy per ASV ID

df_sample_id_raw <- read_excel("data/01_raw/20250331_DNA_extraction_samples.xlsx") |>
  filter(sample_id != "LBH82-0815-2") |>
  mutate(sample_id = ifelse(sample_id == "LBH82-0815-1", "LBH82-0815", sample_id),
         sample_code = ifelse(sample_code == "L820815-1", "L820815", sample_code))

fasta <- readDNAStringSet("data/01_raw/dna-sequences_BOLD_NCBI.fasta") 
df_asv <- data.frame(id = names(fasta),
                     seq = as.character(fasta),
                     stringsAsFactors = FALSE)

# Filtering ----
## Remove ASVs with insufficient reads relative to the sample's total reads ----
df_samples_filtered <- feature_table |>
  select(-taxonomy) |>
  pivot_longer(-ID) |>
  group_by(name) |>
  mutate(total_reads = sum(value),
         prop_reads = value/total_reads) |>
  ungroup() |> 
  filter(prop_reads > 0.0001) |>
  select(-prop_reads)

## Assign taxonomy to filtered ASVs ----
df_taxa_raw <- df_samples_filtered |> 
  left_join(taxonomy |> select(Feature.ID, Taxon), by=c("ID" = "Feature.ID")) |>
  clean_names() |>
  select(id, name, taxon, value, total_reads) |> 
  separate(taxon, 
           into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), 
           sep = ";", 
           fill = "right", 
           remove = TRUE) |> 
  mutate(kingdom = str_replace(kingdom, "d__tax=", ""),
         across(kingdom:species, ~ sub("^[a-z]__*", "", .)),
         across(kingdom:species, ~ na_if(.x, ""))) # Assigns NA if taxa is not assigned (up until species level of resolution)

## Aggregate ASVs under same taxa and fix weird species names ----
df_taxa <- df_taxa_raw |>
  mutate(kingdom = ifelse(kingdom %in% c("Metazoa", "Animalia"), "Animalia", kingdom),
         family = ifelse(family == "Hemiptera", NA, family),
         genus = ifelse(genus == "Hemiptera", NA, genus),
         species = ifelse(species == "sp.", NA, species),
         species = ifelse(species == "Anopheles claviger s.l.", "Anopheles claviger", species),
         species = ifelse(species == "Chironomus nr. pseudothummi", "Chironomus pseudothummi", species),
         species = ifelse(species == "dipterum", "Cloeon dipterum", species),
         species = ifelse(species == "Fannia sp.", NA, species),
         species = ifelse(species == "Molophilus sp.", NA, species),
         species = ifelse(species == "Planuncus s.l.", NA, species),
         species = ifelse(species == "plumosus", "Chironomus plumosus", species),
         species = ifelse(species == "PsychIntGen spACG4148", NA, species),
         species = ifelse(species == "Psychoda cf. satchelli", "Psychoda satchelli", species),
         species = ifelse(species == "Chironomus nr. pseudothummi", "Chironomus pseudothummi", species),
         order = ifelse(family == "Limoniidae", "Diptera", order)) |>
  group_by(name, across(kingdom:species), total_reads) |>
  summarise(value = sum(value), .groups = "drop") 


## Remove unassigned taxa ----
df_taxa_assigned <- df_taxa |> filter(kingdom != "Unassigned")


## Remove samples which cant be assigned to pipistrelles ----
non_pip_samples <- df_taxa_assigned |>
  filter(family == "Vespertilionidae",
         !startsWith(species, "Pipistrellus")) |>
  pull(name)

df_taxa_assigned |> # to verify that these samples dont have any Pip DNA
  filter(name %in% non_pip_samples, # one of the samples mostly contains Pipistrellus pipistrellus and a tiny bit of E. nilssonii, so dont remove this one
         family == "Vespertilionidae")

non_pip_samples_filtered <- non_pip_samples[!non_pip_samples %in% "AMV50730"] # removes 2 urban samples and 3 rural ones (so not the one mentioned above, that has more Ppip sequences)

df_pip <- df_taxa_assigned |>
  filter(!name %in% non_pip_samples_filtered)

## Remove anything not Insecta and Zygentoma, which are most likely contamination ----
df_taxa_filtered <- df_pip |>
  filter(class == "Insecta",
         order != "Zygentoma") # silverfish

## Remove one sample which was sequenced twice (with different feces pallets)
df_cleaned <- df_taxa_filtered |>
  filter(name != "L820815.2") |>
  mutate(name = ifelse(name == "L820815.1", "L820815", name))

## Find rarefaction threshold ----
### Explore
df_plot <- df_cleaned |>
  group_by(name) |>
  summarise(total_reads = sum(value))

df_plot |>
  ggplot(aes(x=total_reads)) +
  geom_histogram() +
  xlim(40000,150000)

set.seed(123)
g_rare_threshold_scatter <- df_plot |> 
  ggplot(aes(x=1, y=total_reads)) +
  geom_jitter() +
  scale_y_log10(labels = scales::comma) +
  geom_hline(yintercept = 8000, colour = "red") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "Total reads (per sample)")
g_rare_threshold_scatter
ggsave("output/suppl_materials/Fig_S5_left.png", g_rare_threshold_scatter, width = 5, height = 4) # Supplementary information


g_rare_threshold_knee <- df_plot |> 
  arrange(total_reads) |>
  ggplot(aes(x = 1:nrow(df_plot), y = total_reads)) +
  geom_line() +
  geom_hline(yintercept = 8000, colour = "red") +
  theme_bw() +
  scale_y_continuous(labels = scales::comma, limits = c(0,100000)) +
  labs(y = "Total reads (per sample)", x = "Sample #")
g_rare_threshold_knee
ggsave("output/suppl_materials/Fig_S5_right.png", g_rare_threshold_knee, width = 5, height = 4) # Supplementary information


df_plot |> # 8k read and above seems good (first 'knee'), without removing too many valuable urban samples
  arrange(total_reads) |>
  print(n = 20)

### Remove samples with less than 8k reads
df <- df_cleaned |>
  group_by(name) |> 
  mutate(total_reads = sum(value)) |>
  filter(total_reads > 8000) |>
  ungroup()

## Add metadata to samples ----
### filter sample metadata
df_sample_id <- df_sample_id_raw |>
  select(sample_id, sample_code)

### Merge with DNA data and extract metadata
df_done <- df |> 
  left_join(df_sample_id, by=c("name" = "sample_code")) |>
  select(name, sample_id, everything()) |>
  rename(sample_id_full = sample_id,
         sample_id = name,
         reads = value) |>
  separate(sample_id_full, into = c("subsite", "date_raw"), remove = F) |>
  mutate(site = substr(subsite, 1, 3),
         date_sampling = as.Date(paste0(date_raw, "24"), format = "%m%d%y"),
         habitat = ifelse(site %in% c("AMV", "LDS", "LEI", "UTR"), "urban", "rural"),
         taxon = paste(kingdom, phylum, class, order, family, genus, species, sep = "_")) |>
  select(sample_id, site, habitat, date_sampling, kingdom, phylum, class, order, family, genus, species, reads, taxon)

### Set classes for later analyses and visualisations
df_donev2 <- df_done |>
  mutate(sample_id = factor(sample_id),
         site = factor(site),
         habitat = factor(str_to_sentence(habitat)),
         taxon = factor(taxon)) |>
  group_by(sample_id) |>
  mutate(total_reads = sum(reads)) |>
  ungroup()

### Standardise dates (because some sampling moments were done in two consecutive days. Here we standardize it for visualisation purposes)
df_stand_dates <- df_donev2 


## Standardise sampling date
df_sampling_bin <- tibble(
  bin = as.Date(seq.Date(
    min(df_donev2$date_sampling), 
    max(df_donev2$date_sampling), 
    by = 14
  ))) |>
  mutate(low_end = bin - 7, 
         high_end = bin + 7)

df_bin <- df_donev2 |> # Difference between actual date and standarised date is mean/min/max: 0/0.33/2
  rowwise() |>
  mutate(
    standard_sampling_date = as.Date(df_sampling_bin$bin[which(date_sampling > df_sampling_bin$low_end &
                                                                 date_sampling <= df_sampling_bin$high_end)])) |>
  ungroup() |>
  mutate(standard_sampling_date_label = format(standard_sampling_date, "%d-%b"),
         standard_sampling_date_label = factor(standard_sampling_date_label, levels = format(sort(unique(standard_sampling_date)), "%d-%b")),
         diff_date = as.numeric(abs(date_sampling - standard_sampling_date)))

# Export ----
write_rds(df_bin |> distinct(sample_id, site, habitat, date_sampling, standard_sampling_date, standard_sampling_date_label, total_reads), "data/02_processed/df_diet_meta.rds")
write_rds(df_bin |> select(sample_id, reads, taxon), "data/02_processed/df_diet.rds")




