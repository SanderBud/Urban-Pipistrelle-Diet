
# Code for bar graph of relative read abundance, divided by habitat.
# For project Urban Bat Diet
# NIOO-KNAW
# By Sander W. J. Buddendorf, 2025-2026

# Prelims ----
source("code/00_source/packages_and_functions.R")

# Import ----
df_meta <- read_rds("data/02_processed/df_diet_meta.rds")
df_raw <- read_rds("data/02_processed/df_diet.rds") 

# Wrangle ----
df_rra <- df_raw |>
  separate(taxon, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = "_") |>
  left_join(df_meta, by = "sample_id") |>
  filter(family != "NA") |>
  group_by(habitat) |>
  mutate(reads_habitat = sum(reads)) |>
  group_by(habitat, order, reads_habitat) |>
  summarise(reads_taxa = sum(reads), .groups = "drop") |>
  mutate(rra = reads_taxa / reads_habitat,
         habitat = ifelse(habitat == "Rural", "Natural", "Urban"))

order_sort <- df_rra |> filter(habitat == "Urban") |> arrange(desc(rra)) |> pull(order)


# Visualise ----
## Colour
base_colors <- c(
  "#9198FF", 
  "#91cfff", 
  "#66d1a8", 
  "#a6d96a", 
  "#d9ef8b",
  "#fee08b",
  "#fcbe61", 
  "#fcae91", 
  "#fc8d59", 
  "#f46d43"
)

## Visualisation
g_rra <- df_rra |>
  ggplot(aes(x = habitat, y = rra, fill = factor(order, levels = order_sort))) +
  geom_col() +
  labs(y = "Relative read abundance (RRA)", 
       x = "Habitat", 
       fill = "Prey order") +
  scale_fill_manual(values = base_colors) +
  theme_classic()
g_rra
ggsave("output/Fig_2.png", g_rra, width = 4.5, height = 4)
