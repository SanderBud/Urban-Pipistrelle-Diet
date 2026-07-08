
# Calculating, analysing, and visualising beta diversity metrics (Jaccard)
# For project Urban Bat Diet
# NIOO-KNAW
# By Sander W. J. Buddendorf, 2025-2026

# Prelims ----
source("code/00_source/packages_and_functions.R")

# Import ----
df_meta <- read_rds("data/02_processed/df_diet_meta.rds")
df <- read_rds("data/02_processed/df_diet.rds")

# Wrangling ----
min_reads <- df |> group_by(sample_id) |> summarise(total_reads = sum(reads)) |> pull(total_reads) |> min()
max_reads <- df |> group_by(sample_id) |> summarise(total_reads = sum(reads)) |> pull(total_reads) |> max()

df_wide <- df |>
  pivot_wider(names_from = taxon, values_from = reads, values_fill = 0) |>
  as.data.frame()

rownames(df_wide) <- df_wide$sample_id
matrix_wide <- as.matrix(df_wide[, -1])

## Jaccard
set.seed(123)
matrix_pa <- decostand(matrix_wide, method = "pa")
dist_jacc <- vegdist(matrix_pa, method = "jaccard", binary = TRUE)

nmds_jacc <- metaMDS(dist_jacc)
df_nmds_jacc <- scores(nmds_jacc) |>
  as_tibble(rownames = "sample_id") |>
  left_join(df_meta, by = "sample_id")

# Visualisation ----
g_nmds_habitat_jacc <- df_nmds_jacc |>
  mutate(habitat = ifelse(habitat == "Rural", "Natural", "Urban")) |>
  ggplot(aes(x=NMDS1, y=NMDS2, colour=habitat)) +
  geom_point() +
  theme_classic() +
  stat_ellipse(show.legend = F) +
  labs(colour = "Habitat") +
  custom_scale_colour
g_nmds_habitat_jacc
ggsave("output/Fig_4.png", g_nmds_habitat_jacc, width = 5.5, height = 4)


# Stats ----
## Site nested in habitat to see if most variance is explained by habitat or site
adonis2(dist_jacc ~ standard_sampling_date + habitat/site, data = df_nmds_jacc, 
        permutations = 9999, by = "terms") # about equal, both significant. Habitat still important

