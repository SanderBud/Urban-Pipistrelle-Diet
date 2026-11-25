
# Calculating, analysing, and visualising beta diversity metrics (Jaccard)
# For project Urban Bat Diet
# NIOO-KNAW
# By Sander W. J. Buddendorf, 2025-2026

# Prelims ----
source("code/00_source/packages_and_functions.R")

# Import ----
df_meta <- read_rds("data/02_processed/df_diet_meta.rds")
df <- read_rds("data/02_processed/df_diet.rds")


# R1
# Wrangling ----
taxa_matrix <- df |>
  select(sample_id, taxon, reads) |>
  pivot_wider(names_from = taxon, values_from = reads, values_fill = 0) |>
  column_to_rownames("sample_id") |>
  as.matrix()

min_reads <- 8000 # see manuscript
max_reads <- df |> group_by(sample_id) |> summarise(total_reads = sum(reads)) |> pull(total_reads) |> max()

## Jaccard, averaged over repeated rarefaction (same depth as alpha diversity)
dist_jacc <- avgdist(taxa_matrix, sample = min_reads, distfun = vegdist, dmethod = "jaccard", binary = TRUE, iterations = 500, meanfun = mean)

meta <- df_meta[match(labels(dist_jacc), df_meta$sample_id), ]
meta$standard_sampling_date <- factor(meta$standard_sampling_date)
stopifnot(identical(as.character(meta$sample_id), labels(dist_jacc)))

# NMDS ----
set.seed(123)
nmds_jacc <- metaMDS(dist_jacc, k = 2, trymax = 100)
nmds_jacc$stress  # is 0.257 which is high and should therefore be interpret with caution
stressplot(nmds_jacc)

nmds_jacc_3d <- metaMDS(dist_jacc, k = 3, trymax = 100)
nmds_jacc_3d$stress  # fall below 0.2: 0.18549
stressplot(nmds_jacc_3d)

# PERMANOVA ----
## Observed effect size (R2 of habitat when given date)
perm_obs <- adonis2(dist_jacc ~  habitat + standard_sampling_date, data = meta, by = "margin", permutations = 99)
perm_obs # be aware that the higher R2 of sampling date is influenced by the higher df. Effect is stronger for habitat (F=7.8133, which is variance explained / degrees of freedom)

### If I report on sampling date (variance explained by sampling over the season)
perm_date <- adonis2(dist_jacc ~ habitat + standard_sampling_date, data = meta, by = "margin", permutations = how(nperm = 9999, blocks = meta$site))
perm_date # valid p-value for sampling date, 

### Interaction between variance explained by habitat and sampling_date?
meta <- meta |> mutate(doy = yday(as.Date(as.character(standard_sampling_date))))

perm_int <- adonis2(dist_jacc ~ habitat * doy, data = meta, by = "margin", permutations = how(nperm = 9999, blocks = meta$site))
perm_int  # with an interaction in the model, "margin" tests only the interaction term

### Check if interaction still stands when filtering out later sampling dates (because very few samples for urban)
keep <- meta$doy <= yday(as.Date("2024-08-28"))
meta_sub <- meta[keep, ]
dist_sub <- as.dist(as.matrix(dist_jacc)[keep, keep])

perm_int_filtered <- adonis2(dist_sub ~ habitat * doy, data = meta_sub, by = "margin", permutations = how(nperm = 9999, blocks = meta_sub$site))
perm_int_filtered

### Does it still stand after adding site (because of uneven sampling per site)
adonis2(dist_sub ~ site + doy + habitat:doy, data = meta_sub, by = "margin", permutations = how(nperm = 9999, blocks = meta_sub$site))

### Visualise this
set.seed(123)
nmds_sub <- metaMDS(dist_sub, k = 3, trymax = 100)

df_traj <- scores(nmds_sub) |>
  as_tibble(rownames = "sample_id") |>
  left_join(meta_sub, by = "sample_id") |>
  group_by(habitat, standard_sampling_date) |>
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), n = n(), .groups = "drop")

ggplot(df_traj, aes(NMDS1, NMDS2, colour = habitat)) +
  geom_path(arrow = arrow(length = unit(2, "mm"))) +
  geom_point(aes(size = n)) +
  theme_classic() +
  custom_scale_colour

## Site-level permutation test: all 70 assignments of 4 sites to urban
sites <- distinct(meta, site, habitat)
combos <- combn(as.character(sites$site), 4, simplify = FALSE)

f_hab <- function(urban_sites) {
  meta$hab_perm <- ifelse(meta$site %in% urban_sites, "Urban", "Natural")
  res <- adonis2(dist_jacc ~ standard_sampling_date + hab_perm, data = meta, by = "margin", permutations = 99)
  res["hab_perm", "F"]
}

f_obs <- f_hab(as.character(sites$site[sites$habitat == "Urban"]))
f_all <- vapply(combos, f_hab, numeric(1))
p_site <- mean(f_all >= f_obs)  # minimum attainable p = 1/35 ≈ 0.0286
p_site

combos[f_all >= f_obs]

# Dispersion ----
bd_hab <- betadisper(dist_jacc, meta$habitat)
permutest(bd_hab, permutations = 9999)
bd_site <- betadisper(dist_jacc, meta$site)
permutest(bd_site, permutations = 9999)

bd_hab  # prints average distance to median per habitat
boxplot(bd_hab)

keep <- !meta$site %in% c("RKO", "LEI") # leave out sites with only one datapoint as the spread is by definition zero or these
bd_site <- betadisper(as.dist(as.matrix(dist_jacc)[keep, keep]), droplevels(factor(meta$site[keep])))
permutest(bd_site, permutations = 9999)
boxplot(bd_site)


# Visualisation NMDS----
library(patchwork)

## Scores from 3D solution
df_nmds_jacc <- scores(nmds_jacc_3d, display = "sites") |>
  as_tibble(rownames = "sample_id") |>
  left_join(df_meta, by = "sample_id") |>
  mutate(habitat = ifelse(habitat == "Rural", "Natural", "Urban"))

stress_lab <- paste0("Stress = ", round(nmds_jacc_3d$stress, 3), " (k = 3)")
perm_lab <- "R² = 0.052, p = 0.029, F = 7.81"

## Plot function for one pair of axes
plot_nmds <- function(df, x, y, labels = TRUE) {
  g <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], colour = habitat)) +
    geom_point(alpha = 0.7) +
    stat_ellipse(show.legend = FALSE) +
    theme_classic() +
    labs(x = x, y = y, colour = "Habitat") +
    custom_scale_colour
  if (labels) {
    g <- g +
      annotate("text", x = Inf, y = Inf, label = stress_lab, hjust = 1.05, vjust = 1.5, size = 3) +
      annotate("text", x = Inf, y = Inf, label = perm_lab, hjust = 1.05, vjust = 3.2, size = 3)
  }
  g
}

## Main figure: axes 1 and 2
g_nmds_12 <- plot_nmds(df_nmds_jacc, "NMDS1", "NMDS2")
g_nmds_12
ggsave("output/Fig_5_R1.png", g_nmds_12, width = 5.5, height = 4)

## Supplementary figure: all axis pairs
g_nmds_all <- plot_nmds(df_nmds_jacc, "NMDS1", "NMDS2") +
  plot_nmds(df_nmds_jacc, "NMDS1", "NMDS3", labels = FALSE) +
  plot_nmds(df_nmds_jacc, "NMDS2", "NMDS3", labels = FALSE) +
  plot_layout(guides = "collect", nrow = 1) +
  plot_annotation(tag_levels = "a")
g_nmds_all
ggsave("output/Fig_S7_nmds_3axes.png", g_nmds_all, width = 12, height = 4)




