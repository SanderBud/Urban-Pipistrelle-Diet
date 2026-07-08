
# Calculating, analysing, and visualising alpha diversity metrics (richness and shannon index)
# For project Urban Bat Diet
# NIOO-KNAW
# By Sander W. J. Buddendorf, 2025-2026

# Prelims ----
source("code/00_source/packages_and_functions.R")

# Import ----
df_meta <- read_rds("data/02_processed/df_diet_meta.rds")
df <- read_rds("data/02_processed/df_diet.rds")


# Richness ----
## Visualisation raw richness data ----
min_reads <- df |> group_by(sample_id) |> summarise(total_reads = sum(reads)) |> pull(total_reads) |> min()
max_reads <- df |> group_by(sample_id) |> summarise(total_reads = sum(reads)) |> pull(total_reads) |> max()

taxa_data <- df |> 
  select(sample_id, taxon, reads) |>
  pivot_wider(names_from = "taxon", values_from = "reads", values_fill = 0) |>
  column_to_rownames("sample_id")

### Calculating rarefied prey richness per sample
df_richness_rarefied <- taxa_data |>
  rarefy(sample = min_reads) |>
  as_tibble(rownames = "sample_id")

### Rarefaction does not really change the outcome
df_richness <- df_richness_rarefied |>
  left_join(df_meta, by=c("sample_id"))

### Visualisation per habitat
set.seed(123)
g_rich_habitat <- df_richness |>
  mutate(habitat = ifelse(habitat == "Rural", "Natural", "Urban")) |>
  ggplot(aes(x=habitat, y=value, fill = habitat)) +
  geom_boxplot(colour = "black", width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.5, ) +
  theme_classic() +
  labs(y = "Prey richness") +
  custom_scale_fill +
  theme(axis.title.x = element_blank(),
        legend.position = "none") 
g_rich_habitat
ggsave("output/Fig_3_left.png", g_rich_habitat, width = 2.8, height = 3.5)


## Models richness ----
df_rich <- df_richness |>
  rename(richness=value)

hist(df_rich$richness)
wilcox.test(richness ~ habitat, data = df_rich)

### Linear model
mod_rich_null <- glmmTMB(richness ~ 1 + (1|standard_sampling_date), data = df_rich, family = gaussian())
mod_rich_habitat <- glmmTMB(richness ~ habitat + (1|standard_sampling_date), data = df_rich, family = gaussian())
mod_rich_site <- glmmTMB(richness ~ site + (1|standard_sampling_date), data = df_rich, family = gaussian())

anova(mod_rich_null, mod_rich_habitat, mod_rich_site)

anova(mod_rich_site, mod_rich_habitat) # habitat does not improve the model, site is more imporant

df_annova_richness <- anova_summary(
  models = list(
    null = mod_rich_null,
    habitat = mod_rich_habitat,
    site = mod_rich_site
  ), 
  subject = "Richness"
)
write_csv(df_annova_richness, "output/Table_1_richness.csv")

# Shannons index ----
## Visualisation raw shannon index ----
### Calculating shannon index from 500 resamplings
shannon_iteration <- function() {
  taxa_data |> 
    rrarefy(sample = min_reads) |>
    diversity()
}

set.seed(123)
df_shannon_rare_iter <- replicate(500, shannon_iteration()) |>
  as_tibble(rownames="sample_id", .name_repair="unique") |>
  pivot_longer(-sample_id)

#### Visualising the stability of rarefied Shannon diversity estimates across 500 subsampling iterations
g_shannon_rarefied_bootstrapped <- df_shannon_rare_iter |>
  group_by(sample_id) |>
  summarise(mean_shannon=mean(value),
            se_shannon=sd(value)) |>
  ungroup() |>
  mutate(sample_id=reorder(sample_id, mean_shannon),
         panel=ifelse(as.integer(sample_id) <= 70, "1", "2")) |>
  ggplot(aes(x=sample_id, y=mean_shannon)) +
  geom_point(size=0.8) +
  geom_errorbar(aes(ymin=mean_shannon-se_shannon, ymax=mean_shannon+se_shannon), width=0.1) +
  facet_wrap(~panel, ncol=1, scales="free") +
  labs(x="Sample", y="Shannon diversity (rarefied)") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        strip.text=element_blank())
g_shannon_rarefied_bootstrapped
ggsave("output/suppl_materials/Fig_S6.png", g_shannon_rarefied_bootstrapped, width=10, height=10)

### Summarising bootstrapped data
df_shannon_rare <- df_shannon_rare_iter |>
  group_by(sample_id) |>
  summarise(shannon = mean(value))|>
  left_join(df_meta, by=c("sample_id"))

### Visualise shannon index per habitat
set.seed(123)
g_shannon_habitat <- df_shannon_rare |>
  mutate(habitat = ifelse(habitat == "Rural", "Natural", "Urban")) |>
  ggplot(aes(x=habitat, y=shannon, fill=habitat)) +
  geom_boxplot(colour = "black", width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.5) +
  theme_classic() +
  labs(y = "Shannon index") +
  custom_scale_fill +
  theme(axis.title.x = element_blank(),
        legend.position = "none") 

g_shannon_habitat
ggsave("output/Fig_3_right.png", g_shannon_habitat, width = 2.8, height = 3.5)

## Model Shannon index ----
ggplot(df_shannon_rare, aes(x = shannon)) +
  geom_histogram() +
  facet_wrap(~habitat)

shapiro.test(df_shannon_rare$shannon[df_shannon_rare$habitat == "Urban"])
shapiro.test(df_shannon_rare$shannon[df_shannon_rare$habitat == "Rural"]) # Normally distributed

#### Linear model
mod_sha_null <- glmmTMB(shannon ~ 1 + (1|standard_sampling_date), data = df_shannon_rare, family = gaussian())
mod_sha_habitat <- glmmTMB(shannon ~ habitat + (1|standard_sampling_date), data = df_shannon_rare, family = gaussian())
mod_sha_site <- glmmTMB(shannon ~ site + (1|standard_sampling_date), data = df_shannon_rare, family = gaussian())

anova(mod_sha_null, mod_sha_habitat, mod_sha_site)
anova(mod_sha_habitat, mod_sha_site)

df_annova_shannon <- anova_summary(
  models = list(
    null = mod_sha_null,
    habitat = mod_sha_habitat,
    site = mod_sha_site
  ), 
  subject = "Shannon"
)
write_csv(df_annova_shannon, "output/Table_1_shannon.csv")





