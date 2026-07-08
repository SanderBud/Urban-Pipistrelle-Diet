
# Calculation of amount of impervious surface within a 2km range of the sites ('urbanness'), for the suppl materials. Note, this script is computationally expensive and slow
# For project Urban Bat Diet
# NIOO-KNAW
# By Sander W. J. Buddendorf, 2025-2026

# Prelims ----
source("code/00_source/packages_and_functions.R")
get_packages(c("terra", "sf", "readxl"))

# Import ----
sites <- read_excel("data/01_raw/sampled_sites_coords.xlsx")
lgn <- rast("data/01_raw/LGN2023.tif") |> project("EPSG:28992")

# Wrangle ----
## Sites to spatial points
sites_sf <- sites |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(28992)   # RD New

## Setting LGN classes
impervious_classes <- c(
  8, # glastuinbouw = greenhouses
  18, # bebouwing primair gebied = built-up primary
  19, # bebouwing secundair gebied = built-up secondary
  26, # bebouwing buitengebied = built-up rural
  251, # Wegen, hoofd = main roads
  252, # wegen, half-verhard = semi-roads
  253 # wegen, smal = small roads
)

# Get buffer ----
sites_buf <- sites_sf |>
  st_buffer(2000) |>     # buffer in meters
  vect()

ext <- extract(lgn, sites_buf)

impervious_prop <- ext |>
  as.data.frame() |>
  filter(!is.na(LGN2023)) |>
  mutate(impervious = LGN2023 %in% impervious_classes) |>
  group_by(ID) |>
  summarise(prop_imp = sum(impervious) / n(), .groups = "drop")

sites$ID <- 1:47

sites_comp <- sites |>
  left_join(impervious_prop, by = "ID")

# Visualise ----
set.seed(2)
g_habitat_point <- sites_comp |>
  group_by(site, habitat) |>
  summarise(mean_prop = mean(prop_imp),
            .groups = "drop") |>
  ggplot(aes(x = habitat, y = mean_prop, colour = habitat)) +
  geom_jitter(width = 0.2) +
  theme_classic() +
  custom_scale_colour +
  labs(y = "Mean proportion impervious surfaces \n(within a radius of 2km)") +
  theme(axis.title.x = element_blank(),
        legend.position = "none")
g_habitat_point
ggsave("output/suppl_materials/Fig_S3.png", g_habitat_point, width = 2.5, height = 3.2)



