
# Code for donut graph, showing the taxonomic composition of the common pipistrelle’s diet. Originally adapted from Curran et al. (2022)
# For project Urban Bat Diet
# NIOO-KNAW
# By Sander W. J. Buddendorf, 2025-2026

# Prelims ----
source("code/00_source/packages_and_functions.R")
get_packages(c("forcats", "showtext"))

# Import data ----
df <- read_rds("data/02_processed/df_diet.rds") |>
  separate(taxon, into = c("kingdom","phylum","class","order","family","genus","species"), sep = "_")

# Wrangle ----
## Prepare presence table and summary at order/family (FOO = number of samples)
order_sort <- df |> filter(family != "NA") |> distinct(sample_id, order, family) |> count(order, name="n_lvl1")

dfb <- df |>
  filter(family != "NA") |>
  distinct(sample_id, order, family) |>
  group_by(lev.1 = order, lev.2 = family) |>
  summarise(n_lvl2 = n(), .groups = "drop") |>
  left_join(order_sort, by = c("lev.1"="order")) |>
  mutate(value = 1)

## Create "empty" root layer (keeps donut center)
lvl0 <- tibble(lev.1 = "Diet", value = 0, level = 0, fill = NA, num = NA, ang = NA, label = NA)

## Inner layer: orders
lvl1 <- dfb |>
  group_by(lev.1, n_lvl1, value) |>
  summarise(value = sum(value), .groups = "drop") |>
  arrange(desc(n_lvl1)) |>
  mutate(level = 1, fill = lev.1, num = NA, ang = NA, label = lev.1)

## Outer layer: families (replacing your previous genus layer)
lvl2 <- dfb |>
  arrange(desc(n_lvl1), desc(n_lvl2)) |>
  mutate(level = 2, fill = lev.2)

## compute ordering numbers and base angles for outer ring entries
lvl2 <- lvl2 |>
  mutate(num = rev(row_number()))   

angcoef <- 360 / nrow(lvl2)
lvl2 <- lvl2 |>
  mutate(ang = num * angcoef,
         label = paste0(lev.2, " ", n_lvl2))

## flip labels on the right-hand side so they remain readable (like your original logic)
mid <- round(nrow(lvl2)/2) + 1
lvl2 <- lvl2 |>
  mutate(
    ang = ifelse(row_number() >= mid, ang + 180, ang),
    label = ifelse(row_number() >= mid, paste0(n_lvl2, " ", lev.2), label)
  )

## Combine layers
comb_d <- bind_rows(lvl0, lvl1, lvl2) |>
  mutate(level = as.factor(level), fill = as.factor(fill))

## Preserve factor order as appears in data (avoid alphabetical reorder)
ord <- unique(comb_d$fill)
comb_d$fill <- factor(comb_d$fill, levels = ord)

## Compute order-level angles as mean of their families' angles so inner labels sit centered
orders <- unique(lvl2$lev.1)
for (o in orders) {
  inds <- which(lvl2$lev.1 == o)
  if (length(inds) > 0) {
    mean_ang <- mean(lvl2$ang[inds])
    comb_d$ang[which(comb_d$level == "1" & comb_d$lev.1 == o)] <- mean_ang
  }
}



## for root (level 0) no label/angle needed
comb_d$label[comb_d$level == "1" & is.na(comb_d$label)] <- comb_d$lev.1[comb_d$level == "1"]

# Create color mapping ----
unique_orders <- unique(lvl1$lev.1)

color_palette <- c()
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

## Manual colour pallette
for (i in seq_along(unique_orders)) {
  order_name <- unique_orders[i]
  base_color <- base_colors[i]
  
  ### Orde gets darker version 
  color_palette[order_name] <- base_color
  
  ### Family lighter version
  families_in_order <- lvl2$lev.2[lvl2$lev.1 == order_name]
  
  for (fam in families_in_order) {
    color_palette[as.character(fam)] <- base_color
  }
}

# Visualise ----
p <- ggplot(comb_d, aes(x = level, y = value, fill = fill, alpha = level))

p_final <- p +
  geom_col(width = 1, color = "gray40", linewidth = 0.25, position = position_stack()) +
  ## order labels (inner ring)
  geom_text(data = subset(comb_d, level == "1"),
            aes(label = ifelse(lev.1 %in% c("Blattodea", "Ephemeroptera"), "", label), 
                angle = ang + 88),
            position = position_stack(vjust = 0.5),
            size = 13) +
  geom_text(data = subset(comb_d, level == "1" & lev.1 == "Ephemeroptera"),
            aes(label = label, angle = ang + 88),
            position = position_stack(vjust = 0.8),
            size = 10) +
  geom_text(data = subset(comb_d, level == "1" & lev.1 == "Blattodea"),
            aes(label = label, angle = ang + 88),
            position = position_stack(vjust = 0.5),
            size = 8) +
  ## family labels (outer ring)
  geom_text(data = subset(comb_d, level == "2"),
            aes(label = label, angle = ang + 88),
            position = position_stack(vjust = 0.5),
            size = 8, alpha = 1) +
  coord_polar(theta = "y", direction = -1) +
  scale_alpha_manual(values = c("0" = 0, "1" = 1, "2" = 0.7), guide = "none") +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  scale_fill_manual(values = color_palette) +
  theme_void() +
  theme(legend.position = "none")
p_final
ggsave("output/Fig_1.png", p_final, width=25, height=25, dpi=100)


