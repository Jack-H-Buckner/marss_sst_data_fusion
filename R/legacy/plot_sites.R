# ---------------------------------------------------------------
# Plot lat/lon points over a coastline, clipped to a bounding box
#
# install.packages(c("ggplot2", "sf", "rnaturalearth", "rnaturalearthdata"))
# For high-resolution coastline (needed at scales under ~100 km):
#   install.packages("rnaturalearthhires",
#                    repos = "https://ropensci.r-universe.dev", type = "source")
# ---------------------------------------------------------------

library(ggplot2)
library(sf)
library(rnaturalearth)
source("R/plot_helpers.R")

# ---- site coordinates -----------------------------------------
site_nms <- c("CB001","SF058","SF017","SF001","SF051",
           "SF0013","SF044","SF047","SF048")

sites <- read.csv("data/raw/salmon_0819_run_points.csv") %>% 
  filter(point_id %in% site_nms) %>% group_by(point_id) %>%
  summarize(lat = mean(lat), lon = mean(lon))

# ---- bounding box ---------------------------------------------
bbox <- bbox_from_points(sites$lon, sites$lat, buffer = 0.250)


# ---- coastline ------------------------------------------------
# scale = "large" (1:10m) needs rnaturalearthhires; falls back to "medium".
land <- tryCatch(
  ne_countries(scale = "large", returnclass = "sf"),
  error = function(e) {
    message("rnaturalearthhires not found; using medium resolution.")
    ne_countries(scale = "medium", returnclass = "sf")
  }
)

# Crop to a slightly padded box so the coastline reaches the plot edges,
# then let coord_sf do the final clipping.
pad  <- 0.25
land <- st_crop(
  st_make_valid(land),
  st_bbox(c(bbox["xmin"] - pad, bbox["ymin"] - pad,
            bbox["xmax"] + pad, bbox["ymax"] + pad),
          crs = st_crs(4326))
)

# ---- plot -----------------------------------------------------
p <- ggplot() +
  geom_sf(data = land, fill = "grey88", colour = "grey55", linewidth = 0.3) +
  geom_text(
    data = sites,
    aes(x = lon, y = lat, label = point_id),
    size = 2.6, alpha = 0.9
  ) +
  scale_colour_viridis_c(name = "TIR scenes", option = "magma", end = 0.85) +
  coord_sf(
    xlim   = c(bbox["xmin"], bbox["xmax"]),
    ylim   = c(bbox["ymin"], bbox["ymax"]),
    expand = FALSE
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.background = element_rect(fill = "aliceblue", colour = NA),
    panel.grid       = element_line(colour = "white", linewidth = 0.25),
    panel.border     = element_rect(colour = "grey40", fill = NA, linewidth = 0.4)
  )

print(p)