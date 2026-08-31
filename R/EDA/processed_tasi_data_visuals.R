library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
source("R/plot_helpers.R")

dat <- readRDS("outputs/tasi_salmon_20_yrs_2026-08-29/clean_data.rds")

out_dir <- "outputs/processed_tasi_data_eda"


# Plot correlation between high resolution TIR instruments and existing 
# level 4 data products
# state with all data 
dat %>% 
  filter(variable %in% c("eco_sst_v002_clean", "mur_sst","insitu_sst", 
                         "modis_sst", "lst_sst_clean")) %>%
  pivot_wider(names_from = "variable", values_from= "value") %>%
  pivot_longer(cols = c("lst_sst_clean","eco_sst_v002_clean"),
               names_to = "tir_instrument", values_to = "tir_sst") %>% 
  pivot_longer(cols = c("modis_sst", "mur_sst","insitu_sst"),
               names_to = "reference", values_to = "ref_sst") %>%
  ggplot(aes(x=ref_sst, y = tir_sst))+
    geom_point(size = 0.5, alpha = 0.25)+theme_classic(base_size = 14)+
    facet_grid(tir_instrument~reference)+
    geom_abline(aes(slope = 1, intercept = 0))+
    xlab("reference SST estiamte") + 
    ylab("High RES SST")+
  ggtitle("All high-res + reference matchups")
  
ggsave("outputs/processed_tasi_data_eda/matchup_correlations.png",
         height = 4, width = 6)


# plot the data at the monitoring station 
dat %>% 
  filter(variable %in% c("eco_sst_v002_clean", "mur_sst","insitu_sst", 
                         "modis_sst", "lst_sst_clean"),
         point_id == "CB002") %>%
  pivot_wider(names_from = "variable", values_from= "value") %>%
  pivot_longer(cols = c("lst_sst_clean","eco_sst_v002_clean"),
               names_to = "tir_instrument", values_to = "tir_sst") %>% 
  pivot_longer(cols = c("modis_sst", "mur_sst","insitu_sst"),
               names_to = "reference", values_to = "ref_sst") %>%
  ggplot(aes(x=ref_sst, y = tir_sst))+
  geom_point(size = 0.5, alpha = 0.25)+theme_classic(base_size = 14)+
  facet_grid(tir_instrument~reference)+
  geom_abline(aes(slope = 1, intercept = 0))+
  xlab("reference SST estiamte") + 
  ylab("High RES SST")+
  ggtitle("High-res reference matchups at monitoring site") 

ggsave("outputs/processed_tasi_data_eda/matchup_correlations_CB002.png",
         height = 4, width = 6)

# plot at a site with MODIS data 
dat %>% 
  filter(variable %in% c("eco_sst_v002_clean", "mur_sst", 
                         "modis_sst", "lst_sst_clean"),
         point_id == "SF041") %>%
  pivot_wider(names_from = "variable", values_from= "value") %>%
  pivot_longer(cols = c("lst_sst_clean","eco_sst_v002_clean"),
               names_to = "tir_instrument", values_to = "tir_sst") %>% 
  pivot_longer(cols = c("modis_sst", "mur_sst"),
               names_to = "reference", values_to = "ref_sst") %>%
  ggplot(aes(x=ref_sst, y = tir_sst))+
  geom_point(size = 0.5, alpha = 0.25)+theme_classic(base_size = 14)+
  facet_grid(tir_instrument~reference)+
  geom_abline(aes(slope = 1, intercept = 0))+
  xlab("reference SST estiamte") + 
  ylab("High RES SST")+
  ggtitle("High-res matchups at \n farm site with MODIS") 

ggsave("outputs/processed_tasi_data_eda/matchup_correlations_SF041.png",
         height = 4, width = 4.5)

# plot at a site without MODIS data 
dat %>% 
  filter(variable %in% c("eco_sst_v002_clean", "mur_sst", "lst_sst_clean"),
         point_id == "SF001") %>%
  pivot_wider(names_from = "variable", values_from= "value") %>%
  pivot_longer(cols = c("lst_sst_clean","eco_sst_v002_clean"),
               names_to = "tir_instrument", values_to = "tir_sst") %>% 
  pivot_longer(cols = c("mur_sst"),
               names_to = "reference", values_to = "ref_sst") %>%
  ggplot(aes(x=ref_sst, y = tir_sst))+
  geom_point(size = 0.5, alpha = 0.25)+theme_classic(base_size = 14)+
  facet_grid(tir_instrument~reference)+
  geom_abline(aes(slope = 1, intercept = 0))+
  xlab("Reference SST estiamte") + 
  ylab("High RES SST")+
  ggtitle("High-res matchups at \n farm site without MODIS") 

ggsave("outputs/processed_tasi_data_eda/matchup_correlations_SF001.png",
         height = 4, width = 3)


# number of observations per site and data product
site_counts <- dat %>%
  filter(variable %in% c("eco_sst_v002_clean", "mur_sst",
                         "lst_sst_clean","modis_sst")) %>%
  group_by(point_id, variable, lat, lon) %>%
  summarize(n = sum(!is.na(value)), .groups = "drop")

# print nu
site_coverage <- site_counts %>%
  group_by(variable) %>%
  summarize(n_site_100 = sum(n > 100),
            n_site_50 = sum(n > 50))

print(site_coverage)
write.csv(site_coverage,
          file.path(out_dir, "site_coverage_summary.csv"),
          row.names = FALSE)


site_counts %>%
  ggplot(aes(x = n))+
  geom_histogram()+theme_classic(base_size = 14)+
  facet_wrap(~variable, scales = "free_x")

ggsave(file.path(out_dir, "observations_per_site_histogram.png"),
       height = 4, width = 6)


# ---- maps of observation counts, over the coastline ------------
# scale = "large" (1:10m) needs rnaturalearthhires; falls back to "medium".
land <- tryCatch(
  ne_countries(scale = "large", returnclass = "sf"),
  error = function(e) {
    message("rnaturalearthhires not found; using medium resolution.")
    ne_countries(scale = "medium", returnclass = "sf")
  }
)

bbox <- bbox_from_points(site_counts$lon, site_counts$lat, buffer = 0.10)

# Crop to a padded box so the coastline reaches the plot edges,
# then let coord_sf do the final clipping.
pad  <- 0.25
land <- st_crop(
  st_make_valid(land),
  st_bbox(c(bbox["xmin"] - pad, bbox["ymin"] - pad,
            bbox["xmax"] + pad, bbox["ymax"] + pad),
          crs = st_crs(4326))
)

map_counts <- function(var_name, title, label = FALSE) {
  d <- filter(site_counts, variable == var_name)
  p <- ggplot() +
    geom_sf(data = land, fill = "grey88", colour = "grey55", linewidth = 0.3) +
    geom_point(data = d, aes(x = lon, y = lat, color = n))

  # Sites are tightly clustered, so repel the labels off the points.
  if (label) {
    p <- p + ggrepel::geom_text_repel(
      data = d, aes(x = lon, y = lat, label = point_id),
      size = 2.2, colour = "grey20", segment.colour = "grey60",
      segment.size = 0.2, min.segment.length = 0, max.overlaps = Inf
    )
  }

  p +
    viridis::scale_color_viridis(name = "n obs.") +
    scale_x_continuous(breaks = seq(-180, 180, by = 1)) +
    coord_sf(xlim   = c(bbox["xmin"], bbox["xmax"]),
             ylim   = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    labs(x = NULL, y = NULL, title = title) +
    theme_classic(base_size = 12)
}

map_counts("modis_sst", "MODIS observations per site")

ggsave(file.path(out_dir, "observations_map_modis.png"),
       height = 4, width = 5)


map_counts("lst_sst_clean", "Landsat observations per site")

ggsave(file.path(out_dir, "observations_map_lst.png"),
       height = 4, width = 5)


# Same map with site labels; Landsat covers every site, so this doubles as a
# site-name reference. Needs a larger canvas to fit the labels.
map_counts("lst_sst_clean", "Landsat observations per site", label = TRUE)

ggsave(file.path(out_dir, "observations_map_lst_labelled.png"),
       height = 8, width = 10)


map_counts("modis_sst", "MODIS observations per site", label = TRUE)

ggsave(file.path(out_dir, "observations_map_modis_labelled.png"),
       height = 8, width = 10)
