

dat <- read.csv("data/raw/tasi_salmon_20_yrs.csv")



source("R/00_preprocessing.R")
data <- dat
  


data %>% 
  filter(variable %in% "modis_sst_aqua") %>%
  group_by(point_id) %>%
  summarize(n = sum(!is.na(value)), .groups = "drop")  %>%
  summarise(sum(n > 1000))


source("R/00_preprocessing.R")
config <- "parameters/preprocessing.R"
df <- run_pipeline(data,config)
df




value_var <- "value"
vars <- list(eco_sst_v002_clean = list(variable = "eco_sst_v002_clean", stat = "nanmean", radius_m = 300),
             eco_wind_speed_era5 = list(variable = "eco_wind_speed_era5", stat = "nanmean", radius_m = 300),
             eco_hour_v002 = list(variable = "eco_hour_v002", stat = "nearest", radius_m = 0),
             lst_sst_clean = list(variable = "lst_sst_clean", stat = "nanmean", radius_m = 300),
             lst_wind_speed_era5 = list(variable = "lst_wind_speed_era5", stat = "nanmean", radius_m = 300),
             lst_hour = list(variable = "lst_hour", stat = "nearest", radius_m = 0),
             modis_sst = list(variable = "modis_sst_aqua", stat = "nanmean", radius_m = 1000),
             mur_sst = list(variable = "mur_sst", stat = "nanmean", radius_m = 1000),
             insitu_sst = list(variable = "insitu_sst", stat = "nanmean", radius_m = 300))



dat_formatted <- get_variable_value_format(data,vars,value_var)
data_units <- convert_from_K_to_C(dat_formatted,
                                  c("eco_sst_v002_clean","lst_sst_clean",
                                    "modis_sst","mur_sst"))
source("R/00_preprocessing.R")
data_phys <- trim_implausibe_values(data_units ,"eco_sst_v002_clean")

outliers <- trim_ecostress_outliers(data_phys$data,"eco_sst_v002_clean",time_var="time",thresholds = 2.58)

outliers_wind_eco <- wind_filter(outliers$data,
                                 list(target = "eco_sst_v002_clean", 
                                      wind = "eco_wind_speed_era5",
                                      time = "eco_hour_v002"),
                                 threshold = 2.0,
                                 early = 12, 
                                 late = 18)

out <- match_up_filter(outliers_wind_eco$data,
                vars = list(target = "eco_sst_v002_clean", reference = "modis_sst"),
                lower   = 2.58,
                upper   = 3.0,)

out <- match_up_filter(out$data,
                       vars = list(target = "eco_sst_v002_clean", reference = "modis_sst"),
                       lower   = 2.58,
                       upper   = 3.0,
)

out <- match_up_filter(out$data,
                       vars = list(target = "lst_sst_clean", reference = "modis_sst"),
                       lower   = 2.58,
                       upper   = 4.0,
)
out <- match_up_filter(out$data,
                       vars = list(target = "lst_sst_clean", reference = "modis_sst"),
                       lower   = 2.58,
                       upper   = 4.0,
)


source("R/00_preprocessing.R")
variance <- match_up_variance(out$data,
                  targets = c("lst_sst_clean","eco_sst_v002_clean"),
                  references = c("insitu_sst","modis_sst"),
                  pairs = "cross")
                              

