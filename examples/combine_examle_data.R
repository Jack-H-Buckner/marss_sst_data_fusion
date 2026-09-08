

dat_tasi_modis <- read.csv("examples/models/tasi_modis_sites/tasi_modis_sites_2026-09-05/states.csv")
dat_tasi_modis$sites <- "Tasi with MODIS"
dat_tasi_modis$converged <- "yes"
dat_east_coast <- read.csv("examples/models/east_coast/east_coast_2026-09-05/states.csv")
dat_east_coast$sites <- "Tasi east cost"
dat_east_coast$converged <- "yes"
dat_storm_bay <- read.csv("examples/models/storm_bay/storm_bay_models_2026-09-05/states.csv")
dat_storm_bay$sites <- "Tasi storm bay"
dat_storm_bay$converged <- "no"
dat_macquarie_harbor <- read.csv("examples/models/macquarie_harbor/macquarie_harbor_2026-09-07/states.csv")
dat_macquarie_harbor$sites <- "Macquarie Harbor"
dat_macquarie_harbor$converged <- "no"
dat_nc <- read.csv("examples/models/north_carolina_clouds/north_carolina_2026-09-04/states.csv")
dat_nc$sites <- "Pamlico Sound North Carolina"
dat_nc$converged <- "yes"
dat_nc_imputed <- read.csv("examples/models/north_carolina_imputed/north_carolina_2026-09-04/states.csv")
dat_nc_imputed$sites <- "Pamlico Sound North Carolina (imputed MODIS)"
dat_nc_imputed$converged <- "yes"

dat_states <- rbind(dat_tasi_modis,dat_east_coast,dat_storm_bay,
             dat_macquarie_harbor,dat_nc,dat_nc_imputed)


write.csv(dat_states,"examples/models/estiamtes_tasi_nc.csv")



dat_mur_nc <- read.csv("data/benchmarks/nc_insitu_and_mur.csv") %>%
  select(point_id, lat, lon, time, variable, value)

dat_tasi_nc <- read.csv("data/benchmarks/tasi_insitu_and_mur.csv")%>%
  select(point_id, lat, lon, time, variable, value)


dat_bench <- rbind(dat_tasi_nc,dat_mur_nc)

names(dat_bench) <- c("site","lat","lon","date","variable","value")
head(dat_bench)


dat_bench %>% filter(site %in% unique(dat_states$site)) %>%
  write.csv("examples/models/benchmark_data_tasi_nc.csv")


