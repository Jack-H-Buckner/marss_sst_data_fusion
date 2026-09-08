

dat_tasi <- read.csv("data/raw/tasi_salmon_farms_and_monitoring.csv")



dat_tasi %>%
  filter(variable %in% c("insitu_sst", "mur_sst"),
         stat == "nanmean", !is.na(value), ) %>%
  mutate(storm_bay_test_site = point_id %in% c("George III Reef", "Mouldies Hole", 
                                     "Whale Head", "Iron Pot", "One Tree Point"),
         storm_bay_calibration = point_id == "One Tree Point",
         east_coast_test_site = point_id %in% c("Wedge Island", "crescent bay", "Tasman Island",
                                                "Deep Glen Bay", "Cape Peron", "Magistrates Point"),
         east_coast_calibration = point_id == "Tasman Island",
         macquarie_test_site = point_id %in% c("Cape Sorell","Table Head_1 215_1","North East Double Cove 220",
                                               "Liberty Point Central Harbour 214", "North East of Hogan Cove 213" ),
         east_coast_calibration = point_id == "Cape Sorell",
         tasi_with_modis = point_id %in% c("Wedge Island", "Magistrates Point", "Iron Pot","George III Reef", "Cape Peron"),
         modis_calibration = point_id == "Iron Pot",
         included_time = as.Date(time) > as.Date("2018-10-01")) %>%
  select(-radius_m) %>%
  write.csv("data/benchmarks/tasi_insitu_and_mur.csv")
  


dat_tasi %>%
  filter(variable %in% c("insitu_sst", "mur_sst"),
         stat == "nanmean", !is.na(value), 
         as.Date(time) > as.Date("2018-10-01"),
         substr(time,9,10) %in% c("01","10","20")) %>%
  mutate(storm_bay_test_site = point_id %in% c("George III Reef", "Mouldies Hole", 
                                               "Whale Head", "Iron Pot", "One Tree Point"),
         storm_bay_calibration = point_id == "One Tree Point",
         east_coast_test_site = point_id %in% c("Wedge Island", "crescent bay", "Tasman Island",
                                                "Deep Glen Bay", "Cape Peron", "Magistrates Point"),
         east_coast_calibration = point_id == "Tasman Island",
         macquarie_test_site = point_id %in% c("Cape Sorell","Table Head_1 215_1","North East Double Cove 220",
                                               "Liberty Point Central Harbour 214", "North East of Hogan Cove 213" ),
         east_coast_calibration = point_id == "Cape Sorell",
         tasi_with_modis = point_id %in% c("Wedge Island", "Magistrates Point", "Iron Pot","George III Reef", "Cape Peron"),
         modis_calibration = point_id == "Iron Pot") %>%
  select(-radius_m) %>%
  write.csv("data/benchmarks/tasi_insitu_and_mur_small.csv")
  

dat_nc <- read.csv("data/raw/coastal_sst_data_nc_clouds.csv")

dat_nc %>%
  filter(variable %in% c("insitu_sst", "mur_sst"),
         stat == "nanmean", !is.na(value)) %>%
  mutate(test_site = point_id %in% c("ism-secoora-noaa_nos_co_ops_8656",
                                     "neuse-river-at-marker-7-modmo", 
                                     "neuse-river-at-marker-15-modm", 
                                     "neuse-river-at-marker-17-modm",
                                     "neuse-river-at-marker-9-modmo",
                                     "neuse-river-at-cm-22-fairfiel-2",
                                     "neuse-river-at-marker-38-modm",
                                     "neuse-river-at-marker-52-a-mo"),
         calibration = point_id == "ism-secoora-noaa_nos_co_ops_8656",
         included_time = T) %>%
  select(-radius_m) %>%
  write.csv("data/benchmarks/nc_insitu_and_mur.csv")


min(dat_nc$time )
max(dat_nc$time )
