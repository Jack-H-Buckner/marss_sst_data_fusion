

dat <- read.csv("data/raw/tasi_salmon_20_yrs.csv")


dat_MI <- dat %>%
  filter(variable=="insitu_sst",
         stat == "nanmean",point_id == "CB001")
  
acf_maria_island <- acf(dat_MI$value, na.action = na.pass)


states_mh <- read.csv("examples/models/macquarie_harbor/macquarie_harbor_sites_2026-08-31/states.csv")


ggplot(states_mh %>% filter(variable == "sst"),
       aes(x = as.Date(date), y = value))+
  facet_wrap(~site)+
  geom_line()

states_mh_mur <- read.csv("examples/models/macquarie_harbor/macquarie_harbor_sites_2026-08-31/states_mur_only.csv")

ggplot(states_mh_mur %>% filter(variable == "mur_sst"),
       aes(x = as.Date(date), y = value))+
  facet_wrap(~site)+
  geom_line()


states_combined <- rbind(states_mh,states_mh_mur)


states_combined %>%
  tidyr::pivot_wider(
              id_cols = c("date","site"),
              names_from = "variable",
              values_from = "value") %>%
  ggplot(aes(x = mur_sst, y = sst))+
  facet_wrap(~site)+
  geom_point()+
  geom_abline(aes(slope = 1, intercept = 0))






states_wsb <- read.csv("examples/models/west_storm_bay/west_storm_bay_sites_2026-08-31/states.csv")


ggplot(states_wsb  %>% filter(variable == "sst"),
       aes(x = as.Date(date), y = value))+
  facet_wrap(~site)+
  geom_line()

states_wsb_mur <- read.csv("examples/models/west_storm_bay/west_storm_bay_sites_2026-08-31/states_mur_only.csv")

ggplot(states_wsb_mur %>% filter(variable == "mur_sst"),
       aes(x = as.Date(date), y = value))+
  facet_wrap(~site)+
  geom_line()


states_combined <- rbind(states_wsb ,states_wsb_mur)


states_combined %>%
  tidyr::pivot_wider(
    id_cols = c("date","site"),
    names_from = "variable",
    values_from = "value") %>%
  ggplot(aes(x = mur_sst, y = sst))+
  facet_wrap(~site)+
  geom_point()+
  geom_abline(aes(slope = 1, intercept = 0))


### ---- plot combined states ---- ###
combined_states <- read.csv("examples/models/states_combined.csv")




dat_CB <- combined_states %>% 
  filter(site %in% c("CB001","CB002")) %>%
  tidyr::pivot_wider(id_cols = c("lat","lon","model","run","date",
                                 "site_name","site_type","variable"),
                     names_from = "site",
                     values_from = "value")


combined_states %>% 
  filter(site %in% c("CB001","CB002")) %>% 
  ggplot(aes(x = as.Date(date), y = value, color = site)) + 
  geom_line()+
  facet_wrap(~model)



obs <- readRDS("examples/data/marss_inputs.rds")

keep <- obs$row_var_keys == "insitu_sst"
keep <- keep & (obs$row_site_keys %in% c("CB001","CB002"))
insitu_cols <- obs$ts_matrix[keep,]

plot(insitu_cols[1,],type = "l", col = "blue")
lines(insitu_cols[2,], col = "red")






