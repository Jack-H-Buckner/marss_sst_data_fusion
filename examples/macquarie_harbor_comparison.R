

full <- read.csv("examples/models/macquarie_harbor/macquarie_harbor_2026-09-07/states.csv")
mur_only <- read.csv("data/benchmarks/tasi_insitu_and_mur.csv")

mur_only_merge <- mur_only %>% 
  filter(stat == "nanmean",
         variable == "mur_sst",
         point_id %in% unique(full$site)) %>%
  mutate(site = point_id,date=time,se=NaN) %>%
  select(date, site, variable, value, se)
  

dat <- rbind(full,mur_only_merge) %>%
  select(-se) %>%
  pivot_wider(names_from = "variable",
              values_from = "value")

head(dat)


ggplot(dat,aes(x = mur_sst-273.15, y = sst))+
  geom_point()+theme_classic(base_size = 14)+
  facet_wrap(~site)+
  geom_abline(aes(slope = 1, intercept = 0 ))+
  xlab("MUR sst data")+
  ylab("High res data fusion model")
ggsave("examples/figures/macquarie_harbor_coparison.png")


ggplot(dat %>% filter(as.Date(date) > as.Date("2006-10-01")),
       aes(x = as.Date(date), y = sst))+
  geom_line()+geom_line(aes(y = mur_sst - 273.15), color = "grey")+
  theme_classic(base_size = 14)+
  facet_wrap(~site, ncol = 2)+
  geom_abline(aes(slope = 1, intercept = 0 ))+
  xlab("MUR sst data")+
  ylab("High res data fusion model")


ggsave("examples/figures/macquarie_harbor_time_series.png")
  


ggplot(dat,aes(x = mur_sst-273.15, y = sst))+
  geom_point()+theme_classic(base_size = 14)+
  facet_wrap(~site)+
  geom_abline(aes(slope = 1, intercept = 0 ))+
  geom_vline(aes(xintercept = 20), color = "red")+
  geom_hline(aes(yintercept = 20), color = "red")+
  xlab("MUR sst data")+
  ylab("High res data fusion model")




ggsave("examples/figures/macquarie_harbor_time_series.png")


