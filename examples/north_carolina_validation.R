library(tidyr)
library(dplyr)
library(ggplot2)

data <- read.csv("data/raw/tasi_salmon_farms_and_monitoring.csv")
states <- read.csv("examples/models/tasi_modis_sites/tasi_modis_sites_2026-09-05/states.csv")
states_mur <- read.csv("examples/models/tasi_modis_sites/tasi_modis_sites_2026-09-05/states_mur_only.csv")

data %>% 
  filter(stat == "nanmean", variable %in% c("modis_sst_aqua","insitu_sst"),
         !is.na(value),point_id %in% c(sites_1,sites_2,sites_3)) %>%
  group_by(point_id,variable) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = variable,
              values_from = n) %>%
  filter(modis_sst_aqua > 300)
  
  
insitu <- data %>% 
  filter(stat == "nanmean", variable == "insitu_sst",
         point_id %in% unique(states$site)) %>%
  mutate(site = point_id, date = time, obs = value) %>%
  select(site,date,obs)

states %>% filter(variable == "sst") %>%
  ggplot(aes(x = as.Date(date), y = value))+
  geom_line() + theme_classic()+  facet_wrap(~site) + 
  geom_point(data = insitu, mapping = aes( y = obs))

states_mur %>% filter(variable == "mur_sst") %>%
  ggplot(aes(x = as.Date(date), y = value))+
  geom_line() + theme_classic()+  facet_wrap(~site) + 
  geom_point(data = insitu, mapping = aes(y = obs))

states$match_date <- states$date
states_mur$match_date <- states_mur$date

states_insitu  <- merge(states,insitu, by = c("site","date")) 
states_mur_insitu  <- merge(states_mur,insitu, by = c("site","date")) 

ggplot(states_insitu,aes(x = obs,y=value))+
  geom_point()+ theme_classic()+  facet_wrap(~site)+
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2)

ggplot(states_mur_insitu,aes(x = obs,y=value))+
  geom_point()+ theme_classic()+  facet_wrap(~site)+
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2)



data %>%
  filter(stat == "nanmean", variable %in% c("insitu_sst","mur_sst"),
         point_id %in% unique(states$site), stat == "nanmean") %>%
  pivot_wider(names_from = variable) %>% 
  group_by(point_id,time) %>%
  summarize(insitu_sst=mean(insitu_sst, na.rm=T), 
            mur_sst = mean(mur_sst, na.rm=T))%>% 
  ggplot(aes(x = mur_sst-273.15, y = insitu_sst))+
  geom_point()+ facet_wrap(~point_id)+
  theme_classic()+
  geom_abline(aes(slope = 1, intercept = 0), linetype = 2)
  

states_insitu %>% 
  group_by(site) %>%
  summarize(bais = mean(value-obs, na.rm = T),
            MSE = mean( (obs-value)^2, na.rm = T))


states_mur_insitu %>% 
  group_by(site) %>%
  summarize(bais = mean(value-obs, na.rm = T),
            MSE = mean( (obs-value)^2, na.rm = T))





full_model_skill$model <- "Full"
full_model_skill











full_model_skill <- states_insitu %>% 
  group_by(site) %>%
  summarize(MSE = mean((obs-value)^2, na.rm = T))
full_model_skill$model <- "Full"
full_model_skill

mur_model_skill <- states_mur_insitu %>% 
  group_by(site) %>%
  summarize(MSE = mean((obs-value)^2, na.rm = T))
mur_model_skill$model <- "MUR"

mur_only_skill <- data %>% 
  filter(stat == "nanmean", variable %in% c("insitu_sst","mur_sst"),
         point_id %in% unique(states$site)) %>%
  pivot_wider(names_from = "variable", values_from = "value") %>% 
  group_by(point_id)%>%
  group_by(point_id,time) %>%
  summarize(n_insitu = sum(!is.na(insitu_sst)),
            n_mur_sst = sum(!is.na(mur_sst)),
            insitu_sst=mean(insitu_sst, na.rm=T),
            mur_sst = mean(mur_sst, na.rm=T))%>%
  group_by(point_id) %>%
  summarize(MSE = mean((insitu_sst - (mur_sst-273.15))^2, na.rm = T)) %>%
  mutate(site = point_id) %>% select(-point_id)

mur_only_skill$model <- "MUR only"

full_model_skill
mur_model_skill
mur_only_skill

rbind(full_model_skill,mur_model_skill,mur_only_skill) %>%
  filter(site != "ism-secoora-noaa_nos_co_ops_8656") %>%
  ggplot(aes(x = model, y = MSE)) +
  geom_boxplot()+geom_point()+
  theme_classic()


data.frame(site = full_model_skill$site,
           mur_model_diff = mur_model_skill$MSE - full_model_skill$MSE,
           mur_only_diff = mur_only_skill$MSE - full_model_skill$MSE) %>%
  filter(site != "ism-secoora-noaa_nos_co_ops_8656") %>%  
  ggplot(aes(x = site, y = mur_model_diff)) +
  geom_point(color = "orange")+ 
  geom_point(mapping=aes(y = mur_only_diff), color = "purple")+
  theme_classic()



# check efect os data avaiabiltiy on skill differnces
dat_skill_diff <- data.frame(site = full_model_skill$site,
                             mur_model_diff = mur_model_skill$MSE - full_model_skill$MSE,
                             mur_only_diff = mur_only_skill$MSE - full_model_skill$MSE) %>%
  filter(site != "ism-secoora-noaa_nos_co_ops_8656")

head(dat_skill_diff)


MODIS <- data %>% 
  filter(stat == "nanmean", 
         point_id %in% dat_skill_diff$site,
         variable == "modis_sst_aqua") %>%
  group_by(point_id) %>%
  summarize(n = sum(!is.na(value)))


dat_skill_diff$n_modis <- MODIS$n


ECO <- data %>% 
  filter(stat == "nanmean", 
         point_id %in% dat_skill_diff$site,
         variable == "eco_sst_v002_clean") %>%
  group_by(point_id) %>%
  summarize(n = sum(!is.na(value)))
dat_skill_diff$n_eco <- ECO$n


LST <- data %>% 
  filter(stat == "nanmean", 
         point_id %in% dat_skill_diff$site,
         variable == "lst_sst_clean") %>%
  group_by(point_id) %>%
  summarize(n = sum(!is.na(value)))
dat_skill_diff$n_lst <- LST$n


dat_skill_diff


pairs(dat_skill_diff %>%
        select( "mur_model_diff","n_modis","n_eco", "n_lst"))



data$variable %>%
  filter(stat == "nanmean",
         variable %in% c("mur_sst","lst_sst_clean" ),
         point_id %in% dat_skill_diff$site) %>%
  group_by(point_id,time,variable) %>%
  summarize(value = mean(value,na.rm = T) ) %>%
  pivot_wider(names_from = "variable",
              values_from = "value") %>%
  summarize(diff = lst_sst_clean - mur_sst) %>%
  ggplot(aes(x = diff))+
  geom_histogram()+theme_classic()+
  facet_wrap(~point_id)+
  geom_vline(aes(xintercept = -4))



data %>%
  filter(stat == "nanmean",
         variable %in% c("mur_sst","lst_sst_clean"),
         point_id %in% dat_skill_diff$site) %>%
  group_by(point_id,time,variable) %>%
  summarize(value = mean(value,na.rm = T) ) %>%
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  mutate(diff = lst_sst_clean - mur_sst) %>%
  group_by(point_id) %>%
  summarize(n = sum(!is.na(lst_sst_clean)),
           p = sum(diff < -4.0, na.rm = T)/n())


data %>%
  filter(stat == "nanmean",
         variable %in% c("mur_sst","lst_sst_clean"),
         point_id %in% dat_skill_diff$site) %>%
  group_by(point_id,time,variable) %>%
  summarize(value = mean(value,na.rm = T) ) %>%
  pivot_wider(names_from = "variable",
              values_from = "value") %>% 
  mutate(diff = lst_sst_clean - mur_sst) %>%
  group_by(point_id) %>%
  summarize(n = sum(!is.na(lst_sst_clean)),
            p = sum(diff < -4.0, na.rm = T)/n(),
            m = min(diff, na.rm = T))