library(dplyr)
library(reshape2)
library(ggplot2)
library(lubridate)
library(slider)

dat_CS <- read.csv("~/github/tas_sst/data/validation/insitu/tasi_temp_monitoring.csv")
dat_SB <- read.csv("~/github/tas_sst/data/validation/insitu/Storm_Bay_profiles.csv")
dat_MI <- read.csv("~/github/tas_sst/data/validation/insitu/Maria_Island_NRS.csv")
# Maria island -42.2090, 145.1781
# Cape Sorell −42.596, 148.233

dat_MI

lat0 <- mean(dat_SB$lat)
lon0 <- mean(dat_SB$lon)

M_PER_DEG_LAT <- 111320                          # near enough at any latitude
m_per_deg_lon <- 111320 * cos(lat0 * pi / 180)

bin_size_m <- 2000

dat_SB_grouped <- dat_SB %>%
  mutate(
    date = round_date(as.Date(datetime)),
    x_m = (lon - lon0) * m_per_deg_lon,
    y_m = (lat - lat0) * M_PER_DEG_LAT
  ) %>%  
  mutate(
    bin_x  = floor(x_m / bin_size_m),
    bin_y  = floor(y_m / bin_size_m),
    bin_id = paste(bin_x, bin_y, sep = "_"),
    bin_lon = lon0 + (bin_x + 0.5) * bin_size_m / m_per_deg_lon,
    bin_lat = lat0 + (bin_y + 0.5) * bin_size_m / M_PER_DEG_LAT
  ) %>%
  group_by(bin_lat,bin_lon) %>%
  mutate(n = length(unique(date))) %>% 
  filter(n > 150, depth_m < 5.0) %>%
  group_by(bin_lat,bin_lon,date) %>%
  summarize(sst = mean(sst_insitu)) 

ggplot(dat_SB_grouped, aes(x=date, y = sst, 
             color = paste(bin_lat,bin_lon))) + geom_point()

dat_SB <- dat_SB_grouped %>%
  ungroup() %>%
  mutate(lat = bin_lat, lon = bin_lon)%>%
  select(lat,lon,date,sst) %>%
  mutate(site = paste("Storm Bay", lat, lon))


dat_MI_filtered_deep <- as.data.frame(dat_MI[(year(dat_MI$datetime) > 2016) &(dat_MI$depth_m > 10),])
dat_MI_filtered_deep $lat <- -42.596 
dat_MI_filtered_deep $lon <- 148.233
dat_MI_filtered_deep $sst <- dat_MI_filtered_deep $sst_insitu
dat_MI_filtered_deep $date <- as.Date(dat_MI_filtered_deep $datetime)
plot(dat_MI_filtered_deep $date, dat_MI_filtered_deep $sst)

dat_MI_deep  <- dat_MI_filtered_deep  %>%
  ungroup()%>%
  select(lat,lon,date,sst)

dat_MI_deep $site <- "Maria Island deep"
dat_MI_filtered <- as.data.frame(dat_MI[(year(dat_MI$datetime) > 2016) &(dat_MI$depth_m < 10),])
dat_MI_filtered$lat <- -42.596 
dat_MI_filtered$lon <- 148.233 
dat_MI_filtered$sst <- dat_MI_filtered$sst_insitu
dat_MI_filtered$date <- as.Date(dat_MI_filtered$datetime)
plot(dat_MI_filtered$date, dat_MI_filtered$sst)

dat_MI <- dat_MI_filtered %>%
  ungroup()%>%
  select(lat,lon,date,sst)

dat_MI$site <- "Maria Island"
ggplot(rbind(dat_MI,dat_MI_deep), aes(x=date,y=sst,color=site))+
  geom_point() + theme_classic()

maria_island <- rbind(dat_MI,dat_MI_deep)
maria_island$date <- as.Date(maria_island$date)
compare <- reshape2::dcast(maria_island,lat+lon+date~site, value.var = "sst")  %>% 
  arrange(date) %>%
  mutate(
    sst_3dat_surface = slide_index_dbl(`Maria Island`, 
                                       date, mean, na.rm = TRUE,
                                      .before = 1, .after = 1),
    sst_3dat_deep = slide_index_dbl(`Maria Island deep`, 
                                       date, mean, na.rm = TRUE,
                                       .before = 1, .after = 1)
  ) 

mod <- lm(data = compare, formula = sst_3dat_surface ~ sst_3dat_deep)
coefs <- coef(mod)
dat_MI_deep <- dat_MI_deep %>%
  mutate(sst = coefs[1] + coefs[2]*sst)

dat_CS <- dat_CS %>% 
  filter(SITE_NAME == "Cape Sorell") %>%
  mutate(date = round_date(as.Date(DATETIME_LOCAL))) %>%
  group_by(date) %>%
  summarize(sst = mean(TEMP_C))

dat_CS$lat <- -42.2090 
dat_CS$lon <- 145.1781
dat_CS$site <- "Cape Sorell"

dat <- rbind(dat_SB,dat_MI,dat_MI_deep,dat_CS)

ggplot(dat, aes(x = lon, y = lat))+
  geom_point()

dat$depth_m <- 1.0
dat <- dat %>%
  group_by(site) %>%
  mutate(lat = mean(lat), lon = mean(lon), .groups = "drop") %>%
  select(lat,lon,date,site,depth_m,sst)

write.csv(dat, "~/github/tas_sst/data/validation/insitu/combined_with_deep.csv")




d <- dat %>% filter(site =="Maria Island" )
d$lat[1]
d$lon[1]