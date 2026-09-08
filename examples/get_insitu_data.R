library(dplyr)

dat <- read.csv("data/raw/tasi_salmon_20_yrs.csv")
dat_states <- read.csv("examples/models/tasi_temps_low_and_high_res.csv")

head(dat)


insitu <- dat %>% 
  filter(variable  == "insitu_sst", point_id %in% c("CB001","CB002"),
         stat == "nanmean") 
  

insitu$model <- "data"
insitu$run <- NaN
insitu$site <- insitu$point_id
insitu$site_name <- insitu$point_id
insitu$site_type <- "fixed station"
insitu$region <- NaN
insitu$date <- insitu$time
insitu$se <- NaN
insitu$notes <- "in situ measurmens from long term monitoring "

insitu <- insitu %>% 
  select(model,run,site,site_name,site_type,date,se,notes) 
estiamtes <- dat_states %>% 
  select(model,run,site,site_name,site_type,date,se,notes) 

combo <- rbind(insitu,estiamtes)

write.csv(combo,"examples/models/tasi_hi_and_low_res_estimates_an_insitu.csv")



