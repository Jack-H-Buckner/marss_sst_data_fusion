library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(lubridate)

dat <- read.csv("data/raw/tasi_salmon_20_yrs.csv")
dat$time <- as.Date(dat$time)


dat %>% 
  filter(stat %in% c("nanmean", "nearest")) %>%
  group_by(point_id,variable) %>%
  summarize(n = sum(!(is.na(value))), .groups = "drop") %>%
  ggplot(aes(x = n))+
  facet_wrap(~variable, scales = "free_x")+
  geom_histogram()




dat %>% 
  filter(point_id == "CB001",
        stat %in% c("nanmean", "nearest"),
         variable %in% c("lst_sst_clean","eco_sst_v002_clean",
                         "modis_sst_aqua", "mur_sst")) %>%
  ggplot(aes(x = time, y = value-273.15, color = point_id))+
    facet_wrap(~variable, scales = "free")+
    geom_point() + 
    theme_classic(base_size = 12)+
    theme(legend.position = "none")








dat %>% 
  filter(stat %in% c("nanmean", "nearest"),
         variable %in% c("eco_sst_v002_clean")) %>%
  ggplot(aes(x = time, y = value, color = point_id))+
  facet_wrap(~variable, scales = "free")+
  geom_point() + 
  theme_classic(base_size = 12)+
  theme(legend.position = "none")

source("R/01_format_data_for_marss.R")
source("R/legacy_helpers.R")   # pre-pipeline trim_ecostress_outliers / convert_from_K_to_C
dat <- read.csv("data/raw/salmon_0819_run_points.csv")
dat$time <- as.Date(dat$time)
dat_trimmed <- trim_ecostress_outliers(dat,"nanmean","eco_sst_v002_clean",time_var="time",threshold = 10.0)
dat_trimmed$plot
dat_trimmed <- trim_ecostress_outliers(dat_trimmed$data,"nanmean","eco_sst_v002_clean",time_var="time",threshold = 3.0)
dat_trimmed$plot



dat_trimmed <- dat_trimmed$data

dat_trimmed %>% 
  filter(stat == "nanmean", variable == "eco_sst_v002_clean") %>%
  ggplot(aes(x = time, y = value, color = point_id))+
  facet_wrap(~variable, scales = "free")+
  geom_point() + 
  theme_classic(base_size = 12)+
  theme(legend.position = "none")

dat %>% 
  filter(stat %in% c("nanmean", "nearest"),
         point_id == "CB001",
         variable %in% c("insitu_sst","mur_sst")) %>%
  dcast(time~variable ) %>% mutate(diff = mur_sst-insitu_sst)%>%
  ggplot(aes(x=diff))+
  geom_histogram()



  ggplot(aes(x = time, y = value, color = point_id))+
  facet_wrap(~variable, scales = "free")+
  geom_point() + 
  theme_classic(base_size = 12)+
  theme(legend.position = "none")


  
  
  
### in situ match ups 
dat <- read.csv("data/raw/tasi_salmon_sst_data.csv")


dat_trimmed <- trim_ecostress_outliers(dat,"nanmean","eco_sst_v002_clean",time_var="time",threshold = 10.0)
dat_trimmed <- trim_ecostress_outliers(dat_trimmed$data,"nanmean","eco_sst_v002_clean",time_var="time",threshold = 3.0)
dat_trimmed$plot
dat_trimmed <- trim_ecostress_outliers(dat_trimmed$data,"nanmean","eco_sst_v002_clean",time_var="time",threshold = 2.58)
dat_trimmed$plot


dat <- dat_trimmed$data

dat %>% 
  filter(variable == "insitu_sst") %>%
  group_by(point_id) %>%
  summarize(n = sum(!is.na(value)))


dat_high_skin_temp <- dat %>% 
  filter(stat %in% c("nearest","nanmean")) %>%
  select(-radius_m,-stat) %>% 
  reshape2::dcast(point_id+lat+lon+aoi+time~variable, 
                  fun.aggregate = mean) %>% 
  mutate(local_time_eco = (eco_hour_v002 + lon/15) %% 24,
         high_zeneth = (18 > local_time_eco) & (local_time_eco > 10),
         eco_sst_valid = if_else(eco_wind_speed_era5 > 2.0, eco_sst_v002_clean, NA,   missing = NA),
         lst_sst_valid = if_else((lst_wind_speed_era5 > 2.0), lst_sst_clean, NA,  missing = NA)) 
#        reshape2::melt(id.vars = c("point_id","lat","lon","aoi","time")) 
head(dat_high_skin_temp)





dat_match_ups <- dat_high_skin_temp %>% 
  filter((!is.na(modis_sst_aqua))& (!is.na(lst_sst_valid)))

sd_lst <- sd(dat_match_ups$lst_sst_valid - dat_match_ups$modis_sst_aqua)

dat_match_ups_filter_1 <- dat_match_ups %>% 
  filter((lst_sst_valid - modis_sst_aqua) > -2.58*sd_lst)


dat_match_ups_filter_1 %>%
  ggplot(aes(x = modis_sst_aqua, y = lst_sst_valid)) +
  geom_point()+theme_classic()+
  geom_abline(aes(intercept = 0, slope = 1)) 


sd_lst <- sd(dat_match_ups_filter_1$lst_sst_valid - dat_match_ups_filter_1$modis_sst_aqua)

dat_match_ups_filter_2 <- dat_match_ups_filter_1 %>% 
  filter((lst_sst_valid - modis_sst_aqua) > -2.58*sd_lst)


dat_match_ups_filter_2 %>%
  ggplot(aes(x = modis_sst_aqua, y = lst_sst_valid)) +
  geom_point()+theme_classic()+
  geom_abline(aes(intercept = 0, slope = 1)) 


dat_match_ups_filter_2 <- dat_match_ups_filter_2 %>% 
  filter((!is.na(insitu_sst))& (!is.na(lst_sst_valid)))


hist(dat_match_ups_filter_2$lst_sst_valid- dat_match_ups_filter_2$insitu_sst - 273.15)

sd(dat_match_ups_filter_2$lst_sst_valid - dat_match_ups_filter_2$insitu_sst, na.rm = T)^2

plot(dat_match_ups_filter_2$insitu_sst, dat_match_ups_filter_2$lst_sst_valid)




#### ecostress cold outliers

dat_match_ups <- dat_high_skin_temp %>% 
  filter((!is.na(modis_sst_aqua))& (!is.na(eco_sst_valid)))


mod <- lm(eco_sst_valid ~ modis_sst_aqua, dat_high_skin_temp)
sd_eco <- sd(mod$residuals)

dat_match_ups_filter_1 <- dat_high_skin_temp %>% 
  filter((eco_sst_valid - coef(mod)[1] - coef(mod)[2]*modis_sst_aqua) > -2.58*sd_eco,
         (eco_sst_valid - coef(mod)[1] - coef(mod)[2]*modis_sst_aqua) < 3.0*sd_eco)


dat_match_ups_filter_1 %>%
  ggplot(aes(x = modis_sst_aqua, y = lst_sst_valid)) +
  geom_point()+theme_classic()+
  geom_abline(aes(intercept = 0, slope = 1)) 


mod <- lm(eco_sst_valid ~ modis_sst_aqua, dat_match_ups_filter_1)
sd_eco <- sd(mod$residuals)

dat_match_ups_filter_2 <- dat_match_ups_filter_1 %>% 
  filter((eco_sst_valid - coef(mod)[1] - coef(mod)[2]*modis_sst_aqua) > -2.58*sd_eco,
         (eco_sst_valid - coef(mod)[1] - coef(mod)[2]*modis_sst_aqua) < 3.0*sd_eco)

dat_match_ups_filter_2 %>%
  ggplot(aes(x = modis_sst_aqua, y = eco_sst_valid)) +
  geom_point()+theme_classic()+
  geom_abline(aes(intercept = 0, slope = 1)) 


dat_match_ups_filter_2 <- dat_match_ups_filter_2 %>% 
  filter((!is.na(insitu_sst))& (!is.na(eco_sst_valid)))


hist(dat_match_ups_filter_2$eco_sst_valid- dat_match_ups_filter_2$insitu_sst - 273.15)

sd(dat_match_ups_filter_2$eco_sst_valid - dat_match_ups_filter_2$insitu_sst, na.rm = T)^2

plot(dat_match_ups_filter_2$insitu_sst, dat_match_ups_filter_2$eco_sst_valid)


### test out filtering functions
dat <- read.csv("data/raw/tasi_salmon_20_yrs.csv")

data <- dat
value_var <- "value"
vars <- list(eco_sst_v002_clean = list(variable = "eco_sst_v002_clean", stat = "nanmean"),
             eco_wind_speed_era5 = list(variable = "eco_wind_speed_era5", stat = "nanmean"),
             eco_hour_v002 = list(variable = "eco_hour_v002", stat = "nearest"))
             
keep <- rep(F,nrow(data))
var_name <- rep("",nrow(data))
nms <- names(vars)
i <- 0
for(var in vars){
  i <-  i+ 1
  keep_var <- rep(T,nrow(data))
  for(nm in names(var)){
    keep_var <- keep_var & (data[nm]==var[nm])
  } 
  var_name[keep_var] <- nms[i]
  keep <- keep | keep_var
  
}
data_vars <- data[keep,]
data_vars <- data_vars[,!(names(data_vars)%in%names(vars[[1]]))]
data_vars$variable <- var_name[keep]
data_vars$value <- data_vars[,value_var]
if(value_var != "value"){
  data_vars <- data_vars[,!value_var]
}

head(data_vars)


#' Collapse variable/stat column pairs into a single variable column
#'
#' Selects rows of `data` matching each specification in `vars`, labels them with
#' the corresponding name from `vars`, drops the columns used for matching, and
#' renames `value_var` to `value`.
#'
#' @param data A long format data frame.
#' @param vars A named list of named lists. Each inner list gives the column
#'   values identifying one variable, e.g.
#'   `list(var_1 = list(variable = "lst_sst", stat = "nanmean"),
#'         var_2 = list(variable = "lst_wind", stat = "nearest"))`.
#' @param value_var Character scalar naming the column holding the values.
#'
#' @return A data frame with `variable` and `value` columns replacing the
#'   matching columns and `value_var`.
#' @export
get_variable_value_format <- function(data, vars, value_var) {
  stopifnot(
    is.data.frame(data),
    is.list(vars), length(vars) > 0,
    !is.null(names(vars)), all(nzchar(names(vars))),
    is.character(value_var), length(value_var) == 1
  )
  
  key_cols <- unique(unlist(lapply(vars, names)))
  missing_cols <- setdiff(c(key_cols, value_var), names(data))
  if (length(missing_cols)) {
    stop("Columns not found in `data`: ", paste(missing_cols, collapse = ", "))
  }
  
  keep     <- rep(FALSE, nrow(data))
  var_name <- rep(NA_character_, nrow(data))
  
  for (i in seq_along(vars)) {
    spec <- vars[[i]]
    keep_var <- rep(TRUE, nrow(data))
    for (nm in names(spec)) {
      keep_var <- keep_var & (data[[nm]] %in% spec[[nm]])
    }
    if (any(keep & keep_var)) {
      warning("Rows match more than one entry of `vars`; later entries win.")
    }
    var_name[keep_var] <- names(vars)[i]
    keep <- keep | keep_var
  }
  
  if (!any(keep)) warning("No rows matched any entry of `vars`.")
  
  data_vars <- data[keep, setdiff(names(data), key_cols), drop = FALSE]
  data_vars$variable <- var_name[keep]
  data_vars$value    <- data[[value_var]][keep]
  if (!identical(value_var, "value")) data_vars[[value_var]] <- NULL
  
  rownames(data_vars) <- NULL
  data_vars
}
get_variable_value_format(data,vars,value_var) %>% head()



data <- data %>%
  filter(variable %in% c(vars$target,vars$wind,vars$time),
         ) %>%
  filter((stat==stats$target)|(variable!=vars$target),
         (stat==stats$wind)|(variable!=vars$wind),
         (stat==stats$time)|(variable!=vars$time))%>%
  select(-radius_m,-stat)


data <- dat
vars <- list(target = "eco_sst_v002_clean", wind = "eco_wind_speed_era5", time = "eco_hour_v002")
stats <- list(target = "nanmean", wind = "nanmean", time = "nearest")
data <- data %>%
  filter(variable %in% c(vars$target,vars$wind,vars$time)) %>%
  filter((stat==stats$target)|(variable!=vars$target),
         (stat==stats$wind)|(variable!=vars$wind),
         (stat==stats$time)|(variable!=vars$time))%>%
  select(-radius_m,-stat)

unique(data$variable)

early = 10
late = 18
threshold = 2.0
local_time = F 
  
dat_var <- data %>% 
  filter(variable %in% c(vars$target,vars$wind,vars$time)) %>%
  tidyr::pivot_wider(names_from = "variable", values_from = c("value")) 

# calcualte local solar time 
if(!local_time){
  solar_time <- (dat_var[vars$time] + dat_var["lon"]/15) %%24
} else{
  solar_time <- dat_var[vars$time] 
}

# calculate time window for the filter
valid_time <- (solar_time < late) & (solar_time > early)

# get valid wind
valid_wind_speed <- (dat_var[vars$wind] > threshold)

# combined valid observations 
valid <- (valid_time & valid_wind_speed)
valid <- if_else(as.vector(valid),T,F,missing=T)  
dat_var_long <- dat_var[valid,] %>% 
  pivot_longer(cols = !c(point_id,lat,lon,aoi,time), 
    names_to = "variable", values_to = "value") 


# add filtered variable back to full data set 
dat_filtered <- data %>% filter(!(variable %in% unique(dat_var_long$variable)))

dat_filtered <- rbind(dat_filtered,dat_var_long)


dat_filtered %>% head()


