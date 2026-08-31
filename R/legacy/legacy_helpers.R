#############################################################
#############################################################
###
### Pre-pipeline helpers, kept for the exploratory scripts.
###
### These are the versions that predate R/00_preprocessing.R.
### They are NOT the ones to use for new work:
###
###   * they expect the raw extract layout, with a `stat`
###     column and the raw `variable` names, rather than the
###     post-`get_variable_value_format()` names the pipeline
###     produces;
###   * `trim_ecostress_outliers()` here takes a `statistic`
###     as its second argument and a single `threshold`, where
###     the pipeline version takes `var` second and a vector
###     of `thresholds`;
###   * `convert_from_K_to_C()` here subtracts 273 rather than
###     273.15, and does no range check.
###
### They live in their own file because they would otherwise
### shadow the working versions in R/00_preprocessing.R for
### anything that sources both, which would silently break
### `run_pipeline()`.
###
### New work should use `trim_ecostress_outliers()` and
### `convert_temperature()` / `convert_from_K_to_C()` from
### R/00_preprocessing.R instead.
###
#############################################################
#############################################################

library(dplyr)
library(ggplot2)
library(lubridate)


trim_ecostress_outliers <- function(data,statistic,var,time_var="time",threshold = 4.0){
  # extract data of interest
  dat_var <- data %>%
    filter(stat==statistic,variable==var)%>%
    filter(!is.na(value))
  # Fit seasonal trend
  dat_var$doy <- yday(as.Date(dat_var[,time_var]))
  lm_mod <- lm(data = dat_var,
               formula = value ~ sin(6.28*doy/365)+cos(6.28*doy/365) + point_id)

  # filter residuals and filter
  dat_var$resid <- scale(lm_mod$residual)

  plt <- ggplot(dat_var,
         aes(x = resid, fill = abs(resid) <threshold))+
    geom_histogram()+
    theme_classic()

  dat_var_filtered <- dat_var %>%
    filter(abs(resid)<threshold)

  # add back in filtered data
  data_out <- data %>% filter( (stat!=statistic) | (variable!=var) )
  print(unique(data_out$variable))
  data_out <- rbind(data_out,dat_var_filtered %>% select(-resid,-doy))
  return(list(data = data_out, data_filtered=dat_var_filtered, plot = plt ))
}


convert_from_K_to_C <- function(data,variables,statistic){
  data_sel <- data %>%
    filter(stat==statistic,variable %in% variables) %>%
    mutate(value = value - 273)
  data_out <- data %>%
    filter((stat!=statistic) | !(variable %in% variables))
  return(rbind(data_out,data_sel))
}
