library(MARSS)
library(dplyr)
library(reshape2)
library(ggplot2)



get_states_mur_factors <- function(model,scales,dates,design, y_matrix){
  
  # extract model paramters 
  params <- coef(model,type="matrix")
  
  ### gets site means and harmonics 
  # define parameters name lists
  site_mean = paste0("mu_", seq_len(design$n_site))
  sin_coef = paste0("c1_", seq_len(design$n_site))
  cos_coef = paste0("c2_", seq_len(design$n_site))
  # extract paramters by names 

  mu <- model$par$A[site_mean,]
  c1 <- model$par$A[sin_coef,]
  c2 <- model$par$A[cos_coef,]
  
  # define parameters name lists
  site_mean_mur = paste0("mu_mur_", seq_len(design$n_site))
  sin_coef_mur = paste0("c1_mur_", seq_len(design$n_site))
  cos_coef_mur = paste0("c2_mur_", seq_len(design$n_site))
  # extract paramters by names 
  mu_mur <- model$par$A[site_mean_mur,]
  c1_mur <- model$par$A[sin_coef_mur,]
  c2_mur <- model$par$A[cos_coef_mur,]
  
  # define harmonics 
  dt <- matrix(model$model$fixed$d,nrow = 2)
  days <- seq(dates$start,dates$end,by=1)
  sin_doy <- dt[1,]
  cos_doy <- dt[2,]
  
  # get design variables
  n <- design$n_site
  factors_inds <- y_matrix$row_var_keys == "mur_sst" 
  site_names <- y_matrix$row_site_keys[factors_inds]
  X <- model$states[1:n,]
  
  print(scales$sigma)
  ### reshape sst data
  print("sst")
  sst <- scales$mu + scales$sigma * (X + mu +  + 
                                              matrix(c1,ncol = 1, nrow = n) %*% matrix(sin_doy, nrow = 1) +
                                              matrix(c2,ncol = 1, nrow = n) %*% matrix(cos_doy, nrow = 1)) 
  
  print(sst)
  rownames(sst) <- site_names
  colnames(sst) <- as.character(days)
  dat_sst <- reshape2::melt(sst)  
  names(dat_sst) <- c("site","date","value")
  dat_sst$variable <- "sst"
  
  
  ### Calculate mur predictions
  print("sst mur")
  sst_mur <- scales$mu + scales$sigma * (X + mu_mur +   
                                           matrix(c1_mur,ncol = 1, nrow = n) %*% matrix(sin_doy, nrow = 1) +
                                           matrix(c2_mur,ncol = 1, nrow = n) %*% matrix(cos_doy, nrow = 1)) 
  
  rownames(sst_mur) <- site_names
  colnames(sst_mur) <- as.character(days)
  dat_mur <- reshape2::melt(sst_mur)  
  names(dat_mur) <- c("site","date","value")
  dat_mur$variable <- "mur_sst"
  
  
  ### harmonics
  print("sst mean")
  sst_mean <- scales$mu + scales$sigma * (mu +
                                            matrix(c1,ncol = 1, nrow = n) %*% matrix(sin_doy, nrow = 1) +
                                            matrix(c2,ncol = 1, nrow = n) %*% matrix(cos_doy, nrow = 1))
  
  rownames(sst_mean) <- site_names
  colnames(sst_mean) <- as.character(days)
  sst_mean <- as.data.frame(reshape2::melt(sst_mean)  )
  names(sst_mean) <- c("site","date","value")
  sst_mean$variable <- "mean_sst"
  
  print("sst mur mean")
  sst_mean_mur <- scales$mu + scales$sigma * (mu_mur +
                                            matrix(c1_mur,ncol = 1, nrow = n) %*% matrix(sin_doy, nrow = 1) +
                                            matrix(c2_mur,ncol = 1, nrow = n) %*% matrix(cos_doy, nrow = 1))
  
  rownames(sst_mean_mur) <- site_names
  colnames(sst_mean_mur) <- as.character(days)
  sst_mean_mur <- as.data.frame(reshape2::melt(sst_mean_mur)  )
  names(sst_mean_mur) <- c("site","date","value")
  sst_mean_mur$variable <- "mean_sst_mur"
  
  
  out <- rbind(dat_sst, dat_mur, sst_mean, sst_mean_mur)
  return(out)
}

model <- readRDS("models/marss_test_mur_fe/fit_chunk_30.rds")
y_matrix <- readRDS("models/marss_test_mur_fe/observations.rds")
scales <- readRDS("models/marss_test_mur_fe/scales.rds")
dates <- readRDS("models/marss_test_mur_fe/dates.rds")
design <- readRDS("models/marss_test_mur_fe/design.rds")
states <- get_states_mur_factors(model,scales,dates,design,y_matrix)


ggplot(states,
       aes(x = as.Date(date), y = value, color = variable))+
  geom_line()+facet_wrap(~site)+
  theme_classic()



ggplot(states,
       aes(x = as.Date(date), y = value, color = site))+
  geom_line()+facet_wrap(~variable)+
  theme_classic()


states %>%
  filter(variable %in% c("sst","mur_sst")) %>%
  reshape2::dcast(date+site ~ variable ) %>%
  mutate(diff = sst- mur_sst) %>%
  ggplot(aes(x = as.Date(date), y = diff))+
  geom_line()+facet_wrap(~site)+
  theme_classic()



