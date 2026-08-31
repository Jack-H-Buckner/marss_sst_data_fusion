library(MARSS)
library(dplyr)
library(reshape2)
library(ggplot2)



get_states_mur_factors <- function(model,scales,dates,design, y_matrix){
  
  # extract model parameters 
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
  
  
  ### reshape mur data
  X_seasonal <- scales$mu + scales$sigma * (X + mu +  + 
                                  matrix(c1,ncol = 1, nrow = n) %*% matrix(sin_doy, nrow = 1) +
                                  matrix(c2,ncol = 1, nrow = n) %*% matrix(cos_doy, nrow = 1)) 
  
  rownames(X_seasonal) <- site_names
  colnames(X_seasonal) <- as.character(days)
  dat_sst <- reshape2::melt(X_seasonal)  
  names(dat_sst) <- c("site","date","value")
  dat_sst$variable <- "sst"

  
  ### Calculate mur predictions
  sst_mur <- (params$Z %*% model$states)[factors_inds,]
  print(factors_inds)
  print(nrow(sst_mur))
  print(ncol(sst_mur))
  sst_mur <- scales$mu + scales$sigma * (sst_mur + mu +   
                                  matrix(c1,ncol = 1, nrow = n) %*% matrix(sin_doy, nrow = 1) +
                                  matrix(c2,ncol = 1, nrow = n) %*% matrix(cos_doy, nrow = 1)) 
  
  rownames(sst_mur) <- site_names
  colnames(sst_mur) <- as.character(days)
  dat_mur <- reshape2::melt(sst_mur)  
  names(dat_mur) <- c("site","date","value")
  dat_mur$variable <- "mur_sst"

  
  ### harmonics
  sst_mean <- scales$mu + scales$sigma * (mu +
    matrix(c1,ncol = 1, nrow = n) %*% matrix(sin_doy, nrow = 1) +
                                         matrix(c2,ncol = 1, nrow = n) %*% matrix(cos_doy, nrow = 1))

  rownames(sst_mean) <- site_names
  colnames(sst_mean) <- as.character(days)
  sst_mean <- as.data.frame(reshape2::melt(sst_mean)  )
  names(sst_mean) <- c("site","date","value")
  sst_mean$variable <- "mean_sst"

  out <- rbind(dat_sst, dat_mur, sst_mean)
  return(out)
}

model <- readRDS("models/marss_test_sites/fit_chunk_27.rds")
y_matrix <- readRDS("models/marss_test_sites/observations.rds")
scales <- readRDS("models/marss_test_sites/scales.rds")
dates <- readRDS("models/marss_test_sites/dates.rds")
design <- readRDS("models/marss_test_sites/design.rds")
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
  mutate(diff = sst - mur_sst) %>%
  ggplot(aes(x = as.Date(date), y = diff))+
  geom_line()+facet_wrap(~site)+
  theme_classic()


