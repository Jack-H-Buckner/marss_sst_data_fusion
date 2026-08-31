library(ggplot2)
library(dplyr)
library(reshape2)
library(lubridate)

source("R/01_format_data_for_marss.R")
source("R/marss_matrix_functions.R")
sites <- c("CB001","SF058","SF017","SF001","SF051")
dat <- read.csv("data/raw/salmon_small_eco_points_9.csv") %>%
  filter(point_id %in% sites)
#filter(point_id != "CB002")

ggplot(dat %>% group_by(point_id) %>%
         summarize(lat = mean(lat), lon = mean(lon)), 
       aes(y = lat, x = lon, label = point_id)) +
  geom_text()

obs_nms <- c("insitu_sst","lst_sst_clean",
             "eco_sst_v002_clean","modis_sst_aqua",
             "mur_sst")

obs_covar_nms <- c("eco_wind_speed_era5","eco_hour_v002",
                   "lst_hour","lst_wind_speed_era5")


####################################################
####################################################
###
### Format obserations for MARSS
###
####################################################
####################################################
# build observations ts matrix
obs_mat_y <- convert_long_format_marss(dat, obs_nms, 1,value.var = "value",
                                       variable.var = "variable",site.var = "point_id",
                                       time.var = "time", statistic.var = "stat", statistic = c("nanmean"))

y <- obs_mat_y$ts_matrix
y_scale <- t(scale(t(y)))

####################################################
####################################################
###
### Define observation model
###
####################################################

insitu_rows <- which(obs_mat_y$row_var_keys == "insitu_sst")
lst_rows <- which(obs_mat_y$row_var_keys == "lst_sst_clean")
eco_rows <- which(obs_mat_y$row_var_keys == "eco_sst_v002_clean")
modis_rows <- which(obs_mat_y$row_var_keys == "modis_sst_aqua")
mur_rows <- which(obs_mat_y$row_var_keys == "mur_sst" )

####################################################
### -------------------------------------------- ###
####################################################
### Define observation matrix 
####################################################
k_obs <- nrow(obs_mat_y$ts_matrix)
n_site <- length(unique(obs_mat_y$row_site_index))

# site level component
I_sites <- matrix(0, nrow=n_site,ncol = n_site)
diag(I_sites) <- 1
chi_obs <- I_sites[obs_mat_y$row_site_index,]

# full obs matrix 
Z <- chi_obs

####################################################
### Define observation constants 
####################################################
a <- rep(0,k_obs)
# in situ
mu_site <- paste0("mu_", seq_len(n_site))
insitu_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="insitu_sst" ]
a[insitu_rows] <- mu_site[insitu_inds]

# landsat
mu_lst <- paste0("mu_lst+mu_", seq_len(n_site))
lst_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="lst_sst_clean" ]
a[lst_rows] <- mu_lst[lst_inds]

# ecostress
mu_eco <- paste0("mu_eco+mu_", seq_len(n_site))
eco_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="eco_sst_v002_clean" ]
a[eco_rows] <- mu_eco[eco_inds]

# modis - independent intercepts for the low res data sources
mu_modis <- paste0("mu_modis+mu_", seq_len(n_site))
modis_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="modis_sst_aqua" ]
a[modis_rows] <- mu_modis[modis_inds]

# mur
mu_mur <- paste0("mu_mur_", seq_len(n_site))
mur_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="mur_sst" ]
a[mur_rows] <- mu_mur[mur_inds]


####################################################
### Define observation covariates
####################################################
### build covariates ts matrix
dat_single <- dat %>% filter(point_id == "SF051")
obs_covar_mat_d <- convert_long_format_marss(dat_single,
                                             obs_covar_nms, 1,value.var = "value",
                                             variable.var = "variable",site.var = "point_id",
                                             time.var = "time", statistic.var = "stat", 
                                             statistic = c("nearest", "nanmean"))

# calculate harmonics and functional forms for skin effects 
sin_doy <- sin(6.28*yday(as.Date(obs_covar_mat_d$col_dates))/365)
cos_doy <- cos(6.28*yday(as.Date(obs_covar_mat_d$col_dates))/365)
#wind_effect_eco <-log(obs_covar_mat_d$ts_matrix[1,]+1)
#wind_effect_lst <- log(obs_covar_mat_d$ts_matrix[4,]+1)
#day_time_eco <- pmax(0,cos(6.28*(obs_covar_mat_d$ts_matrix[2,]-14)/24))*exp(-0.5*obs_covar_mat_d$ts_matrix[1,])
#day_time_lst <- pmax(0,cos(6.28*(obs_covar_mat_d$ts_matrix[3,]-14)/24))*exp(-0.5*obs_covar_mat_d$ts_matrix[4,])
d <- matrix(data = NA, nrow = 2, ncol = length(sin_doy))
d[1,] <- sin_doy; d[2,] <- cos_doy
# d[3,] <- wind_effect_eco; d[4,] <- wind_effect_lst
# d[5,] <- day_time_eco; d[6,] <- day_time_lst
d[is.na(d)] <- 0

####################################################
### Define covariates effects
####################################################
D <- matrix(list(0),ncol = 2, nrow = k_obs)
# in situ
c1_site <- paste0("c1_", seq_len(n_site))
c2_site <- paste0("c2_", seq_len(n_site))
zero_site <- rep(list(0), n_site)
D_mat_insitu <- cbind(c1_site,c2_site)# %>% 
# cbind(zero_site) %>% cbind(zero_site) %>% 
# cbind(zero_site) %>% cbind(zero_site)
insitu_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="insitu_sst" ]
D[insitu_rows,] <- D_mat_insitu[insitu_inds,]

# Landsat
a1_lst_site <- rep("a1_lst", n_site)
a2_lst_site <- rep("a2_lst", n_site)
D_mat_lst <- cbind(c1_site,c2_site)# ,zero_site,a1_lst_site,zero_site,a2_lst_site)
lst_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="lst_sst_clean" ]
D[lst_rows,] <- D_mat_lst[lst_inds,]


# ECOSTRESS
a1_eco_site <- rep("a1_eco", n_site)
a2_eco_site <- rep("a2_eco", n_site)
D_mat_eco <- cbind(c1_site,c2_site) #,a1_eco_site,zero_site,a2_eco_site,zero_site ) 
eco_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="eco_sst_v002_clean" ]
D[eco_rows,] <- D_mat_eco[eco_inds,]


# MODIS 
c1_site <- paste0("c1_", seq_len(n_site))
c2_site <- paste0("c2_", seq_len(n_site))
zero_site <- rep(list(0), n_site)
D_mat_modis <- cbind(c1_site,c2_site) #,zero_site,zero_site,zero_site,zero_site ) 
modis_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="modis_sst_aqua" ]
D[modis_rows,] <- D_mat_modis[modis_inds,]



# MUR - add independent seasonality for low res
c1_site <- paste0("c1_mur_", seq_len(n_site))
c2_site <- paste0("c2_mur_", seq_len(n_site))
zero_site <- rep(list(0), n_site)
D_mat_mur <- cbind(c1_site,c2_site)#,zero_site,zero_site,zero_site,zero_site ) 
mur_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="mur_sst" ]
D[mur_rows,] <- D_mat_mur[mur_inds,]



####################################################
### Define observation covariance matrix
####################################################
R <- matrix(list(0),ncol=k_obs, nrow=k_obs)
diag(R)[obs_mat_y$row_var_keys=="insitu_sst"] <- "sigma_2_insitu"
diag(R)[obs_mat_y$row_var_keys=="lst_sst_clean"] <- "sigma_2_lst"
diag(R)[obs_mat_y$row_var_keys=="eco_sst_v002_clean"] <- "sigma_2_eco"
diag(R)[obs_mat_y$row_var_keys=="modis_sst_aqua"] <- "sigma_2_modis"
diag(R)[obs_mat_y$row_var_keys=="mur_sst"] <- "sigma_2_mur"
# for (i in lst_rows) for (j in lst_rows) if (i != j) R[[i, j]] <- "cov_lst"
# for (i in eco_rows) for (j in eco_rows) if (i != j) R[[i, j]] <- "cov_eco"




####################################################
####################################################
###
### Define state model
###
####################################################
####################################################
### -------------------------------------------- ###
####################################################
### Define state matrix 
####################################################
n_chi <- n_site
n_total <- n_chi 
B <- matrix(list(0), ncol = n_total, nrow = n_total)

# set rho_fast (site leve auto correlation) at a fixedvalue
diag(B)[1:n_chi] <- rep("rho_chi",n_chi)


####################################################
### State innovations covariance 
####################################################
Q <- matrix(list(0), ncol = n_total, nrow = n_total)
diag(Q)[1:n_chi] <- rep("tau_1",n_chi)


####################################################
### State means 
####################################################
u = rep(list(0),n_total)


####################################################
### State initial conditions 
####################################################
x0 = rep(list(0),n_total)
V0 = matrix(list(0), ncol = n_total, nrow = n_total)
diag(V0) <- 1.0


# Z = Z, A = A, R = R, B = B, U = U,
# Q = Q, x0 = x0, V0 = V0
marss_model_list <- list(
  Z = Z, A = matrix(a, ncol = 1), R = R, B = B, U = matrix(u, ncol = 1),
  Q = Q, D = D, x0 = matrix(x0,ncol=1), V0 = V0,
  d = d
)

library(MARSS)
# mod <- MARSS(y_scale,marss_model_list, control = list(trace=1))
controls = list(
  trace     = 1,
  maxit     = 20,
  minit     = 5,
  abstol    = 0.001,      # loglik convergence tolerance
  conv.test.slope.tol = 0.5,
  safe      = TRUE        # slower, more robust — useful while debugging
)

fit <- NULL
plot_loglik <- c()
for (chunk in 1:30) {
  fit <- MARSS(y_scale, marss_model_list,
               inits   = if (is.null(fit)) NULL else coef(fit, type = "list"),
               control = controls)
  plot_loglik <- append(plot_loglik,fit$logLik)
  message(sprintf("chunk %2d | iters %5d | logLik %.3f",
                  chunk, fit$numIter, fit$logLik))
  saveRDS(fit, sprintf("models/marss_test/fit_chunk_%02d.rds", chunk))
  if (fit$convergence == 0) break
  plot(plot_loglik)
}
plot(plot_loglik)







summary(fit)



