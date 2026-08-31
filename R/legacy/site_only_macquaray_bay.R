library(ggplot2)
library(dplyr)
library(reshape2)
library(lubridate)

source("R/01_format_data_for_marss.R")
source("R/legacy_helpers.R")   # pre-pipeline trim_ecostress_outliers / convert_from_K_to_C
source("R/marss_matrix_functions.R")
sites <- c("CB002","SF039","SF033","SF040","SF0043","SF034")

dat <- read.csv("data/raw/tasi_salmon_20_yrs.csv") %>% filter(point_id %in% sites)

dat_trimmed <- trim_ecostress_outliers(dat,"nanmean","eco_sst_v002_clean",time_var="time",threshold = 10.0)
dat_trimmed <- trim_ecostress_outliers(dat_trimmed$data,"nanmean","eco_sst_v002_clean",time_var="time",threshold = 3.0)
dat <- dat_trimmed$data

read.csv("data/raw/tasi_salmon_20_yrs.csv") %>%
  filter(!(point_id %in% c("CB002","CB001","SF032")))%>% 
  group_by(point_id) %>%
  summarize(lat = mean(lat), lon = mean(lon)) %>%
  ggplot(aes(y = lat, x = lon, label = point_id)) +
  geom_text(size = 3)

dat %>% group_by(point_id) %>%
  summarize(lat = mean(lat), lon = mean(lon)) %>%
  ggplot(aes(y = lat, x = lon, label = point_id)) +
  geom_text()

obs_nms <- c("insitu_sst","lst_sst_clean",
             "eco_sst_v002_clean","modis_sst_aqua",
             "mur_sst")

obs_covar_nms <- c("eco_hour_v002","lst_hour")


obs_nms <- c("insitu_sst","lst_sst_clean",
             "eco_sst_v002_clean","modis_sst_aqua",
             "mur_sst")

dat <- convert_from_K_to_C(dat,c("lst_sst_clean","eco_sst_v002_clean",
                                 "modis_sst_aqua","mur_sst"),"nanmean")


dat %>%  filter(variable %in% obs_nms,
                stat == "nanmean") %>%
  ggplot(aes( x = time, y = value, color = point_id)) +
  geom_point() + facet_grid(point_id~variable)+
  theme_classic()

####################################################
####################################################
###
### Format observations for MARSS
###
####################################################
####################################################
# build observations ts matrix
obs_mat_y <- convert_long_format_marss(dat, obs_nms, 1,value.var = "value",
                                       variable.var = "variable",site.var = "point_id",
                                       time.var = "time", statistic.var = "stat", statistic = c("nanmean"))



y <- obs_mat_y$ts_matrix
scales <- list(mu = mean(y[1,], na.rm = T), sigma = sd(y[1,], na.rm = T))
dates <- list(start=min(as.Date(dat$time)), end=max(as.Date(dat$time)), dt=1)
y_scaled <- (y - scales$mu)/scales$sigma

rownames(y_scaled) <- obs_mat_y$row_site_keys

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
m_factors_l4 <- 3
n_site <- length(unique(obs_mat_y$row_site_index))

# site level component
I_sites <- matrix(0, nrow=n_site,ncol = n_site)
diag(I_sites) <- 1
chi_obs <- I_sites[obs_mat_y$row_site_index,]

# level 4 error factors
Lambda <- dfa_loadings(n_site,m_factors_l4,"lambda")
eta_obs <- matrix(list(0), nrow = k_obs, ncol = m_factors_l4)
eta_obs[mur_rows,] <- Lambda[obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="mur_sst"],]

# full obs matrix 
Z <- cbind(chi_obs,eta_obs) #,f_obs

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
mu_mur <- paste0("mu_", seq_len(n_site))
mur_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="mur_sst" ]
a[mur_rows] <- mu_mur[mur_inds]


####################################################
### Define observation covariates
####################################################
### build covariates ts matrix
dat_single <- dat %>% group_by(time,point_id,variable,stat) %>%
  summarize(value = mean(value)) %>%
  filter(time %in% obs_mat_y$col_dates)
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
#day_time_eco <- 500*pmax(0,cos(6.28*(obs_covar_mat_d$ts_matrix[2,]-14)/24))*exp(-0.5*obs_covar_mat_d$ts_matrix[1,])
#day_time_lst <- pmax(0,cos(6.28*(obs_covar_mat_d$ts_matrix[3,]-14)/24))*exp(-0.5*obs_covar_mat_d$ts_matrix[4,])
d <- matrix(data = NA, nrow = 2, ncol = length(sin_doy))
d[1,] <- sin_doy; d[2,] <- cos_doy
# d[3,] <- wind_effect_eco; d[4,] <- wind_effect_lst
#d[3,] <- day_time_eco ; #d[4,] <- day_time_lst
d[is.na(d)] <- 0
rownames(d) <- c("sin_doy","cos_doy")
image(d)
####################################################
### Define covariates effects
####################################################
D <- matrix(list(0),ncol = 2, nrow = k_obs)
# in situ
c1_site <- paste0("c1_", seq_len(n_site))
c2_site <- paste0("c2_", seq_len(n_site))
zero_site <- rep(list(0), n_site)
D_mat_insitu <- cbind(c1_site,c2_site)# %>% 
insitu_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="insitu_sst" ]
D[insitu_rows,] <- D_mat_insitu[insitu_inds,]

# Landsat
a1_lst_site <- rep("a1_lst", n_site)
a2_lst_site <- rep("a2_lst", n_site)
D_mat_lst <- cbind(c1_site,c2_site)
lst_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="lst_sst_clean" ]
D[lst_rows,] <- D_mat_lst[lst_inds,]


# ECOSTRESS
a1_eco_site <- rep("a1_eco", n_site)
a2_eco_site <- rep("a2_eco", n_site)
D_mat_eco <- cbind(c1_site,c2_site) 
eco_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="eco_sst_v002_clean" ]
D[eco_rows,] <- D_mat_eco[eco_inds,]


# MODIS 
c1_site <- paste0("c1_", seq_len(n_site))
c2_site <- paste0("c2_", seq_len(n_site))
zero_site <- rep(list(0), n_site)
D_mat_modis <- cbind(c1_site,c2_site) # ,zero_site
modis_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="modis_sst_aqua" ]
D[modis_rows,] <- D_mat_modis[modis_inds,]



# MUR - add independent seasonality for low res
c1_site <- paste0("c1_", seq_len(n_site))
c2_site <- paste0("c2_", seq_len(n_site))
zero_site <- rep(list(0), n_site)
D_mat_mur <- cbind(c1_site,c2_site) 
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
for (i in lst_rows) for (j in lst_rows) if (i != j) R[[i, j]] <- "cov_lst"
for (i in eco_rows) for (j in eco_rows) if (i != j) R[[i, j]] <- "cov_eco"




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
n_eta <- m_factors_l4
n_total <- n_chi + n_eta
B <- matrix(list(0), ncol = n_total, nrow = n_total)

# set rho_fast (site leve auto correlation) at a fixedvalue
diag(B)[1:n_chi] <- rep("rho_chi",n_chi)
diag(B)[(n_chi+1):(n_chi+n_eta)] <- rep(0.5,n_eta) #"rho_f"


####################################################
### State innovations covariance 
####################################################
Q <- matrix(list(0), ncol = n_total, nrow = n_total)
diag(Q)[1:n_chi] <- rep("tau_1",n_chi)
diag(Q)[(n_chi+1):(n_chi+n_eta)] <- rep(list(1),n_eta) #"rho_f"


####################################################
### State means 
####################################################
u = rep(list(0),n_total)


####################################################
### State initial condtions 
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


### Get initial parameter estimates for means and seasonality
seasonal_inits <- function(y, obs_mat_y, d) {
  n <- length(unique(obs_mat_y$row_site_index))
  out <- data.frame(site = 1:n, mu = NA, c1 = NA, c2 = NA)
  ins <- which(obs_mat_y$row_var_keys == "mur_sst")
  
  for (r in ins) {
    s  <- obs_mat_y$row_site_index[r]
    yy <- y[r, ]
    ok <- !is.na(yy)
    if (sum(ok) < 30) next
    fit <- lm(yy[ok] ~ d[1, ok] + d[2, ok])
    out[s, 2:4] <- coef(fit)
  }
  out
}
pars <- seasonal_inits(y_scaled,obs_mat_y,d)

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



fit <- MARSS(y_scaled, marss_model_list,control = list(minit = 1, maxit = 2,safe      = TRUE  ))
init_vec <- MARSSvectorizeparam(fit)   # named, correct ordering

# 2. overwrite the seasonal entries with your lm estimates
for (i in seq_len(n_site)) {
  init_vec[paste0("A.c1_", i)] <- pars$c1[i]
  init_vec[paste0("A.c2_", i)] <- pars$c2[i]
  init_vec[paste0("A.mu_",  i)] <- pars$mu[i]
}

# 3. convert back to the inits list
inits <- MARSSvectorizeparam(fit, init_vec) 
plot_loglik <- c()
fit <- MARSS(y_scaled, marss_model_list, inits = coef(inits),
             control = controls)
plot_loglik <- append(plot_loglik,fit$logLik)

plot(plot_loglik)
saveRDS(obs_mat_y, "models/marss_test_sites_macquaray/observations.rds")
saveRDS(scales, "models/marss_test_sites_macquaray/scales.rds")
saveRDS(dates, "models/marss_test_sites_macquaray/dates.rds")
saveRDS(list(n_site=n_site,m_factors=m_factors_l4), "models/marss_test_sites_macquaray/design.rds")
for (chunk in 1:75) {
  fit <- MARSS(y_scaled, marss_model_list,
               inits   = if (is.null(fit)) NULL else coef(fit, type = "list"),
               control = controls)
  plot_loglik <- append(plot_loglik,fit$logLik)
  message(sprintf("chunk %2d | iters %5d | logLik %.3f",
                  chunk, fit$numIter, fit$logLik))
  saveRDS(fit, sprintf("models/marss_test_sites_macquaray/fit_chunk_%02d.rds", chunk))
  if (fit$convergence == 0) break
  plot(plot_loglik)
}
plot(plot_loglik)