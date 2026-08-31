library(ggplot2)
library(dplyr)
library(reshape2)
library(lubridate)

source("R/01_format_data_for_marss.R")
source("R/legacy_helpers.R")   # pre-pipeline trim_ecostress_outliers / convert_from_K_to_C
source("R/marss_matrix_functions.R")
sites <- c("CB001","SF017","SF0013","SF001","SF048","SF051")

dat <- read.csv("data/raw/tasi_salmon_20_yrs.csv") %>% filter(point_id %in% sites)

dat_trimmed <- trim_ecostress_outliers(dat,"nanmean","eco_sst_v002_clean",time_var="time",threshold = 10.0)
dat_trimmed <- trim_ecostress_outliers(dat_trimmed$data,"nanmean","eco_sst_v002_clean",time_var="time",threshold = 3.0)
dat <- dat_trimmed$data


obs_nms <- c("insitu_sst","mur_sst")
dat <- convert_from_K_to_C(dat,c("mur_sst"),"nanmean")
  
  
obs_covar_nms <- c("lst_wind_speed_era5")

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
print("Format observations for MARSS")
# build observations ts matrix
obs_mat_y <- convert_long_format_marss(dat, obs_nms, 1,value.var = "value",
                                       variable.var = "variable",site.var = "point_id",
                                       time.var = "time", statistic.var = "stat", statistic = c("nanmean"))

y <- obs_mat_y$ts_matrix
scales <- list(mu = mean(y[1,], na.rm = T),
               sigma = sd(y[1,], na.rm = T))
y_scaled <- (y - scales$mu)/scales$sigma
dates <- list(start=min(as.Date(dat$time)), end=max(as.Date(dat$time)), dt=1)


####################################################
####################################################
###
### Define observation model
###
####################################################
print("Define observation model")
insitu_rows <- which(obs_mat_y$row_var_keys == "insitu_sst")
mur_rows <- which(obs_mat_y$row_var_keys == "mur_sst" )

####################################################
### -------------------------------------------- ###
####################################################
### Define observation matrix 
####################################################
print("Define observation matrix")
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
print("Define observation constants")
a <- rep(0,k_obs)
# in situ
mu_site <- paste0("mu_", seq_len(n_site))
insitu_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="insitu_sst" ]
a[insitu_rows] <- mu_site[insitu_inds]

# mur
mu_mur <- paste0("mu_", seq_len(n_site))
mur_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="mur_sst" ]
a[mur_rows] <- mu_mur[mur_inds]


####################################################
### Define observation covariates
####################################################
print("Define observation covariates")
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
d <- matrix(data = NA, nrow = 2, ncol = length(sin_doy))
d[1,] <- sin_doy; d[2,] <- cos_doy
d[is.na(d)] <- 0

####################################################
### Define covariates effects
####################################################
print("Define covariates effects")
D <- matrix(list(0),ncol = 2, nrow = k_obs)
# in situ
c1_site <- paste0("c1_", seq_len(n_site))
c2_site <- paste0("c2_", seq_len(n_site))
zero_site <- rep(list(0), n_site)
D_mat_insitu <- cbind(c1_site,c2_site)
insitu_inds <-obs_mat_y$row_site_index[obs_mat_y$row_var_keys=="insitu_sst" ]
D[insitu_rows,] <- D_mat_insitu[insitu_inds,]

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
print("Define observation covariance matrix")
R <- matrix(list(0),ncol=k_obs, nrow=k_obs)
diag(R)[obs_mat_y$row_var_keys=="insitu_sst"] <- "sigma_2_insitu"
diag(R)[obs_mat_y$row_var_keys=="mur_sst"] <- "sigma_2_mur"




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
print("Define state matrix ")
n_chi <- n_site
n_total <- n_chi 
B <- matrix(list(0), ncol = n_total, nrow = n_total)

# set rho_fast (site leve auto correlation) at a fixedvalue
diag(B)[1:n_chi] <- rep("rho_chi",n_chi)


####################################################
### State innovations covariance 
####################################################
print("State innovations covariance ")
Q <- matrix(list(0), ncol = n_total, nrow = n_total)
diag(Q)[1:n_chi] <- rep("tau_1",n_chi)


####################################################
### State means 
####################################################
print("State means ")
u = rep(list(0),n_total)


####################################################
### State initial conditions 
####################################################
print("State initial conditions ")
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


print("seasonal_inits ")
# get seasonal parameters 
seasonal_inits <- function(y, obs_mat_y, d) {
  print(d)
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


# initialize parameters 
print("initialize parameters ")
fit <- MARSS(y_scaled, marss_model_list,
             control = list(minit = 1, maxit = 2, safe = TRUE ))
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

print("Get 'er going")
saveRDS(obs_mat_y, "models/marss_test_mur_only/observations.rds")
saveRDS(scales, "models/marss_test_mur_only/scales.rds")
saveRDS(dates, "models/marss_test_mur_only/dates.rds")
saveRDS(list(n_site=n_site), "models/marss_test_mur_only/design.rds")

for (chunk in 1:30) {
  fit <- MARSS(y_scaled, marss_model_list,
               inits   = if (is.null(fit)) NULL else coef(fit, type = "list"),
               control = controls)
  plot_loglik <- append(plot_loglik,fit$logLik)
  message(sprintf("chunk %2d | iters %5d | logLik %.3f",
                  chunk, fit$numIter, fit$logLik))
  saveRDS(fit, sprintf("models/marss_test_mur_only/fit_chunk_%02d.rds", chunk))
  if (fit$convergence == 0) break
  plot(plot_loglik)
}
plot(plot_loglik)




