library(ggplot2)
library(dplyr)
library(reshape2)
library(lubridate)
library(MARSS)
library(marssTMB)

dense_rank <- function(x) {
  if (all(is.na(x))) return(as.integer(x))
  match(x, sort(unique(x)))
}

source("R/01_format_data_for_marss.R")
source("R/marss_matrix_functions.R")
sites <- c("CB001","CB008","SF001","SF0017","SF025","SF050","SF054")
#sites <- c("CB001","CB008","SF001","SF017","SF025","SF050","SF054","SF041")
dat <- readRDS("outputs/tasi_salmon_20_yrs_2026-08-29/marss_inputs.rds")
rows_to_keep <- vapply(dat$row_site_keys, function(x){x %in% sites}, logical(1))
row_vars <- dat$row_var_keys[rows_to_keep]
row_site <- dat$row_site_keys[rows_to_keep]
row_site_index <- dense_rank(dat$row_site_index[rows_to_keep])
keep_cols <- dat$col_dates > "2013-01-01"
dates <- dat$col_dates[keep_cols]

#### Row indeces ########
# instrument indexes 
insitu_rows  <- which(row_vars == "insitu_sst")
lst_rows     <- which(row_vars == "lst_sst_clean")
eco_rows     <- which(row_vars == "eco_sst_v002_clean")
modis_rows   <- which(row_vars == "modis_sst")
mur_rows     <- which(row_vars == "mur_sst")

# site indexes 
insitu_sites_ind  <- row_site_index[insitu_rows]
lst_sites_ind     <- row_site_index[lst_rows]
eco_sites_ind    <- row_site_index[eco_rows]
modis_sites_ind  <- row_site_index[modis_rows]
mur_sites_ind    <- row_site_index[mur_rows]

# get observations
y     <- dat$ts_matrix[rows_to_keep,keep_cols]
mu    <- apply(y[mur_rows,],FUN=function(x){mean(x,na.rm=T)}, MARGIN = 1)
sigma <- apply(y[mur_rows,],FUN=function(x){sd(x,na.rm=T)}, MARGIN = 1)

y[insitu_rows,] <- (y[insitu_rows,]-mu[insitu_sites_ind])/sigma[insitu_sites_ind]
y[lst_rows,]    <- (y[lst_rows,]-mu[lst_sites_ind])/sigma[lst_sites_ind]
y[eco_rows,]    <- (y[eco_rows,]-mu[eco_sites_ind])/sigma[eco_sites_ind]
y[modis_rows,]  <- (y[modis_rows,]-mu[modis_sites_ind])/sigma[modis_sites_ind]
y[mur_rows,]    <- (y[mur_rows,]-mu[mur_sites_ind])/sigma[mur_sites_ind]


scales <- list(mu=mu,sigma=sigma)
image(y)
plot(y[2,])
####################################################
####################################################
###
### Define observation model
###
####################################################
### -------------------------------------------- ###
####################################################
### Define observation matrix 
####################################################
k_obs <- nrow(y)
m_factors_l4 <- 2
n_site <- length(unique(row_site))

# site level component
I_sites <- matrix(0, nrow=n_site,ncol = n_site)
zero_sites <- matrix(0, nrow=n_site,ncol = n_site)
diag(I_sites) <- 1
chi_obs <- I_sites[row_site_index[-mur_rows],]
chi_obs <- rbind(chi_obs,zero_sites[mur_sites_ind,])
image(chi_obs)


# level 4 error factors
Lambda <- dfa_loadings(n_site,m_factors_l4,"lambda")
eta_obs <- matrix(list(0), nrow = k_obs, ncol = m_factors_l4)
eta_obs[insitu_rows,] <- Lambda[insitu_sites_ind,]
eta_obs[lst_rows,]    <- Lambda[lst_sites_ind,]
eta_obs[eco_rows,]    <- Lambda[eco_sites_ind,]
eta_obs[modis_rows,]  <- Lambda[modis_sites_ind,]
eta_obs[mur_rows,]    <- Lambda[mur_sites_ind,]

# full obs matrix 
Z <- cbind(chi_obs,eta_obs) #,f_obs

####################################################
### Define observation constants 
####################################################
a <- rep(0,k_obs)
# in situ
mu_site <- paste0("mu_", seq_len(n_site))
a[insitu_rows] <- mu_site[insitu_sites_ind]

# landsat
mu_lst <- paste0("mu_lst+mu_", seq_len(n_site))
a[lst_rows] <- mu_lst[lst_sites_ind]

# ecostress
mu_eco <- paste0("mu_eco+mu_", seq_len(n_site))
a[eco_rows] <- mu_eco[eco_sites_ind]

# modis - independent intercepts for the low res data sources
mu_modis <- paste0("mu_modis+mu_", seq_len(n_site))
a[modis_rows] <- mu_modis[modis_sites_ind]

# mur
mu_mur <- paste0("mu_mur_", seq_len(n_site))
a[mur_rows] <- mu_mur[mur_sites_ind]
a

####################################################
### Define observation covariates
####################################################
### build covariates ts matrix
d <- dat$harmonics[,keep_cols]
image(d)

####################################################
### Define covariates effects
####################################################
D <- matrix(0,ncol = 4, nrow = k_obs)
# in situ
c1_site <- paste0("c1_", seq_len(n_site))
c2_site <- paste0("c2_", seq_len(n_site))
c3_site <- paste0("c3_", seq_len(n_site))
c4_site <- paste0("c4_", seq_len(n_site))
harminics_constants <- cbind(c1_site,c2_site,c3_site,c4_site)
D[insitu_rows,] <- harminics_constants[insitu_sites_ind,]

# Landsat
D[lst_rows,] <- harminics_constants[lst_sites_ind,]

# ECOSTRESS
D[eco_rows,] <- harminics_constants[eco_sites_ind,]

# MODIS 
D[modis_rows,] <- harminics_constants[modis_sites_ind,]

# MUR - add independent seasonality for low res
c1_site <- paste0("c1_mur_", seq_len(n_site))
c2_site <- paste0("c2_mur_", seq_len(n_site))
c3_site <- paste0("c3_mur_", seq_len(n_site))
c4_site <- paste0("c4_mur_", seq_len(n_site))
D_mat_mur <- cbind(c1_site,c2_site,c3_site,c4_site) 
D[mur_rows,] <- D_mat_mur[mur_sites_ind,]
D


####################################################
### Define observation covariance matrix
####################################################
R <- matrix(list(0),ncol=k_obs, nrow=k_obs)
diag(R)[insitu_rows] <- "sigma_2_insitu"
diag(R)[lst_rows] <- "sigma_2_lst"
diag(R)[eco_rows] <- "sigma_2_eco"
diag(R)[modis_rows] <- "sigma_2_modis"
diag(R)[mur_rows] <- "sigma_2_mur"
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
n_eta <- m_factors_l4
n_total <- n_chi + n_eta
B <- matrix(list(0), ncol = n_total, nrow = n_total)

# set rho_fast (site leve auto correlation) at a fixedvalue
diag(B)[1:n_chi] <- rep(0.75,n_chi) #rep("rho_chi",n_chi)
diag(B)[(n_chi+1):(n_chi+n_eta)] <- rep(0.99,n_eta) #rep("rho_eta",n_eta) #"rho_f"



####################################################
### State innovations covariance 
####################################################
Q <- matrix(list(0), ncol = n_total, nrow = n_total)
diag(Q)[1:n_chi] <- rep("tau_1",n_chi)
diag(Q)[(n_chi+1):(n_chi+n_eta)] <- rep(list(1),n_eta) #"rho_f"
Q

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


### Get initial parameter estimates for means and seasonality
seasonal_inits <- function(y,d,row_site_index,ins) {
  n <- length(unique(row_site_index))
  out <- data.frame(site = 1:n, mu = NA, c1 = NA, c2 = NA, c3 = NA, c4 = NA)
  for (r in ins) {
    s  <- row_site_index[r]
    yy <- y[r, ]
    ok <- !is.na(yy)
    if (sum(ok) < 30) next
    fit <- lm(yy[ok] ~ d[1, ok] + d[2, ok] + d[3, ok] + d[4, ok])
    out[s, 2:6] <- coef(fit)
  }
  out
}
pars <- seasonal_inits(y,d,row_site_index,mur_rows)

library(MARSS)
controls = list(
  trace     = 0,
  maxit     = 20,
  minit     = 5,
  abstol    = 0.001,      # loglik convergence tolerance
  conv.test.slope.tol = 0.5,
  safe      = F       # slower, more robust — useful while debugging
)


fit <- MARSS(y, marss_model_list, method = "TMB", control = list(trace = 10))

saveRDS(fit, "models/marss_base/fit_TMB.rds")
saveRDS(scales, "models/marss_base/scales.rds")
saveRDS(dates, "models/marss_base/dates.rds")
saveRDS(list(n_site=n_site,m_factors=m_factors_l4), "models/marss_base/design.rds")



fit <- MARSS(y, marss_model_list,  control = list(minit = 1, maxit = 2,safe = F  ))
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
fit <- MARSS(y, marss_model_list, inits = coef(inits),
             control = controls)
plot_loglik <- append(plot_loglik,fit$logLik)

plot(plot_loglik)
saveRDS(scales, "models/marss_base/scales.rds")
saveRDS(list(n_site=n_site,m_factors=m_factors_l4), "models/marss_test_sites/design.rds")
for (chunk in 1:75) {
  fit <- MARSS(y, marss_model_list,
               method = "TMB",
               inits   = if (is.null(fit)) NULL else coef(fit, type = "list"),
               control = controls)
  plot_loglik <- append(plot_loglik,fit$logLik)
  message(sprintf("chunk %2d | iters %5d | logLik %.3f",
                  chunk, fit$numIter, fit$logLik))
  saveRDS(fit, sprintf("models/marss_base/fit_chunk_%02d.rds", chunk))
  if (fit$convergence == 0) break
  plot(plot_loglik)
}
plot(plot_loglik)