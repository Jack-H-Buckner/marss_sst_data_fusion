library(MARSS)
library(dplyr)
library(reshape2)
library(ggplot2)



source("R/marss_matrix_functions.R")
model <- readRDS("models/marss_base/fit_TMB.rds")
obs <- readRDS("outputs/tasi_salmon_20_yrs_2026-08-29/marss_inputs.rds")
scales <- readRDS("models/marss_base/scales.rds")
dates <- readRDS("models/marss_base/dates.rds")
design <- readRDS("models/marss_base/design.rds")


harmonics <- matrix(model$model$fixed$d, nrow = 4)

sin1_coef = paste0("c1_", seq_len(design$n_site))
cos1_coef = paste0("c2_", seq_len(design$n_site))
sin2_coef = paste0("c3_", seq_len(design$n_site))
cos2_coef = paste0("c4_", seq_len(design$n_site))
mu_coef = paste0("mu_", seq_len(design$n_site))


# extract parameters by names 
c1 <- model$par$A[sin1_coef,]
c2 <- model$par$A[cos1_coef,]
c3 <- model$par$A[sin2_coef,]
c4 <- model$par$A[cos2_coef,]
mu_i <- model$par$A[mu_coef,]

trends <- matrix(c(c1,c2,c3,c4), nrow = design$n_site) %*% harmonics


# Build observations - factor loading matrix 
Lambda <- tri_from_names(model$par$Z) 
I <- matrix(0,nrow=design$n_site,ncol=design$n_site)
diag(I) <- 1
Z <- cbind(I,Lambda)

# extract states and reconstruct factors 
x <- Z %*% model$states
x.se <- sqrt(Z^2 %*% model$states.se^2)
sst <- scales$sigma*(x + trends + mu_i ) + scales$mu
sst.se <- sqrt(scales$sigma^2*x.se^2)

sst


sin1_mur_coef = paste0("c1_mur_", seq_len(design$n_site))
cos1_mur_coef = paste0("c2_mur_", seq_len(design$n_site))
sin2_mur_coef = paste0("c3_mur_", seq_len(design$n_site))
cos2_mur_coef = paste0("c4_mur_", seq_len(design$n_site))


c1_mur <- model$par$A[sin1_mur_coef,]
c2_mur <- model$par$A[cos1_mur_coef,]
c3_mur <- model$par$A[sin2_mur_coef,]
c4_mur <- model$par$A[cos2_mur_coef,]

mu_mur_coef = paste0("mu_mur_", seq_len(design$n_site))
mu_mur_i <- model$par$A[mu_mur_coef,]


trends_mur <- matrix(c(c1_mur,c2_mur,c3_mur,c4_mur), nrow = design$n_site) %*% harmonics

x <- Lambda %*% model$states[(design$n_site+1):nrow(model$states),]
x.se <- sqrt(Lambda^2 %*% model$states.se[(design$n_site+1):nrow(model$states),]^2)
sst_mur <- scales$sigma*(x + trends_mur + mu_mur_i ) + scales$mu
sst_mur.se <- sqrt(scales$sigma^2*x.se^2)


plot(sst_mur[3,],sst[3,])
abline(a=0,b=1)

ggplot(mapping = aes(x = 1:101, y = sst[3,4000:4100] - sst_mur[3,4000:4100]))+
  geom_line()


