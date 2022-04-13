options(mc.cores = parallel::detectCores())

set.seed(1234)
library(dplyr)
library(tidyverse)
library(rstan)
rstan_options(auto_write = TRUE)

long.data <- read.csv("dataset/sure_pd/items.csv", stringsAsFactors = F)

N <- length(unique(long.data$PATNO))
ID.table <- as.numeric(table(long.data$PATNO))
ID <- rep(1:N, ID.table)
long.data$ID <- ID

f1.items <- 1:23
f2.items <- 24:33
n.items.f1 <- length(f1.items)
n.items.f2 <- length(f2.items)
n.items <- n.items.f1 + n.items.f2


## Consolidate items ##
## Items for part 3 data##
items <- long.data[, 16:(15+n.items)]
items <- items[, c(f1.items, f2.items)] # Re-arrange Items 
items.ncat <- sapply(1:n.items, FUN = function(x) length(table(items[, x])))
items.nobs <- lapply(1:n.items, FUN = function(x) table(items[, x]))  ## no need for consolidating

new.items <- items
for (i in 1:n.items){
  tmp.item <- items[, i]
  tmp.ncat <- items.ncat[i]
  tmp.nobs <- items.nobs[[i]]
  tmp.cat <- 0:(tmp.ncat-1)  ## Categories for i-th item: 0, 1, 2, ..., tmp.ncat - 1
  if (tmp.nobs[tmp.ncat]<10){
    tmp.item[which(tmp.item==tmp.ncat - 1)] <- tmp.ncat - 2
    new.items[, i] <- tmp.item
  }
  tmp.level <- as.numeric(names(tmp.nobs))
  if (tmp.level[tmp.ncat] - tmp.level[tmp.ncat-1]!=1){  ## if no observed score for (max level - 1)
    tmp.item[which(tmp.item==tmp.level[tmp.ncat])] <- tmp.level[tmp.ncat-1]
    new.items[, i] <- tmp.item
  }
}

items <- new.items
items.ncat <- sapply(1:n.items, FUN = function(x) length(table(items[, x])))
items.nobs <- sapply(1:n.items, FUN = function(x) table(items[, x]))

## Item 28, 32 consolidate
new.items <- items

for (i in c(28, 32)){
  tmp.item <- items[, i]
  tmp.ncat <- items.ncat[i]
  tmp.nobs <- items.nobs[[i]]
  tmp.cat <- 0:(tmp.ncat-1)  ## Categories for i-th item: 0, 1, 2, ..., tmp.ncat - 1
  tmp.item[which(tmp.item==tmp.ncat - 1)] <- tmp.ncat - 2
  new.items[, i] <- tmp.item
}

items <- new.items
items.ncat <- sapply(1:n.items, FUN = function(x) length(table(items[, x])))
items.nobs <- sapply(1:n.items, FUN = function(x) table(items[, x]))

## Create observed category index matrix ##
obs_cat <- matrix(NA, nrow(items), n.items)
for (i in 1:nrow(items)){
  which.cat <- 0
  for (k in 1:n.items){
    for (l in 0:(items.ncat[k]-1)){
      which.cat <- which.cat + 1
      if (items[i, k]==l){
        obs_cat[i, k] <- which.cat
      }
    }
  }
}

#### Multidimensional Analysis ####

md = stan_model('source_code/surepd_multidim_slope_loglik.stan')
stan_dat <- list(n = N, nobs = nrow(long.data),
                 items_f1 = n.items.f1, items_f2 = n.items.f2, items = n.items,
                 id = long.data$ID,
                 time = long.data$VisitYr, trt = long.data$TRT_ACTV,
                 ncat = items.ncat, sum_ncat = sum(items.ncat), 
                 delta_ncat = sum(items.ncat) - n.items,
                 Y = as.matrix(items), obs_cat = obs_cat, 
                 zero = as.array(rep(0, 4)))
inits1 <- list(a_temp = rep(1, n.items), 
               delta = rep(1,  sum(items.ncat) - n.items),
               b = rep(0.2, n.items),
               beta_f1 = rep(0.1, 2), 
               beta_f2 = rep(0.1, 2),
               sigma_u1 = 1, 
               sigma_u2 = 1,
               rho_01 = 0.2, rho_02 = 0.2, rho_03 = 0.2,
               rho_12 = 0.2, rho_13 = 0.2, rho_23 = 0.2)
inits2 <- list(a_temp = rep(-1, n.items), 
               delta = rep(0.5,  sum(items.ncat) - n.items),
               b = rep(0.5, n.items),
               beta_f1 = rep(0.3, 2), 
               beta_f2 = rep(0.3, 2),
               sigma_u1 = 0.5, 
               sigma_u2 = 0.5,
               rho_01 = 0.5, rho_02 = 0.5, rho_03 = 0.5,
               rho_12 = 0.5, rho_13 = 0.5, rho_23 = 0.5)
inits <- list(c1 = inits1, c2 = inits2)
pars <- c('a_temp', 'delta','location', 'b', 'beta_f1', 'beta_f2', 'sigma_u1', 
          'sigma_u2', 'rho_01', 'rho_02', 'rho_03', 
          'rho_12', 'rho_13', 'rho_23', 'Sigma_u', 'U', 'll')
fitStan <- sampling(md, data = stan_dat, iter = 3000, warmup = 2000, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 123,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
save(list = 'fitStan', file = 'RData/SurePD_multidim_slope.RData')
summ <- summary(fitStan)$summary[1:(sum(items.ncat)*2+2*2+2+6+4*4), c(1, 3, 4, 8, 9, 10)]
write.csv(summ, file = 'result/SurePD_multidim_slope.csv')

## Extract location and discrimination parameters ##

location_sample <- extract(fitStan, pars = 'location')
b_sample <- extract(fitStan, pars = 'b')

location_mean <- colMeans(location_sample$location)
b_mean <- colMeans(b_sample$b)

discri_diff_est <- matrix(NA, nrow = n.items, ncol = 5)
discri_diff_est[, 1] <- b_mean

start.location <- 0
end.location <- 0
for (k in 1:n.items){
  start.location <- end.location+1
  end.location <- end.location+items.ncat[k]-1
  discri_diff_est[k, 2:items.ncat[k]] <- location_mean[start.location:end.location]
}

# rownames(discri_diff_est) <- colnames(items.rename)
colnames(discri_diff_est) <- c('Discrim', '0 to 1', '1 to 2', '2 to 3', '3 to 4')
write.csv(discri_diff_est, 'result/diff_discri_SurePD_multidim_slope.csv')

## Extract samples ##
a_temp_sample <- extract(fitStan, 'a_temp')
delta_sample <- extract(fitStan, 'delta')
location_sample <- extract(fitStan, 'location')
b_sample <- extract(fitStan, 'b')
beta_f1_sample <- extract(fitStan, 'beta_f1')
beta_f2_sample <- extract(fitStan, 'beta_f2')
U_sample <- extract(fitStan, 'U')
ll_sample <- extract(fitStan, 'll')

save(list = c('beta_f1_sample', 'beta_f2_sample'), file = 'RData/SurePD_multidim_slope_beta.RData')

beta_f1 <- beta_f1_sample$beta_f1
beta_f2 <- beta_f2_sample$beta_f2

kk <- beta_f2[, 1] + beta_f2[, 2]
mean(kk)
quantile(kk, c(0.025, 0.975))

## Compute posterior mean
a_temp <- colMeans(a_temp_sample$a_temp)
delta <- colMeans(delta_sample$delta)
location <- colMeans(location_sample$location)
b <- colMeans(b_sample$b)
beta_f1 <- colMeans(beta_f1_sample$beta_f1)
beta_f2 <- colMeans(beta_f2_sample$beta_f2)
U <- colMeans(U_sample$U)
ll <- colMeans(ll_sample$ll)

Dbar <- -2*sum(ll)  ## Dbar

## Dhat
id <- long.data$ID
ncat <- items.ncat
sum_ncat <- sum(ncat)
nobs <- nrow(long.data)
time <- long.data$VisitYr
trt <- long.data$TRT_ACTV

obs_cat_prob <- matrix(NA, nobs, n.items)

inv_logit <- function(x){
  return(1/(1+exp(-x)))
}

for (i in 1:nobs){
  psi_prob <- cat_prob <- rep(NA, sum_ncat)
  which_cat <- 0
  which_location <- 0
  theta_f1 = beta_f1[1]*time[i] + beta_f1[2]*time[i]*trt[i] + U[id[i], 1] + U[id[i], 2]*time[i]
  theta_f2 = beta_f2[1]*time[i] + beta_f2[2]*time[i]*trt[i] + U[id[i], 3] + U[id[i], 4]*time[i]
  
  for (k in 1:n.items.f1){
    for (l in 1:(ncat[k]-1)){ 
      which_cat <- which_cat+1
      which_location <- which_location+1
      if (l==1){
        psi_prob[which_cat] = inv_logit(location[which_location] - b[k]*theta_f1)
        cat_prob[which_cat] = psi_prob[which_cat]
      } else {
        psi_prob[which_cat] = inv_logit(location[which_location] - b[k]*theta_f1)
        cat_prob[which_cat] = psi_prob[which_cat] - psi_prob[which_cat-1]
      }
    }
    which_cat = which_cat+1
    psi_prob[which_cat] = 1
    cat_prob[which_cat] = 1 - psi_prob[which_cat-1]
    obs_cat_prob[i, k] = cat_prob[obs_cat[i, k]]  
  }
  
  for (k in (n.items.f1+1):n.items){
    for (l in 1:(ncat[k]-1)){
      which_cat <- which_cat+1
      which_location <- which_location+1
      if (l==1){
        psi_prob[which_cat] = inv_logit(location[which_location] - b[k]*theta_f2)
        cat_prob[which_cat] = psi_prob[which_cat]
      } else {
        psi_prob[which_cat] = inv_logit(location[which_location] - b[k]*theta_f2)
        cat_prob[which_cat] = psi_prob[which_cat] - psi_prob[which_cat-1]
      }
    }
    which_cat = which_cat+1
    psi_prob[which_cat] = 1 
    cat_prob[which_cat] = 1 - psi_prob[which_cat-1]
    obs_cat_prob[i, k] = cat_prob[obs_cat[i, k]] 
  }
}

ll_posterior <- rowSums(log(obs_cat_prob))

## Compute DIC
np <- sum(items.ncat) + 2*2 + 2+6 ## number of parameters
log_predictive_probability <- sum(ll_posterior) + (np/2)*log(2*pi)
Dhat <- -2*sum(ll_posterior)

DIC <- 2*Dbar - Dhat

save(list = c('DIC', 'Dbar', 'Dhat', 'np', 'log_predictive_probability'), 
     file = 'RData/DIC_SurePD_multidim_slope.RData')



#### Unidimensional analysis ####
md = stan_model('source_code/surepd_unidim_slope_loglik.stan')
stan_dat <- list(n = N, nobs = nrow(long.data),
                 items = n.items,
                 id = long.data$ID,
                 time = long.data$VisitYr, trt = long.data$TRT_ACTV,
                 ncat = items.ncat, sum_ncat = sum(items.ncat), 
                 delta_ncat = sum(items.ncat) - n.items,
                 Y = as.matrix(items), obs_cat = obs_cat, 
                 zero = as.array(rep(0, 2)))

inits1 <- list(a_temp = rep(1, n.items), 
               delta = rep(1,  sum(items.ncat) - n.items),
               b = rep(0.2, n.items),
               beta = rep(0.1, 2), 
               sigma_u1 = 1, 
               rho_01 = 0.2)
inits2 <- list(a_temp = rep(-1, n.items), 
               delta = rep(0.5,  sum(items.ncat) - n.items),
               b = rep(0.5, n.items),
               beta = rep(0.3, 2), 
               sigma_u1 = 0.5, 
               rho_01 = 0.5)
inits <- list(c1 = inits1, c2 = inits2)
pars <- c('a_temp', 'delta','location', 'b', 'beta', 'sigma_u1', 'rho_01', 'U', 'll')
fitStan <- sampling(md, data = stan_dat, iter = 3000, warmup = 2000, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 123,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
save(list = 'fitStan', file = 'RData/SurePD_unidim_slope.RData')
summ <- summary(fitStan)$summary[1:(sum(items.ncat)*2+2+1+1), c(1, 3, 4, 8, 9, 10)]
write.csv(summ, file = 'result/SurePD_unidim_slope.csv')

## Extract location and discrimination parameters ##

location_sample <- extract(fitStan, pars = 'location')
b_sample <- extract(fitStan, pars = 'b')

location_mean <- colMeans(location_sample$location)
b_mean <- colMeans(b_sample$b)

discri_diff_est <- matrix(NA, nrow = n.items, ncol = 5)
discri_diff_est[, 1] <- b_mean

start.location <- 0
end.location <- 0
for (k in 1:n.items){
  start.location <- end.location+1
  end.location <- end.location+items.ncat[k]-1
  discri_diff_est[k, 2:items.ncat[k]] <- location_mean[start.location:end.location]
}

# rownames(discri_diff_est) <- colnames(items.rename)
colnames(discri_diff_est) <- c('Discrim', '0 to 1', '1 to 2', '2 to 3', '3 to 4')
write.csv(discri_diff_est, 'result/diff_discri_SurePD_unidim_slope.csv')

## Extract samples ##
a_temp_sample <- extract(fitStan, 'a_temp')
delta_sample <- extract(fitStan, 'delta')
location_sample <- extract(fitStan, 'location')
b_sample <- extract(fitStan, 'b')
beta_sample <- extract(fitStan, 'beta')
U_sample <- extract(fitStan, 'U')
ll_sample <- extract(fitStan, 'll')

save(list = c('beta_sample'), file = 'RData/SurePD_unidim_slope_beta.RData')

## Compute posterior mean
a_temp <- colMeans(a_temp_sample$a_temp)
delta <- colMeans(delta_sample$delta)
location <- colMeans(location_sample$location)
b <- colMeans(b_sample$b)
beta <- colMeans(beta_sample$beta)
U <- colMeans(U_sample$U)
ll <- colMeans(ll_sample$ll)

Dbar <- -2*sum(ll)  ## Dbar

## Dhat 

obs_cat_prob <- matrix(NA, nobs, n.items)

for (i in 1:nobs){
  psi_prob <- cat_prob <- rep(NA, sum_ncat)
  which_cat <- 0
  which_location <- 0
  theta = beta[1]*time[i] + beta[2]*time[i]*trt[i] + U[id[i], 1] + U[id[i], 2]*time[i]
  
  for (k in 1:n.items){
    for (l in 1:(ncat[k]-1)){ 
      which_cat <- which_cat+1
      which_location <- which_location+1
      if (l==1){
        psi_prob[which_cat] = inv_logit(location[which_location] - b[k]*theta)
        cat_prob[which_cat] = psi_prob[which_cat]
      } else {
        psi_prob[which_cat] = inv_logit(location[which_location] - b[k]*theta)
        cat_prob[which_cat] = psi_prob[which_cat] - psi_prob[which_cat-1]
      }
    }
    which_cat = which_cat+1
    psi_prob[which_cat] = 1
    cat_prob[which_cat] = 1 - psi_prob[which_cat-1]
    obs_cat_prob[i, k] = cat_prob[obs_cat[i, k]]  
  }
}

ll_posterior <- rowSums(log(obs_cat_prob))

## Compute DIC
np <- sum(items.ncat) + 2 +1+1 ## number of parameters
log_predictive_probability <- sum(ll_posterior) + (np/2)*log(2*pi)
Dhat <- -2*sum(ll_posterior)

DIC <- 2*Dbar - Dhat

save(list = c('DIC', 'Dbar', 'Dhat', 'np', 'log_predictive_probability'), 
     file = 'RData/DIC_SurePD_unidim_slope.RData')
