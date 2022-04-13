## Mirror plot for observed item level proportion vs. simulated proportion
## For all visits ##
library(dplyr)
library(tidyverse)
library(gridExtra)
library(grid)
library(extraDistr)
library(rstan)

non_items <- 23 ## non-tremor items
tre_items <- 10 ## tremor items
n.items <- non_items + tre_items

long.data <- read.csv("dataset/sure_pd/items_levodopa.csv", stringsAsFactors = F)
N <- length(unique(long.data$PATNO))
ID.table <- as.numeric(table(long.data$PATNO))
ID <- rep(1:N, ID.table)
long.data$ID <- ID

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
}

items <- new.items
items.ncat <- sapply(1:n.items, FUN = function(x) length(table(items[, x])))
items.nobs <- sapply(1:n.items, FUN = function(x) table(items[, x]))

# Consolidate<10
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
}

items <- new.items
items.ncat <- sapply(1:n.items, FUN = function(x) length(table(items[, x])))
items.nobs <- sapply(1:n.items, FUN = function(x) table(items[, x]))

## Extract samples and compute posterior mean ##
load('RData/SurePD_post_multidim_slope.RData')
beta_1 <- colMeans(extract(fitStan, par = 'beta_f1')$beta_f1)
beta_2 <- colMeans(extract(fitStan, par = 'beta_f2')$beta_f2)
location <- colMeans(extract(fitStan, par = 'location')$location)
b <- colMeans(extract(fitStan, par = 'b')$b)
U <- colMeans(extract(fitStan, par = 'U')$U)
beta_non <- c(beta_1[1], beta_2[1])
beta_tre <- c(beta_1[2], beta_2[2])

inv_logit <- function(x){
  return(1/(1+exp(-x)))
}
ncat <- items.ncat
sum_ncat <- sum(ncat)

## Compute observed proportions for each item
proportion.obs <- matrix(NA, nrow = n.items, ncol = 5)
for (k in 1:n.items){
  tmp.items.table <- table(new.items[, k])/nrow(new.items)
  proportion.obs[k, 1:ncat[k]] <- tmp.items.table
}

## Simulate item-specific proportions: repeating 100 times ##
ID <- long.data$ID
VisitYr <- long.data$VisitYr
trt <- long.data$TRT_ACTV

proportion.sim.list <- list()
y.sim.list <- list()
for (q in 1:100){
  set.seed(2020*q)
  y <- matrix(NA, nrow = nrow(long.data), ncol = n.items)  ## SImulated response
  
  for (i in 1:nrow(long.data)){
    theta.non = beta_non[1]*VisitYr[i] + beta_non[2]*VisitYr[i]*trt[i] + U[ID[i], 1] + U[ID[i], 2]*VisitYr[i]
    theta.tre = beta_tre[1]*VisitYr[i] + beta_tre[2]*VisitYr[i]*trt[i] + U[ID[i], 3] + U[ID[i], 4]*VisitYr[i]
    
    which_location <- 0
    for (k in 1:non_items){
      psi_prob <- cat_prob <- rep(NA, ncat[k])  
      for (l in 1:(ncat[k]-1)){
        which_location <- which_location + 1
        psi_prob[l] <- inv_logit(location[which_location] - b[k]*theta.non)
        if (l==1) cat_prob[l] <- psi_prob[l] else cat_prob[l] <- psi_prob[l] - psi_prob[l-1]
      }
      psi_prob[ncat[k]] <- 1
      cat_prob[ncat[k]] <- 1 - psi_prob[ncat[k]-1]
      
      y[i, k] <- rcat(1, cat_prob) - 1  ## Category starts from 0
    }  ## End loop of k
    
    for (k in (non_items+1):n.items){
      psi_prob <- cat_prob <- rep(NA, ncat[k]) 
      for (l in 1:(ncat[k]-1)){
        which_location <- which_location + 1
        psi_prob[l] <- inv_logit(location[which_location] - b[k]*theta.tre)
        if (l==1) cat_prob[l] <- psi_prob[l] else cat_prob[l] <- psi_prob[l] - psi_prob[l-1]
      }
      psi_prob[ncat[k]] <- 1
      cat_prob[ncat[k]] <- 1 - psi_prob[ncat[k]-1]
      
      y[i, k] <- rcat(1, cat_prob) - 1  ## Category starts from 0
    }  ## End loop of k
  }
  
  y.sim.list[[q]] <- y
  
  proportion.sim <- matrix(NA, nrow = n.items, ncol = 5)
  for (k in 1:n.items){
    tmp.items.table <- table(y[, k])/nrow(new.items)
    proportion.sim[k, 1:ncat[k]] <- tmp.items.table
  }
  
  proportion.sim[20, 2] <- 1 - proportion.sim[20, 1]
  
  proportion.sim.list[[q]] <- proportion.sim
}

save(list = c('y.sim.list', 'proportion.sim.list'), file = 'RData/SurePD_post_sim_final_check_y_proportion.RData')

## Create Mirror Plot ##
category.data <- data.frame(x = NA, y = NA, ymin = NA, ymax = NA, group = NA, item = NA, item.name = NA)
category.data <- category.data[0, ]
item.name <- c('3.1-Speech', '3.2-Facial', 
               '3.3a-RigidNeck', '3.3b-RigidRUE', '3.3c-RigidLUE', '3.3d-RigidRLE', '3.3e-RigidLLE', 
               '3.4a-Finger_Right', '3.4b-Finger_Left', 
               '3.5a-Hand_Right', '3.5b-Hand_Left', 
               '3.6a-Pronation_Right', '3.6b-Pronation_Left',
               '3.7a-Toe_Right', '3.7b-Toe_Left', 
               '3.8a-Leg_Right', '3.8b-Leg_Left', 
               '3.9-Arising', '3.10-Gait', 
               '3.11-Freezing_Gait', '3.12-Post_Stab', 
               '3.13-Postural', '3.14-Global', 
               '3.15a-PosTrem_Right', '3.15b-PosTrem_Left', 
               '3.16a-KinTrem_Right', '3.16b-KinTrem_Left', 
               '3.17a-ResTrem_RUE', '3.17b-ResTrem_LUE', '3.17c-ResTrem_RLE', '3.17d-ResTrem_LLE', '3.17e-ResTrem_LIP',
               '3.18-Constancy')
for (k in 1:n.items){
  tmp.dat <- data.frame(x = 0:(ncat[k]-1), y = proportion.obs[k, 1:ncat[k]], 
                        ymin = NA, ymax = NA, 
                        group = rep('Obs', ncat[k]), 
                        item = rep(k, ncat[k]), 
                        item.name = rep(item.name[k], ncat[k]))
  tmp.proportion.sim <- matrix(NA, nrow = 100, ncol = ncat[k])
  for (q in 1:100){
    tmp.proportion.sim[q, ] <- proportion.sim.list[[q]][k, 1:ncat[k]]
  }
  tmp.mean <- colMeans(tmp.proportion.sim)
  tmp.sd <- sapply(1:ncat[k], FUN = function(x) sd(tmp.proportion.sim[, x]))
  tmp.dat.sim <- data.frame(x = 0:(ncat[k]-1), y = tmp.mean, 
                            ymin = tmp.mean - tmp.sd, ymax = tmp.mean + tmp.sd, 
                            group = rep('Sim', ncat[k]), 
                            item = rep(k, ncat[k]), 
                            item.name = rep(item.name[k], ncat[k]))
  tmp.dat.full <- rbind(tmp.dat, tmp.dat.sim)
  category.data <- rbind(category.data, tmp.dat.full)
}

item.sub <- list(1:9, 10:18, 19:27, 28:33)
for (j in 1:length(item.sub)){
  example.item <- item.sub[[j]]
  category.data.example <- category.data[which(category.data$item %in% example.item), ]
  
  p1 <- ggplot(data = category.data.example, 
               aes(x = x, y = y, ymin = ymin, ymax = ymax, fill = group)) + 
    geom_bar(stat = 'identity', width = 0.5, position = 'dodge') + 
    geom_errorbar(position = 'dodge', width = 0.2, aes(x = x+0.08)) + 
    scale_x_continuous(breaks = 0:(ncat[k]-1)) + 
    xlab('Category') + 
    ylab('Proportion') + 
    ylim(0, 1) + 
    theme(legend.title = element_blank()) + 
    scale_fill_manual(values = c('blue', 'red')) + 
    facet_wrap(~item.name, ncol = 3) 
  
  fname <- paste0('plot/SurePD_post_mirror_plot_', j, '.png')
  ggsave(filename = fname, plot = p1, 
         width = 8, height = 8)
}

category.data$factor <- with(category.data, 
                             ifelse(item %in% 1:non_items, 'Non-tremor', 'Tremor'))
category.data.obs <- category.data[which(category.data$group=='Obs'), ]
category.data.sim <- category.data[which(category.data$group=='Sim'), ]
category.data.melt <- data.frame(x = category.data.obs$y, 
                                 y = category.data.sim$y, 
                                 level = as.factor(category.data.obs$x), 
                                 type = category.data.obs$factor)
p2 <- ggplot(category.data.melt, aes(x = x, y = y, shape = level, color = type)) + 
  geom_point(size = 3) + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_color_manual(values = c('blue', 'red')) + 
  xlab('Observed Proportions') + 
  ylab('Simulated Proportions from 100 Replications') + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


ggsave(filename = 'plot/SurePD_post_ICC_plot.png', plot = p2, 
       width = 8, height = 8)

## Get the residual's correlation heatmap plot ##
residual.corrmat <- matrix(0, n.items, n.items)
for (q in 1:100){
  residual.mat <- new.items - y.sim.list[[q]]
  residual.corrmat <- residual.corrmat + cor(residual.mat)
}
residual.corrmat <- residual.corrmat/100

residual.corrmat[lower.tri(residual.corrmat)] <- NA
colnames(residual.corrmat) <- rownames(residual.corrmat) <- 1:n.items
residual.corrmat.melted <- reshape2::melt(residual.corrmat, na.rm = T)

value <- residual.corrmat.melted$value
value1 <- value[value!=1]
min(value1)
max(value1)

p3 <- ggplot(data = residual.corrmat.melted, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() + 
  scale_x_continuous(breaks = seq(1, n.items, by = 4)) + 
  scale_y_continuous(breaks = seq(1, n.items, by = 4)) + 
  xlab('Items') + 
  ylab('Items') + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggsave(filename = 'plot/SurePD_post_Residual_corr_plot.png', plot = p3, 
       width = 8, height = 8)


## VPC for sum scores and individual items ##
VisitYr.Nominal <- NULL
for (i in 1:N){
  tmp.long <- long.data[which(long.data$ID==i), ]
  tmp.VisitYr <- rep(NA, nrow(tmp.long))
  for (j in 1:nrow(tmp.long)){
    # if (tmp.long$EVENT_ID[j]=='BL') tmp.VisitYr[j] <- 0
    # if (tmp.long$EVENT_ID[j]=='V01') tmp.VisitYr[j] <- 0.75/12
    # if (tmp.long$EVENT_ID[j]=='V02') tmp.VisitYr[j] <- 1.5/12
    # if (tmp.long$EVENT_ID[j]=='V03') tmp.VisitYr[j] <- 3/12
    # if (tmp.long$EVENT_ID[j]=='V04') tmp.VisitYr[j] <- 6/12
    # if (tmp.long$EVENT_ID[j]=='V05') tmp.VisitYr[j] <- 9/12
    # if (tmp.long$EVENT_ID[j]=='V06') tmp.VisitYr[j] <- 1
    # if (tmp.long$EVENT_ID[j]=='V07') tmp.VisitYr[j] <- 15/12
    # if (tmp.long$EVENT_ID[j]=='V08') tmp.VisitYr[j] <- 1.5
    # if (tmp.long$EVENT_ID[j]=='V09') tmp.VisitYr[j] <- 21/12
    # if (tmp.long$EVENT_ID[j]=='V10') tmp.VisitYr[j] <- 2
    tmp.VisitYr[j] <- round(tmp.long$VisitYr[j]*4)/4
  }
  
  VisitYr.Nominal <- c(VisitYr.Nominal, tmp.VisitYr)
}

sum_score_all <- rowSums(new.items)
sum_score_non <- rowSums(new.items[, 1:non_items])
sum_score_tre <- rowSums(new.items[, (non_items+1):n.items])

sum_scores.dat <- data.frame(x = VisitYr.Nominal, 
                             y_all = sum_score_all, 
                             y_non = sum_score_non, 
                             y_tre = sum_score_tre)

compute_quantile <- function(q, curr.time, tmp.sum.score){
  l <- rep(NA, 5)
  l[1] <- q
  l[2] <- curr.time
  l[3] <- median(tmp.sum.score)
  l[c(4, 5)] <- quantile(tmp.sum.score, c(0.025, 0.975))
  return(l)
}

compute_quantile_v2 <- function(vec){
  l <- rep(NA, 3)
  l[1] <- median(vec)
  l[c(2, 3)] <- quantile(vec, c(0.025, 0.975))
  return(l)
}

sim.sum_scores_all <- sim.sum_scores_non <- sim.sum_scores_tre <- 
  data.frame(q = q, x = NA, median = NA, q025 = NA, q975 = NA)
sim.sum_scores_all <- sim.sum_scores_non <- sim.sum_scores_tre <- sim.sum_scores_all[0, ]
index <- 0
for (q in 1:100){
  tmp.y <- y.sim.list[[q]]
  for (curr.time in unique(VisitYr.Nominal)){
    index <- index + 1
    tmp.y.curr.time <- tmp.y[which(VisitYr.Nominal==curr.time), ]
    tmp.sum_score_all <- rowSums(tmp.y.curr.time)
    tmp.sum_score_non <- rowSums(tmp.y.curr.time[, 1:non_items])
    tmp.sum_score_tre <- rowSums(tmp.y.curr.time[, (non_items+1):n.items])
    sim.sum_scores_all[index, ] <- compute_quantile(q, curr.time, tmp.sum_score_all)
    sim.sum_scores_non[index, ] <- compute_quantile(q, curr.time, tmp.sum_score_non)
    sim.sum_scores_tre[index, ] <- compute_quantile(q, curr.time, tmp.sum_score_tre)
  }
}

all.sim.sum_scores <- non.sim.sum_scores <- tre.sim.sum_scores <- 
  data.frame(x = unique(VisitYr.Nominal),
             median_median = NA, median_q025 = NA, median_q975 = NA, 
             q025_median = NA, q025_q025 = NA, q025_q975 = NA,
             q975_median = NA, q975_q025 = NA, q975_q975 = NA)                           
index <- 0
for (curr.time in unique(VisitYr.Nominal)){
  index <- index + 1
  tmp.all <- sim.sum_scores_all[sim.sum_scores_all$x==curr.time, ]
  all.sim.sum_scores[index, 2:4] <- compute_quantile_v2(tmp.all[, 3])
  all.sim.sum_scores[index, 5:7] <- compute_quantile_v2(tmp.all[, 4])
  all.sim.sum_scores[index, 8:10] <- compute_quantile_v2(tmp.all[, 5])
  
  tmp.non <- sim.sum_scores_non[sim.sum_scores_non$x==curr.time, ]
  non.sim.sum_scores[index, 2:4] <- compute_quantile_v2(tmp.non[, 3])
  non.sim.sum_scores[index, 5:7] <- compute_quantile_v2(tmp.non[, 4])
  non.sim.sum_scores[index, 8:10] <- compute_quantile_v2(tmp.non[, 5])
  
  tmp.tre <- sim.sum_scores_tre[sim.sum_scores_tre$x==curr.time, ]
  tre.sim.sum_scores[index, 2:4] <- compute_quantile_v2(tmp.tre[, 3])
  tre.sim.sum_scores[index, 5:7] <- compute_quantile_v2(tmp.tre[, 4])
  tre.sim.sum_scores[index, 8:10] <- compute_quantile_v2(tmp.tre[, 5])
}


p.all <- ggplot() + 
  geom_point(data = sum_scores.dat, aes(x = x, y = y_all), color = 'blue') +
  geom_jitter(data = sum_scores.dat, aes(x = x, y = y_all), color = 'blue', width = 0.1) +
  geom_ribbon(data = all.sim.sum_scores, aes(x = x,  ymin = median_q025, ymax = median_q975), alpha = 0.5, fill = "pink") + 
  geom_line(data = all.sim.sum_scores, aes(x = x, y = median_median), color = 'purple', size = 1) + 
  geom_ribbon(data = all.sim.sum_scores, aes(x = x,  ymin = q025_q025, ymax = q025_q975), alpha = 0.5, fill = "grey70") + 
  geom_line(data = all.sim.sum_scores, aes(x = x, y = q025_median), color = 'orange', size = 1) + 
  geom_ribbon(data = all.sim.sum_scores, aes(x = x,  ymin = q975_q025, ymax = q975_q975), alpha = 0.5, fill = "grey70") + 
  geom_line(data = all.sim.sum_scores, aes(x = x, y = q975_median), color = 'orange', size = 1) + 
  theme_classic() + 
  xlab('Visit Year') + 
  ylab('Sum Scores (All Items)')

p.non <- ggplot() + 
  geom_point(data = sum_scores.dat, aes(x = x, y = y_non), color = 'blue') +
  geom_jitter(data = sum_scores.dat, aes(x = x, y = y_non), color = 'blue', width = 0.1) +
  geom_ribbon(data = non.sim.sum_scores, aes(x = x,  ymin = median_q025, ymax = median_q975), alpha = 0.5, fill = "pink") + 
  geom_line(data = non.sim.sum_scores, aes(x = x, y = median_median), color = 'purple', size = 1) + 
  geom_ribbon(data = non.sim.sum_scores, aes(x = x,  ymin = q025_q025, ymax = q025_q975), alpha = 0.5, fill = "grey70") + 
  geom_line(data = non.sim.sum_scores, aes(x = x, y = q025_median), color = 'orange', size = 1) + 
  geom_ribbon(data = non.sim.sum_scores, aes(x = x,  ymin = q975_q025, ymax = q975_q975), alpha = 0.5, fill = "grey70") + 
  geom_line(data = non.sim.sum_scores, aes(x = x, y = q975_median), color = 'orange', size = 1) + 
  theme_classic() + 
  xlab('Visit Year') + 
  ylab('Sum Scores (Non-Tremor Items)')

p.tre <- ggplot() + 
  geom_point(data = sum_scores.dat, aes(x = x, y = y_tre), color = 'blue') +
  geom_jitter(data = sum_scores.dat, aes(x = x, y = y_tre), color = 'blue', width = 0.1) +
  geom_ribbon(data = tre.sim.sum_scores, aes(x = x,  ymin = median_q025, ymax = median_q975), alpha = 0.5, fill = "pink") + 
  geom_line(data = tre.sim.sum_scores, aes(x = x, y = median_median), color = 'purple', size = 1) + 
  geom_ribbon(data = tre.sim.sum_scores, aes(x = x,  ymin = q025_q025, ymax = q025_q975), alpha = 0.5, fill = "grey70") + 
  geom_line(data = tre.sim.sum_scores, aes(x = x, y = q025_median), color = 'orange', size = 1) + 
  geom_ribbon(data = tre.sim.sum_scores, aes(x = x,  ymin = q975_q025, ymax = q975_q975), alpha = 0.5, fill = "grey70") + 
  geom_line(data = tre.sim.sum_scores, aes(x = x, y = q975_median), color = 'orange', size = 1) + 
  theme_classic() + 
  xlab('Visit Year') + 
  ylab('Sum Scores (Tremor Items)')

cairo_ps(filename='plot/SurePD_post_VPC_plot.eps', height=10, width = 6)
grid.arrange(p.all, p.non, p.tre, nrow = 3)
dev.off()