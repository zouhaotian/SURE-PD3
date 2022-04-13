library(dplyr)
library(tidyverse)
library(nlme)


long.data <- read.csv("dataset/sure_pd/items.csv", stringsAsFactors = F)

N <- length(unique(long.data$PATNO))
ID.table <- as.numeric(table(long.data$PATNO))
ID <- rep(1:N, ID.table)
long.data$ID <- ID

non_items <- 23 ## non-tremor items
tre_items <- 10 ## tremor items
n.items <- non_items + tre_items

items <- long.data[, 16:(15+n.items)]  ## 33 items
long.data$sum_scores_all <- rowSums(items)
long.data$sum_scores_non <- rowSums(items[, 1:non_items])
long.data$sum_scores_tre <- rowSums(items[, (non_items+1):n.items])

long.data2 <- long.data
long.data2$trt <- long.data2$TRT_ACTV

## All items ##
lme.all <- lme(data = long.data2, sum_scores_all ~ VisitYr + trt + I(VisitYr*trt), random = ~ 1+VisitYr|ID)
lme.all.reduced <- lme(data = long.data2, sum_scores_all ~ VisitYr, random = ~ 1+VisitYr|ID)
summ.all <- summary(lme.all)
summary.all <- summ.all$tTable[, c(1, 2, 4, 5)]
varcorr.all <- VarCorr(lme.all)

ll.all <- logLik(lme.all, REML = F)
ll.all.reduced <- logLik(lme.all.reduced, REML = F)
ll.diff <- -2*(ll.all.reduced - ll.all)
p.value.all <- 1 - pchisq(ll.diff, df = 2)

## Non-tremor items 
lme.non <- lme(data = long.data2, sum_scores_non ~ VisitYr + trt + I(VisitYr*trt), random = ~ 1+VisitYr|ID)
lme.non.reduced <- lme(data = long.data2, sum_scores_non ~ VisitYr, random = ~ 1+VisitYr|ID)
summ.non <- summary(lme.non)
summary.non <- summ.non$tTable[, c(1, 2, 4, 5)]
varcorr.non <- VarCorr(lme.non)
ll.non <- logLik(lme.non, REML = F)
ll.non.reduced <- logLik(lme.non.reduced, REML = F)
ll.diff <- -2*(ll.non.reduced - ll.non)
p.value.non <- 1 - pchisq(ll.diff, df = 2)

## tre-tremor items 
lme.tre <- lme(data = long.data2, sum_scores_tre ~ VisitYr + trt + I(VisitYr*trt), random = ~ 1+VisitYr|ID)
lme.tre.reduced <- lme(data = long.data2, sum_scores_tre ~ VisitYr, random = ~ 1+VisitYr|ID)
summ.tre <- summary(lme.tre)
summary.tre <- summ.tre$tTable[, c(1, 2, 4, 5)]
varcorr.tre <- VarCorr(lme.tre)
ll.tre <- logLik(lme.tre, REML = F)
ll.tre.reduced <- logLik(lme.tre.reduced, REML = F)
ll.diff <- -2*(ll.tre.reduced - ll.tre)
p.value.tre <- 1 - pchisq(ll.diff, df = 2)

save(list = c('summary.all', 'summary.non', 'summary.tre',
              'p.value.all', 'p.value.non', 'p.value.tre', 
              'varcorr.all', 'varcorr.non', 'varcorr.tre'), 
     file = 'RData/surepd_linear_mixed_model.RData')
