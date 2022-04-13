library(data.table)
library(tidyverse)
library(rstan)

## baseline characteristics for all (pre+post)
dat <- read.csv('dataset/sure_pd/final_dat.csv')
baseline.dat <- dat[which(dat$EVENT_ID=='BL'), ]
baseline.dat <- as.data.table(baseline.dat)
table(baseline.dat$MALE)
baseline.dat[, .( mean(AGE_BASE), sd(AGE_BASE), min(AGE_BASE), max(AGE_BASE))]
sum(is.na(baseline.dat$C_EDUCCATC))
sum(baseline.dat$C_EDUCCATC==22, na.rm = T)
education_cat <- baseline.dat$C_EDUCCATC[which(!is.na(baseline.dat$C_EDUCCATC) & baseline.dat$C_EDUCCATC!=22)]
sum(education_cat<=13)
sum(education_cat<=13)/length(education_cat)
sum(education_cat>13 & education_cat<=17)
sum(education_cat>13 & education_cat<=17)/length(education_cat)
sum(education_cat==18)
sum(education_cat==18)/length(education_cat)
sum(education_cat>=19)
sum(education_cat>=19)/length(education_cat)
baseline.dat[, .( mean(C_EDUCYRS, na.rm = T), sd(C_EDUCYRS, na.rm = T), 
                 min(C_EDUCYRS, na.rm = T), max(C_EDUCYRS, na.rm = T))]
table(baseline.dat$TRT_ACTV)
table(baseline.dat$MAOBI_AB)

## Baseline characteristics table ##
dat2 <- baseline.dat %>% mutate(ID = PATNO, Trt = TRT_ACTV, Age = AGE_BASE, 
                                Sex = MALE, 
                                Education = C_EDUCYRS, Education_cat = C_EDUCCATC, 
                                MAOBI_AB = MAOBI_AB, 
                                HY = C_NHY, UPDRS_II = updrs3, Urate = urate_valn, 
                                PD_Med = C_PDMEDYN) %>% 
  select(ID, Trt, Age, Sex, Education, MAOBI_AB, 
         HY, Urate, PD_Med, UPDRS_II)

UPDRS_III_Scores = baseline.dat[, 16:(16+32)]
UPDRS_III = rowSums(UPDRS_III_Scores)
Resting_Tremor = rowSums(UPDRS_III_Scores[, 28:32])
dat2 <- cbind(dat2, UPDRS_III = UPDRS_III, Resting_Tremor = as.numeric(Resting_Tremor>0))
active.index = which(dat2$Trt==1)
placebo.index = which(dat2$Trt==0)

dat3 <- dat2
cat_var = c(4, 6, 9, 12)
summary_stat = data.frame(g1 = rep('', 10), g2 = rep('', 10))
for (i in 3:12){
  dat4 = cbind(dat3[, 1:2], dat3[, i])
  colnames(dat4)[3] <- 'V3'
  if (i %in% cat_var){
    g1 <- dat4[placebo.index, 3]
    tmp = paste0(sum(g1==1, na.rm = T), ' (', round(sum(g1==1, na.rm = T)/length(g1)*100, digits = 1), '%)')
    summary_stat[i-2, 1] <- tmp
    
    g2 <- dat4[active.index, 3]
    tmp = paste0(sum(g2==1, na.rm = T), ' (', round(sum(g2==1, na.rm = T)/length(g2)*100, digits = 1), '%)')
    summary_stat[i-2, 2] <- tmp
  } else {
    g1 <- dat4[placebo.index, 3]
    tmp <- paste0(round(mean(g1, na.rm = T), digits = 1) , "pm", round(sd(g1, na.rm = T), digits = 1))
    summary_stat[i-2, 1] <- tmp
    
    g2 <- dat4[active.index, 3]
    tmp <- paste0(round(mean(g2, na.rm = T), digits = 1) , "pm", round(sd(g2, na.rm = T), digits = 1))
    summary_stat[i-2, 2] <- tmp
  }
}

rownames(summary_stat) <- c('Age (years)', 'Gender (Male)', 'Education (years)', 
                            'MAO-B Inhibitor', 'H&Y Scale', 'Urate Level (mg/dl)', 
                            'PD Medication', 'UPDRS-II Score', 'UPDRS-III Score', 
                            'Resting Tremor')
colnames(summary_stat) <- c('Placebo (149, 50%)', 'Inosine (149, 50%)')

write.csv(summary_stat, 'report/surepd_bl_summ.csv', row.names = T)

# 
load('RData/SurePD_multidim.RData')
b <- extract(fitStan, 'b')$b
b_mean <- colMeans(b)
b_tre <- b[, 24:33]
cor.mat.tre <- cor(b_tre)
rownames(cor.mat.tre) <- colnames(cor.mat.tre) <- 
  c('3.15a', '3.15b', '3.16a', '3.16b', '3.17a', '3.17b', 
    '3.17c', '3.17d', '3.17e', '3.18')
write.csv(cor.mat.tre, file = 'result/surepd_multidim_cor_mat_tre.csv')


# 
load('RData/SurePD_unidim.RData')
b <- extract(fitStan, 'b')$b
b_mean <- colMeans(b)
b_tre <- b[, 24:33]
cor.mat.tre <- cor(b_tre)
rownames(cor.mat.tre) <- colnames(cor.mat.tre) <- 
  c('3.15a', '3.15b', '3.16a', '3.16b', '3.17a', '3.17b', 
    '3.17c', '3.17d', '3.17e', '3.18')
write.csv(cor.mat.tre, file = 'result/surepd_unidim_cor_mat_tre.csv')

# Unidim

load('RData/SurePD_unidim_beta.RData')
beta <- beta_sample$beta
quantity <- beta[, 1] + beta[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
load('RData/DIC_SurePD_unidim.RData')

# Unidim random slope 

load('RData/SurePD_unidim_slope_beta.RData')
beta <- beta_sample$beta
quantity <- beta[, 1] + beta[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
load('RData/DIC_SurePD_unidim_slope.RData')

# Unidim Spline

load('RData/SurePD_unidim_spline_beta.RData')
beta <- beta_sample$beta
quantity <- beta[, 1] + beta[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
quantity <- beta[, 1] + beta[, 3]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
quantity <- beta[, 1] + beta[, 2] + beta[, 3] + beta[, 4]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
load('RData/DIC_SurePD_unidim_spline.RData')

# Multidim

load('RData/SurePD_multidim_beta.RData')
beta_f1 <- beta_f1_sample$beta_f1
beta_f2 <- beta_f2_sample$beta_f2
quantity <- beta_f1[, 1] + beta_f1[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
quantity <- beta_f2[, 1] + beta_f2[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
load('RData/DIC_SurePD_multidim.RData')

# Multidim with random slope 

load('RData/SurePD_multidim_slope_beta.RData')
beta_f1 <- beta_f1_sample$beta_f1
beta_f2 <- beta_f2_sample$beta_f2
quantity <- beta_f1[, 1] + beta_f1[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
quantity <- beta_f2[, 1] + beta_f2[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
load('RData/DIC_SurePD_multidim_slope.RData')

load('RData/SurePD_multidim_slope.RData')
rho_01 <- mean(extract(fitStan, par = 'rho_01')$rho_01)
rho_02 <- mean(extract(fitStan, par = 'rho_02')$rho_02)
rho_03 <- mean(extract(fitStan, par = 'rho_03')$rho_03)
rho_12 <- mean(extract(fitStan, par = 'rho_12')$rho_12)
rho_13 <- mean(extract(fitStan, par = 'rho_13')$rho_13)
rho_23 <- mean(extract(fitStan, par = 'rho_23')$rho_23)

corr.mat <- matrix(c(1, rho_01, rho_02, rho_03, 
                     rho_01, 1, rho_12, rho_13, 
                     rho_02, rho_12, 1, rho_23, 
                     rho_03, rho_13, rho_23, 1), nrow = 4, byrow = T)
write.csv(corr.mat, file = 'result/surepd_multidim_slope_corr_mat.csv', row.names = F)


# Multidim Spline 
load('RData/SurePD_multidim_spline_beta.RData')
beta_f1 <- beta_f1_sample$beta_f1
beta_f2 <- beta_f2_sample$beta_f2
quantity <- beta_f1[, 1] + beta_f1[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
quantity <- beta_f1[, 3] + beta_f1[, 4]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
quantity <- beta_f1[, 1] + beta_f1[, 3]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
quantity <- beta_f1[, 1] + beta_f1[, 2] + beta_f1[, 3] + beta_f1[, 4]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
load('RData/DIC_SurePD_multidim_spline.RData')

# Post-levodopa unidim
load('RData/SurePD_post_unidim_beta.RData')
beta <- beta_sample$beta
quantity <- beta[, 1] + beta[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
load('RData/DIC_SurePD_post_unidim.RData')


# Unidim random slope 

load('RData/SurePD_post_unidim_slope_beta.RData')
beta <- beta_sample$beta
quantity <- beta[, 1] + beta[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
load('RData/DIC_SurePD_post_unidim_slope.RData')

# Post-levodopa Multidim

load('RData/SurePD_post_multidim_beta.RData')
beta_f1 <- beta_f1_sample$beta_f1
beta_f2 <- beta_f2_sample$beta_f2
quantity <- beta_f1[, 1] + beta_f1[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
quantity <- beta_f2[, 1] + beta_f2[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
load('RData/DIC_SurePD_post_multidim.RData')


# Multidim with random slope 

load('RData/SurePD_post_multidim_slope_beta.RData')
beta_f1 <- beta_f1_sample$beta_f1
beta_f2 <- beta_f2_sample$beta_f2
quantity <- beta_f1[, 1] + beta_f1[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
quantity <- beta_f2[, 1] + beta_f2[, 2]
mean(quantity)
quantile(quantity, c(0.025, 0.975))
load('RData/DIC_SurePD_post_multidim_slope.RData')

load('RData/SurePD_post_multidim_slope.RData')
rho_01 <- mean(extract(fitStan, par = 'rho_01')$rho_01)
rho_02 <- mean(extract(fitStan, par = 'rho_02')$rho_02)
rho_03 <- mean(extract(fitStan, par = 'rho_03')$rho_03)
rho_12 <- mean(extract(fitStan, par = 'rho_12')$rho_12)
rho_13 <- mean(extract(fitStan, par = 'rho_13')$rho_13)
rho_23 <- mean(extract(fitStan, par = 'rho_23')$rho_23)

corr.mat <- matrix(c(1, rho_01, rho_02, rho_03, 
                     rho_01, 1, rho_12, rho_13, 
                     rho_02, rho_12, 1, rho_23, 
                     rho_03, rho_13, rho_23, 1), nrow = 4, byrow = T)
write.csv(corr.mat, file = 'result/surepd_post_multidim_slope_corr_mat.csv', row.names = F)


