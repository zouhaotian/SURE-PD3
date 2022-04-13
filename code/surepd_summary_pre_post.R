library(data.table)
library(tidyverse)
library(rstan)

## baseline characteristics for Pre- dataset
dat <- read.csv('dataset/sure_pd/items.csv')
baseline.dat <- dat[which(dat$EVENT_ID=='BL'), ]
patno <- unique(dat$PATNO)
bl.patno <- baseline.dat$PATNO

notbl.patno = patno[which(!patno %in% bl.patno)]

baseline.dat <- rbind(baseline.dat, dat[which(dat$PATNO==notbl.patno), ][1, ])
baseline.dat <- baseline.dat %>% arrange(PATNO, EVENT_ID)

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
colnames(summary_stat) <- c('Placebo (149, 50.3%)', 'Inosine (147, 49.7%)')

write.csv(summary_stat, 'report/surepd_pre_bl_summ.csv', row.names = T)


## Post dataset ##
dat <- read.csv('dataset/sure_pd/items_levodopa.csv')
baseline.dat <- dat[0, ]
uID <- unique(dat$PATNO)

for (i in 1:length(uID)){
  tmp.long <- dat[which(dat$PATNO==uID[i]), ]
  baseline.dat <- rbind(baseline.dat, tmp.long[1, ])
}

## 
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
colnames(summary_stat) <- c('Placebo (88, 58.3%)', 'Inosine (63, 41.7%)')

write.csv(summary_stat, 'report/surepd_post_bl_summ.csv', row.names = T)
