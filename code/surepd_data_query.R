library(data.table)
library(tidyverse)
library(gridExtra)

long.dat <- read.csv('dataset/sure_pd/surepd3_goetz_20210222.csv', stringsAsFactors = F)
#long.dat <- as.data.table(long.dat)
#long.dat[, .(.N), by = .(SCHED_VIS)]
table(long.dat$EVENT_ID)

## Remove DF (unplanned visits), SC2 (visits before baseline), SV (washout visit)
## Remove U01, U02, U03, U04 
long.dat2 <- long.dat[-which(long.dat$EVENT_ID %in% c('DF', 'SC2', 'SV', 'U01', 'U02', 'U03', 'U04')), ]
#long.dat2 <- as.data.table(long.dat2)
long.dat2.table <- as.data.table(long.dat2)
long.dat2.table[, .(.N), by = .(SCHED_VIS)]


PATNO <- unique(long.dat2$PATNO)
N <- length(PATNO)  ## number of patients ##


bl.dat.full <- long.dat2[which(long.dat2$EVENT_ID=='BL'), ]
# bl.dat <- as.data.table(bl.dat)
baseline.date.full <- rep(bl.dat.full$UPDRS_DY, table(long.dat2$PATNO))
long.dat.full <- long.dat2
long.dat.full$VisitDay <- long.dat.full$UPDRS_DY - baseline.date.full

## Baseline Characteristics ##
long.dat.full2 <- as.data.table(long.dat.full)
ans <- long.dat.full2[, .(.N, mean(VisitDay), sd(VisitDay), min(VisitDay), max(VisitDay)), by = .(SCHED_VIS)]
ans2 <- long.dat.full2[, .(.N, mean(VisitDay), sd(VisitDay), min(VisitDay), max(VisitDay)), by = .(TRT_ACTV, SCHED_VIS)]


## 
long.information.full <- matrix(NA, 1, ncol = 3) %>% as.data.frame()
colnames(long.information.full) <- c('ID', 'Total_Visits', 'Max_Visit_Days')
for (x in 1:length(PATNO)){
  tmp.long <- long.dat.full[which(long.dat.full$PATNO==PATNO[x]), ]
  long.information.full[x, ] <- rep(NA, 3)
  long.information.full$ID[x] <- PATNO[x]
  long.information.full$Total_Visits[x] <- nrow(tmp.long)
  long.information.full$Max_Visit_Days[x] <- tmp.long$VisitDay[nrow(tmp.long)]
}
mean(long.information.full$Max_Visit_Days)
sd(long.information.full$Max_Visit_Days)




## Keep visits until levodopa treatment
k <- rep(0, 2)
k1 <- NA
long.dat3 <- long.dat2[0, ]
long.dat.levodopa <- long.dat2[0, ]
for (i in 1:N){
  tmp.long <- long.dat2[which(long.dat2$PATNO==PATNO[i]), ]
  tmp.levodopa <- which(tmp.long$C_ONLDOPA==1)
  ## All NA: no levodopa;
  if (length(tmp.levodopa)==0){ 
    long.dat3 <- rbind(long.dat3, tmp.long)
  } else {
    long.dat3 <- rbind(long.dat3, tmp.long[1:(min(tmp.levodopa)), ])
    if (tmp.long$TRT_ACTV[1]==0) k[1] <- k[1]+1 else k[2] <- k[2]+1
    k1 <- c(k1, tmp.long[min(tmp.levodopa), ]$SCHED_VIS)
    long.dat.levodopa <- rbind(long.dat.levodopa, tmp.long[(min(tmp.levodopa)):(nrow(tmp.long)), ])
  }
}

table(k1[-1])

bl.dat <- long.dat3[which(long.dat3$EVENT_ID=='BL'), ]
# bl.dat <- as.data.table(bl.dat)
baseline.date <- rep(bl.dat$UPDRS_DY, table(long.dat3$PATNO))
long.dat3$VisitDay <- long.dat3$UPDRS_DY - baseline.date
write.csv(long.dat3, file ='dataset/sure_pd/final_dat.csv', row.names = F)

## determine baseline date for levodopa dataset #
levodopa.ID <- unique(long.dat.levodopa$PATNO)
baseline.date.levodopa <- NULL
long.dat.levodopa$VisitDay <- rep(NA, nrow(long.dat.levodopa))
for (i in 1:length(levodopa.ID)){
  tmp.long <- long.dat.levodopa[which(long.dat.levodopa$PATNO==levodopa.ID[i]), ]
  baseline.date.levodopa <- c(baseline.date.levodopa, rep(tmp.long$UPDRS_DY[1], nrow(tmp.long)))
}
long.dat.levodopa$VisitDay <- long.dat.levodopa$UPDRS_DY - baseline.date.levodopa

## Check with Table 1 ##
bl.dat <- as.data.table(bl.dat)
bl.dat[, .(.N, mean(AGE_BASE), sd(AGE_BASE)), by = .(TRT_ACTV)]
bl.dat[, .(.N, mean(1-MALE)), by = .(TRT_ACTV)]

bl.dat[, .(.N, mean(C_EDUCYRS, na.rm = T), sd(C_EDUCYRS, na.rm = T)), by = .(TRT_ACTV)]

## Baseline Characteristics ##
long.dat4 <- as.data.table(long.dat3)
ans <- long.dat4[, .(.N, mean(VisitDay), sd(VisitDay), min(VisitDay), max(VisitDay)), by = .(SCHED_VIS)]
ans2 <- long.dat4[, .(.N, mean(VisitDay), sd(VisitDay), min(VisitDay), max(VisitDay)), by = .(TRT_ACTV, SCHED_VIS)]


## 
long.information <- matrix(NA, 1, ncol = 3) %>% as.data.frame()
colnames(long.information) <- c('ID', 'Total_Visits', 'Max_Visit_Days')
for (x in 1:length(PATNO)){
  tmp.long <- long.dat3[which(long.dat3$PATNO==PATNO[x]), ]
  long.information[x, ] <- rep(NA, 3)
  long.information$ID[x] <- PATNO[x]
  long.information$Total_Visits[x] <- nrow(tmp.long)
  long.information$Max_Visit_Days[x] <- tmp.long$VisitDay[nrow(tmp.long)]
}
mean(long.information$Max_Visit_Days)
sd(long.information$Max_Visit_Days)

long.information2 <- matrix(NA, 1, ncol = 3) %>% as.data.frame()
colnames(long.information2) <- c('ID', 'Total_Visits', 'Max_Visit_Days')
levodopa.ID <- unique(long.dat.levodopa$PATNO)
for (x in 1:length(levodopa.ID)){
  tmp.long <- long.dat.levodopa[which(long.dat.levodopa$PATNO==levodopa.ID[x]), ]
  long.information2[x, ] <- rep(NA, 3)
  long.information2$ID[x] <- levodopa.ID[x]
  long.information2$Total_Visits[x] <- nrow(tmp.long)
  long.information2$Max_Visit_Days[x] <- tmp.long$VisitDay[nrow(tmp.long)]
}
mean(long.information2$Max_Visit_Days)
sd(long.information2$Max_Visit_Days)


## Keep complete item scores only ##
n.items <- 33
n.items.non <- 23
n.items.tre <- 10
items <- long.dat3[, 16:(15+n.items)]
not_complete <- long.dat3[which(!complete.cases(items)), ]
not_complete_id <- long.dat3[which(!complete.cases(items)), ]$PATNO
not_complete_table <- long.dat3[which(long.dat3$PATNO %in% not_complete_id), ]
long.dat5 <- long.dat3[which(complete.cases(items)), ]
length(unique(long.dat5$PATNO))
items <- long.dat5[, 16:(15+n.items)]

long.dat5$sum_scores_non <- rowSums(items[, 1:n.items.non])
long.dat5$sum_scores_tre <- rowSums(items[, (n.items.non+1):n.items])
long.dat5$sum_scores_all <- rowSums(items)
long.dat6 <- as.data.table(long.dat5)

ans3 <- long.dat6[, .(.N, mean(sum_scores_all), sd(sum_scores_all), min(sum_scores_all), max(sum_scores_all)), by = .(SCHED_VIS)]
ans4 <- long.dat6[, .(.N, mean(sum_scores_all), sd(sum_scores_all), min(sum_scores_all), max(sum_scores_all)), by = .(TRT_ACTV, SCHED_VIS)]

## Loess plot ##
long.dat5$VisitYr <- long.dat5$VisitDay/365.25
write.csv(long.dat5, file = 'dataset/sure_pd/items.csv', row.names = F)
long.dat5$Treatment <- factor(long.dat5$TRT_ACTV, levels = c(0, 1), labels = c('Placebo', 'Active'))
p <- ggplot(data = long.dat5, aes(x = VisitYr, y = sum_scores_all, group = Treatment, color = Treatment)) +
  geom_smooth(method = 'loess', show.legend = T) + 
  labs(x = 'Visit Year', y = 'Sum of scores for MDS-UPDRS Part III') + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal',
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c('blue', 'red'))

# cairo_ps(filename='plot/LOWESS_SurePD_All.eps', height=4, width = 4)
# grid.arrange(p, nrow = 1)
# dev.off()

p1 <- ggplot(data = long.dat5, aes(x = VisitYr, y = sum_scores_non, group = Treatment, color = Treatment)) +
  geom_smooth(method = 'loess', show.legend = T) + 
  labs(x = 'Visit Year', y = 'Sum of non-tremor items scores for MDS-UPDRS Part III') + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal',
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c('blue', 'red'))

p2 <- ggplot(data = long.dat5, aes(x = VisitYr, y = sum_scores_tre, group = Treatment, color = Treatment)) +
  geom_smooth(method = 'loess', show.legend = T) + 
  labs(x = 'Visit Year', y = 'Sum of tremor items scores for MDS-UPDRS Part III') + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal',
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c('blue', 'red'))

# cairo_ps(filename='plot/LOWESS_SurePD_SEP.eps', height=4, width = 8)
# grid.arrange(p1, p2, nrow = 1)
# dev.off()

## Levodopa visits ##

items.levodopa <- long.dat.levodopa[, 16:(15+n.items)]
not_complete.levodopa <- long.dat.levodopa[which(!complete.cases(items.levodopa)), ]
not_complete_id.levodopa <- long.dat.levodopa[which(!complete.cases(items.levodopa)), ]$PATNO
not_complete_table.levodopa <- long.dat.levodopa[which(long.dat.levodopa$PATNO %in% not_complete_id.levodopa), ]
long.dat5.levodopa <- long.dat.levodopa[which(complete.cases(items.levodopa)), ]
items.levodopa <- long.dat5.levodopa[, 16:(15+n.items)]

long.dat5.levodopa$sum_scores_non <- rowSums(items.levodopa[, 1:n.items.non])
long.dat5.levodopa$sum_scores_tre <- rowSums(items.levodopa[, (n.items.non+1):n.items])
long.dat5.levodopa$sum_scores_all <- rowSums(items.levodopa)
long.dat6.levodopa <- as.data.table(long.dat5.levodopa)

ans3 <- long.dat6.levodopa[, .(.N, mean(sum_scores_all), sd(sum_scores_all), min(sum_scores_all), max(sum_scores_all)), by = .(SCHED_VIS)]
ans4 <- long.dat6.levodopa[, .(.N, mean(sum_scores_all), sd(sum_scores_all), min(sum_scores_all), max(sum_scores_all)), by = .(TRT_ACTV, SCHED_VIS)]

## Loess plot ##
long.dat5.levodopa$VisitYr <- long.dat5.levodopa$VisitDay/365.25
write.csv(long.dat5.levodopa, file = 'dataset/sure_pd/items_levodopa.csv', row.names = F)
long.dat5.levodopa$Treatment <- factor(long.dat5.levodopa$TRT_ACTV, levels = c(0, 1), labels = c('Placebo', 'Active'))
p4 <- ggplot(data = long.dat5.levodopa, aes(x = VisitYr, y = sum_scores_all, group = Treatment, color = Treatment)) +
  geom_smooth(method = 'loess', show.legend = T) + 
  labs(x = 'Visit Year', y = 'Sum of scores for MDS-UPDRS Part III') + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal',
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c('blue', 'red'))

# cairo_ps(filename='plot/LOWESS_SurePD_All_Levodopa.eps', height=4, width = 4)
# grid.arrange(p, nrow = 1)
# dev.off()

p5 <- ggplot(data = long.dat5.levodopa, aes(x = VisitYr, y = sum_scores_non, group = Treatment, color = Treatment)) +
  geom_smooth(method = 'loess', show.legend = T) + 
  labs(x = 'Visit Year', y = 'Sum of non-tremor items scores for MDS-UPDRS Part III') + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal',
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c('blue', 'red'))

p6 <- ggplot(data = long.dat5.levodopa, aes(x = VisitYr, y = sum_scores_tre, group = Treatment, color = Treatment)) +
  geom_smooth(method = 'loess', show.legend = T) + 
  labs(x = 'Visit Year', y = 'Sum of tremor items scores for MDS-UPDRS Part III') + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal',
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c('blue', 'red'))

# cairo_ps(filename='plot/LOWESS_SurePD_SEP_Levodopa.eps', height=4, width = 8)
# grid.arrange(p1, p2, nrow = 1)
# dev.off()

cairo_ps(filename='plot/LOWESS_SurePD.eps', height=8, width = 12)
grid.arrange(p, p1, p2, p4, p5, p6, nrow = 2, ncol = 3)
dev.off()