library(rstan)
library(tidyverse)
library(grid)
library(gridExtra)

## Pre-dopaminergic dataset 
# Unidim #

load('RData/SurePD_unidim_slope_beta.RData')
beta <- beta_sample$beta

placebo_rate <- mean(beta[, 1])
active_rate <- mean(beta[, 1] + beta[, 2])
time <- c(0, 0.5, 1, 1.5, 2)
dat <- data.frame(time = rep(time, 2), trt = rep(c(1, 0), each = length(time)), 
                  theta = c(active_rate*time, placebo_rate*time))
Treatment <- rep(c('Active', 'Placebo'), each = length(time))

p1 <- ggplot(data = dat, aes(x = time, y = theta, color = Treatment)) + 
  geom_point() + 
  geom_line() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        legend.title = element_blank(),
        legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal',
        plot.title = element_text(hjust = 0.5)) + 
  xlab('Visit Year') + 
  ylab('Theta') + 
  xlim(0, 2) + 
  ylim(0, 0.6) + 
  scale_color_manual(values = c('red', 'blue')) + 
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) + 
  ggtitle('Unidimensional for Pre-dopaminergic')

## multidim #
load('RData/SurePD_multidim_slope_beta.RData')
beta_f1 <- beta_f1_sample$beta_f1
beta_f2 <- beta_f2_sample$beta_f2

placebo_rate_non <- mean(beta_f1[, 1])
active_rate_non <- mean(beta_f1[, 1] + beta_f1[, 2])
placebo_rate_tre <- mean(beta_f2[, 1])
active_rate_tre <- mean(beta_f2[, 1] + beta_f2[, 2])
dat2 <- data.frame(time = rep(time, 4), trt = rep(c(1, 0), each = 2*length(time)), 
                   theta = c(active_rate_non*time, active_rate_tre*time
                             ,placebo_rate_non*time, placebo_rate_tre*time))
Color <- rep(c('Active', 'Placebo'), each = 2*length(time))
Lty = rep(c('Non-tremor', 'Tremor', 'Non-tremor', 'Tremor'), each = length(time))

p2 <- ggplot(data = dat2, aes(x = time, y = theta, linetype = Lty, color = Color)) + 
  geom_point() + 
  geom_line() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        legend.title = element_blank(),
        legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal',
        plot.title = element_text(hjust = 0.5)) + 
  xlab('Visit Year') + 
  ylab('Theta') + 
  xlim(0, 2) + 
  ylim(0, 0.6) + 
  scale_color_manual(values = c('red', 'blue')) + 
  scale_linetype_manual(values = c('longdash', 'dotted')) + 
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) + 
  ggtitle('Multidimensional for Pre-dopaminergic')


## dopaminergic dataset 

load('RData/SurePD_post_unidim_slope_beta.RData')
beta <- beta_sample$beta

placebo_rate <- mean(beta[, 1])
active_rate <- mean(beta[, 1] + beta[, 2])
time <- c(0, 0.5, 1, 1.5, 2)
dat3 <- data.frame(time = rep(time, 2), trt = rep(c(1, 0), each = length(time)), 
                   theta = c(active_rate*time, placebo_rate*time))
Treatment <- rep(c('Active', 'Placebo'), each = length(time))

p3 <- ggplot(data = dat3, aes(x = time, y = theta, color = Treatment)) + 
  geom_point() + 
  geom_line() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        legend.title = element_blank(),
        legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal',
        plot.title = element_text(hjust = 0.5)) + 
  xlab('Visit Year') + 
  ylab('Theta') + 
  xlim(0, 2) + 
  ylim(-0.6, 0.3) + 
  scale_color_manual(values = c('red', 'blue')) + 
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) + 
  ggtitle('Unidimensional for Dopaminergic')

## multidim #
load('RData/SurePD_post_multidim_slope_beta.RData')
beta_f1 <- beta_f1_sample$beta_f1
beta_f2 <- beta_f2_sample$beta_f2

placebo_rate_non <- mean(beta_f1[, 1])
active_rate_non <- mean(beta_f1[, 1] + beta_f1[, 2])
placebo_rate_tre <- mean(beta_f2[, 1])
active_rate_tre <- mean(beta_f2[, 1] + beta_f2[, 2])
dat4 <- data.frame(time = rep(time, 4), trt = rep(c(1, 0), each = 2*length(time)), 
                   theta = c(active_rate_non*time, active_rate_tre*time
                             ,placebo_rate_non*time, placebo_rate_tre*time))
Color <- rep(c('Active', 'Placebo'), each = 2*length(time))
Lty = rep(c('Non-tremor', 'Tremor', 'Non-tremor', 'Tremor'), each = length(time))

p4 <- ggplot(data = dat4, aes(x = time, y = theta, linetype = Lty, color = Color)) + 
  geom_point() + 
  geom_line() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(face = "bold", color = "#993333"),
        axis.text.y = element_text(face = "bold", color = "#993333"),
        axis.line.x = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        axis.line.y = element_line(arrow = arrow(length = unit(0.3, "cm"))),
        legend.title = element_blank(),
        legend.position='top', 
        legend.justification='left',
        legend.direction='horizontal',
        plot.title = element_text(hjust = 0.5)) + 
  xlab('Visit Year') + 
  ylab('Theta') + 
  xlim(0, 2) + 
  ylim(-0.6, 0.3) + 
  scale_color_manual(values = c('red', 'blue')) + 
  scale_linetype_manual(values = c('longdash', 'dotted')) + 
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5, 2)) + 
  ggtitle('Multidimensional for Dopaminergic')


cairo_ps(filename='plot/surepd_progression_plot_rev.eps', height=8, width = 10)
grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
dev.off()

