library(tidyverse)
library(growthrates)
blk2data <- read.csv("data/block2data.csv")%>%
  mutate(across(c(diversity, replicate, competition, clone),as.factor))
view(blk2data)

# user-specified model
logistic_alpha_function <- function(time, parms) {
  y0    <- parms[["y0"]]
  r     <- parms[["r"]]
  alpha <- parms[["alpha"]]
  
  y <- (r * y0) / (alpha * y0 + (r - alpha * y0) * exp(-r * time))
  
  cbind(time = time, y = y)
}

grow_logistic_alpha <- growthmodel(
  logistic_alpha_function,
  pnames = c("y0", "r", "alpha")
)

plot_models <- function(models) {
  coefs <- as.data.frame(coef(models))
  model_ids <- c(1:length(coefs$y0))
  
  columns <- floor(sqrt(length(model_ids)))
  rows <- columns + 1
  par(mfrow = c(columns, rows), mar = c(1, 1, 1, 1))
  
  for (i in seq_along(model_ids)) {
    plot(models[[i]])
    
    this_r <- round(coefs$r[i], 6)
    this_a <- round(coefs$alpha[i], 6)
    this_model <- rownames(coefs)[i]
    
    mtext(this_model, side = 3, adj = 0.1, line = -1, cex = 0.5)
    mtext(paste("r:", this_r), side = 3, adj = 0.1, line = -2, cex = 0.5)
    mtext(paste("a:", this_a), side = 3, adj = 0.1, line = -3, cex = 0.5)
  }
}

rotifer_growth_models_2 <- all_growthmodels(
  rotifers ~ day | diversity + competition + clone + replicate, 
  data  = blk2data,
  p     = c(y0 = 5, r = 0.1, alpha = 0.0005), # initial
  upper = c(y0 = 7, r = 8, alpha = 1),
  lower = c(y0 = 0, r = 0.003, alpha = 0.0001),
  FUN   = grow_logistic_alpha
)

growth_summary <- results(rotifer_growth_models) %>%
  mutate(across(c(diversity, competition, clone), as.factor))
view(growth_summary)
#stuff for anova test 
library(moments)
library(car)

qqp(growth_summary$r)
qqp(log(growth_summary$r))

model_log <- lm(log(r) ~ diversity + competition + diversity:competition, data = growth_summary)
anova(model_log)

library(lme4)
library(lmerTest)

model_mix_log <- lmer(log(r) ~ diversity * competition + (1|clone), data = growth_summary)

# alpha

hist(growth_summary$alpha)
hist(log(growth_summary$alpha))
shapiro.test(growth_summary$alpha)
shapiro.test(log(growth_summary$alpha))
shapiro.test(sqrt(growth_summary$alpha))

qqp(sqrt(growth_summary$alpha))

alpha_model <- lm(sqrt(alpha) ~ diversity + competition + diversity:competition, data = growth_summary)


# EQ
growth_summary <- growth_summary %>%
  mutate(eq = r / alpha)

hist(growth_summary$eq)
shapiro.test(growth_summary$eq)
shapiro.test(sqrt(growth_summary$eq))
#THIS IS NOT TURNING NORMAL
qqp(growth_summary$eq)

k_model <- lm(eq ~ diversity*competition, data = growth_summary)
