library(tidyverse)

mydata <- read.csv("data/fulldata.csv") %>%
  mutate(across(c(block, diversity, replicate, competition, clone), factor))

plot(rotifers~day, mydata)

library(growthrates)
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

rotifer_growth_models <- all_growthmodels(
  rotifers ~ day | block + diversity + competition + clone + replicate, 
  data  = mydata,
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
library(lme4)
library(lmerTest)
library(car)

model1 <- lmer(log(r) ~ diversity * competition + (1|block) + (1|block:clone), growth_summary)

summary(model1)
Anova(model1, type = "III")

library(emmeans)

emm <- emmeans(model1, ~ diversity | competition)
pairs(emm)

graphdata <- growth_summary %>%
  group_by(block, diversity, competition) %>%
  summarize(mean = mean((r)),
            se = sd(r)/sqrt(n()))

ggplot(graphdata, aes(x = competition, y = mean, fill = diversity)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(stat = "identity", position = "dodge", aes(ymin = mean - se, ymax = mean + se)) +
  facet_wrap(~block)

# investigate
# 

block3 <- mydata %>%
  filter(block=="3")

ggplot(block3, aes(x = day, y = rotifers, group = interaction(clone, replicate), color = replicate)) +
  geom_point(stat = "identity") +
  geom_line(stat = "identity") +
  facet_grid(competition~diversity) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(linewidth = .3))

block2 <- mydata %>%
  filter(block=="2")

ggplot(block2, aes(x = day, y = rotifers, group = interaction(replicate, clone), color = replicate)) +
  geom_point(stat = "identity") +
  geom_line(stat = "identity") +
  facet_grid(competition~diversity) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(linewidth = .3))


block1 <- mydata %>%
  filter(block=="1")

ggplot(block1, aes(x = day, y = rotifers, group = interaction(replicate, clone), color = replicate)) +
  geom_point(stat = "identity") +
  geom_line(stat = "identity") +
  facet_grid(competition~diversity) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(linewidth = .3))
