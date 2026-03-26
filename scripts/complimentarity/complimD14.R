library(tidyverse)

complim_data <- read.csv("data/complim_data_d16.csv") %>%
  mutate(across(c(diversity, replicate, competition, clone), as.factor))

rotifer_popn_df <- complim_data %>%
  mutate(clone = as.factor(clone),
         competition = case_when(competition == TRUE ~ "Present",
                                 competition == FALSE ~ "Absent")) %>%
  group_by(clone, replicate)

ggplot(rotifer_popn_df, aes(x = day, y = rotifers, group = interaction(clone, replicate), color = competition)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~clone) +
  scale_color_manual(values = c("grey55", "coral2")) +
  labs(x = "Day", y = "Rotifers", color = "Competition") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.2)
  )

# 1: model logistic growth

library(growthrates)

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
    mtext(paste("a:", this_a), side = 3, adj = 0.1, line = -2, cex = 0.5)
  }
}

rotifer_growth_models <- all_growthmodels(
  rotifers ~ day | diversity + competition + clone + replicate, 
  data  = complim_data,
  p     = c(y0 = 5, r = 0.1, alpha = 0.0005), # initial
  upper = c(y0 = 7, r = 8, alpha = 1),
  lower = c(y0 = 0, r = 0.003, alpha = 0.0001),
  FUN   = grow_logistic_alpha
)

growth_summary <- results(rotifer_growth_models) %>%
  mutate(across(c(diversity, competition, clone), as.factor))

# 2: compare growth rate (r) with/without competition at high/low diversity
# Two-way ANOVA

model <- lm(r ~ diversity + competition, data = growth_summary)
plot(model)

library(moments)
library(car)
qqp(growth_summary$r)
qqp(log(growth_summary$r))

model_log <- lm(log(r) ~ diversity + competition + diversity:competition, data = growth_summary)

anova(model_log)

library(lme4)
library(lmerTest)

model_mix <- lmer(r ~ diversity + competition + diversity:competition + (1|clone), data = growth_summary)

# rotifers that experienced competition grew more slowly than those that did not. The effect of competition was not dependent on diversity. 27% of the total variance among replicates was attributed to differences among clones.

sub_summary <- growth_summary %>%
  group_by(competition, diversity, clone) %>%
  summarize(mean_r = mean(r), se_r = sd(r)/sqrt(length(r)))

r_comparison <- sub_summary %>%
  filter(competition == TRUE)
  
mix_mean <- 0.799
mix_se <- 0.263

clone_mean <- 0.853
clone_se <- 0.105

complim_graph <- data.frame(
  group = c("mix", "clone"),
  mean  = c(0.799, 0.853),
  se    = c(0.263, 0.105)
)

ggplot(data = complim_graph, aes(x = group, y = mean, fill = group)) +
  geom_bar(stat = 'identity') +
  geom_errorbar(stat = 'identity', aes(ymin = mean - se, ymax = mean + se),
                width = 0.5) +
  scale_fill_manual(values = c("hotpink", "purple")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 0.5)
  )
