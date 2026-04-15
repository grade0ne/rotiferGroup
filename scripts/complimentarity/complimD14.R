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
    mtext(paste("a:", this_a), side = 3, adj = 0.1, line = -3, cex = 0.5)
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


library(moments)
library(car)

## growth rate

qqp(growth_summary$r)
qqp(log(growth_summary$r))

model_log <- lm(log(r) ~ diversity + competition + diversity:competition, data = growth_summary)
anova(model_log)

library(lme4)
library(lmerTest)

model_mix <- lmer(log(r) ~ competition + (1|clone), data = growth_summary)
model_mix_log <- lmer(log(r) ~ diversity * competition + (1|clone), data = growth_summary)

# alpha

hist(growth_summary$alpha)
hist(log(growth_summary$alpha))

qqp(log(growth_summary$alpha))

alpha_model <- lm(log(alpha) ~ diversity + competition + diversity:competition, data = growth_summary)


# EQ
growth_summary <- growth_summary %>%
  mutate(eq = r / alpha)

hist(growth_summary$eq)
qqp(growth_summary$eq)

k_model <- lm(eq ~ diversity*competition, data = growth_summary)

# obs - exp for percent deviation:

clone_means <- growth_summary %>%
  group_by(competition, clone) %>%
  summarize(mean_r = mean(r),
            ci = sd(log(r))/sqrt(length(r))*1.96) %>%
    ungroup()

observed <- clone_means %>%
  filter(clone == "mix") %>%
  mutate(logr = log(mean_r)) %>%
  rename(observed = logr)

expected <- clone_means %>%
  filter(clone != "mix") %>%
  group_by(competition) %>%
  summarize(expected = log(mean(mean_r))) %>%
    ungroup()

obs_exp <- observed %>%
  left_join(expected, by = "competition") %>%
  mutate(
    obs_minus_exp = observed - expected,
    pct_diff = 100 * obs_minus_exp / expected
  )

# t tests

mix_data <- growth_summary %>%
  filter(clone == "mix") %>%
  select(competition, r)

mix_diff <- mix_data %>%
  left_join(expected, by = "competition") %>%
  mutate(diff = r - expected)


t.test(mix_diff$diff[mix_diff$competition == "FALSE"], mu = 0)

t.test(mix_diff$diff[mix_diff$competition == "TRUE"], mu = 0)

ggplot(clone_means, aes(x = clone, y = log(mean_r), color = clone == "mix")) +
  geom_hline(data = expected, aes(yintercept = expected), 
             linetype = "dashed", inherit.aes = FALSE) +
  geom_point(size = 3) +
  geom_errorbar(data = subset(clone_means, clone == "mix"), 
                stat = "identity", aes(ymin = log(mean_r) - ci, ymax = log(mean_r) + ci),
                width = 0.2) +
  facet_wrap(~competition) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red", guide = "none")) +
  labs(x = "Clone", y = "Intrinsic growth rate, ln(r)") +
  theme_classic() +
  theme(legend.location = "none")


new_graph_data <- clone_means %>%
  group_by(competition) %>%
  summarize(
    mean_mono = mean(mean_r[clone != "mix"]),
    sd_mono = sd(mean_r[clone != "mix"]),
    n_mono = length((mean_r[clone != "mix"])),
    ci_mono = sd_mono / sqrt(n_mono) * 1.96,

    mean_mix = mean_r[clone == "mix"],
    sd_mix = sd(mean_r[clone == "mix"]),
    n_mix = length((mean_r[clone == "mix"])),
    ci_mix = sd_mix / sqrt(n_mix) * 1.96
  )

ggplot(clone_means, aes(x = clone, y = log(mean_r), color = clone == "mix")) +
  geom_hline(data = expected, aes(yintercept = expected), 
             linetype = "dashed", inherit.aes = FALSE) +
  geom_point(size = 3) +
  geom_errorbar( 
                stat = "identity", aes(ymin = log(mean_r) - ci, ymax = log(mean_r) + ci),
                width = 0.2) +
  facet_wrap(~competition) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red", guide = "none")) +
  labs(x = "Clone", y = "Intrinsic growth rate, ln(r)") +
  theme_classic(base_size = 16) +
  theme(legend.location = "none")

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

ggplot(clone_means, aes(x = clone, y = log(mean_r), color = clone == "mix")) +
  geom_hline(data = expected, aes(yintercept = expected), 
             linetype = "dashed", inherit.aes = FALSE) +
  geom_point(size = 3) +
  geom_point(data = growth_summary, aes(x = clone, y = log(r))) +
  geom_errorbar(data = subset(clone_means, clone == "mix"), 
                stat = "identity", aes(ymin = log(mean_r) - ci, ymax = log(mean_r) + ci),
                width = 0.2) +
  facet_wrap(~competition) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red", guide = "none")) +
  labs(x = "Clone", y = "Intrinsic growth rate, ln(r)") +
  theme_classic() +
  theme(legend.location = "none")
library(dplyr)
library(ggplot2)

plot_data <- growth_summary %>%
  mutate(
    competition = factor(competition),
    group = if_else(clone == "mix", "mix", "a-e"),
    log_r = log(r)
  ) %>%
  group_by(competition, group) %>%
  summarize(
    mean_log_r = mean(log_r, na.rm = TRUE),
    sd_log_r = sd(log_r, na.rm = TRUE),
    n = sum(!is.na(log_r)),
    se = sd_log_r / sqrt(n),
    ci = qt(0.975, df = n - 1) * se,
    .groups = "drop"
  ) %>%
  mutate(group = factor(group, levels = c("a-e", "mix")))

ggplot(plot_data, aes(x = group, y = mean_log_r, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_log_r - ci, ymax = mean_log_r + ci),
                width = 0.15) +
  facet_wrap(~competition) +
  scale_color_manual(values = c("a-e" = "black", "mix" = "red")) +
  labs(x = "Clone group", y = "Intrinsic growth rate, ln(r)") +
  theme_classic() +
  theme(legend.position = "none")




library(dplyr)
library(ggplot2)

plot_data <- growth_summary %>%
  mutate(
    competition = factor(competition),
    group = if_else(clone == "mix", "mix", "a-e"),
    log_r = log(r)
  ) %>%
  group_by(competition, group) %>%
  summarize(
    mean_log_r = mean(log_r, na.rm = TRUE),
    sd_log_r = sd(log_r, na.rm = TRUE),
    n = sum(!is.na(log_r)),
    se = sd_log_r / sqrt(n),
    ci = qt(0.975, df = n - 1) * se,
    .groups = "drop"
  ) %>%
  mutate(group = factor(group, levels = c("a-e", "mix")))

ggplot(plot_data, aes(x = group, y = mean_log_r, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_log_r - ci, ymax = mean_log_r + ci),
                width = 0.15) +
  facet_wrap(~competition) +
  scale_color_manual(values = c("a-e" = "black", "mix" = "red")) +
  labs(x = "Clone group", y = "Intrinsic growth rate, ln(r)") +
  theme_classic() +
  theme(legend.position = "none")



library(dplyr)
library(ggplot2)

plot_data <- growth_summary %>%
  mutate(
    competition = factor(competition),
    group = if_else(clone == "mix", "mix", "a-e"),
    log_r = log(r)
  ) %>%
  group_by(competition, group) %>%
  summarize(
    mean_log_r = mean(log_r, na.rm = TRUE),
    sd_log_r = sd(log_r, na.rm = TRUE),
    n = sum(!is.na(log_r)),
    se = sd_log_r / sqrt(n),
    ci = qt(0.975, df = n - 1) * se,
    .groups = "drop"
  ) %>%
  mutate(group = factor(group, levels = c("a-e", "mix")))

ggplot(plot_data, aes(x = group, y = mean_log_r, color = group)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_log_r - ci, ymax = mean_log_r + ci),
                width = 0.15) +
  facet_wrap(~competition) +
  scale_color_manual(values = c("a-e" = "black", "mix" = "red")) +
  labs(x = "Clone group", y = "Intrinsic growth rate, ln(r)") +
  theme_classic() +
  theme(legend.position = "none")
