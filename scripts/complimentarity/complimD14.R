library(tidyverse)

complim_data <- read.csv("data/complim_data_d16.csv")

rotifer_popn_df <- complim_data %>%
  filter(day != 12) %>%
  mutate(clone = as.factor(clone)) %>%
  group_by(clone, replicate)

ggplot(rotifer_popn_df, aes(x = day, y = rotifers, group = interaction(clone, replicate), color = competition, shape = diversity)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~clone) +
  theme_minimal()
