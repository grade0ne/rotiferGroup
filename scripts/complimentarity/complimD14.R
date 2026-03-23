library(tidyverse)

complim_data <- read.csv("data/complim_data_d14.csv")

rotifer_popn_df <- complim_data %>%
  filter(day != 12) %>%
  mutate(clone = as.factor(clone)) %>%
  group_by(clone, replicate)

ggplot(rotifer_popn_df, aes(x = day, y = rotifers, group = interaction(clone, replicate), color = clone, shape = diversity)) + 
  geom_point() +
  geom_line() +
  theme_minimal()
