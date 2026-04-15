library(tidyverse)

data <- read.csv("data/block2_prelim1.csv")

plot(rotifers~Day, data)

ggplot(data, aes(x = Day, y = rotifers, group = interaction(clone, Replicates), color = clone)) +
  geom_point() +
  geom_line() +
  facet_wrap(~Treatment)
