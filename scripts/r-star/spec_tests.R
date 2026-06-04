absdata <- read.csv("data/r-star/spec_test_1.csv")

library(tidyverse)
library(stringr)

abs_summary <- absdata %>%
  group_by(treatment) %>%
  summarize(mean = mean(abs),
            se = sd(abs)/sqrt(length(abs)))


ggplot(abs_summary, aes(x = treatment, y = mean)) +
  geom_bar(stat="identity") +
  geom_point(data = absdata, aes(x = treatment, y = abs)) +
  geom_errorbar(stat = "identity", aes(ymin = mean - se, ymax = mean + se))


X1 <- absdata %>%
  filter(treatment == "X1")
hist(X1$abs, breaks = 10)
