
##### libraries #####

library(tidyverse)
library(magrittr)
library(ggpubr)
# library(RColorBrewer)
library(colorspace)


##### colors #####

cbPalette <- c("#E69F00", # orange
               "#56B4E9", # lightblue
               "#009E73", # green
               "#F0E442", # yellow
               "#0072B2", # darkblue
               "#D55E00", # red
               "#CC79A7") # cyan

hcl_palettes(plot = TRUE)


##### data #####

data <- 
  read_delim("output_hysteresis.csv",
             delim = ";")


##### figures #####

# data %>%
#   filter(run == 1) %>%
#   filter(time > 1000) %>%
#   mutate(species = ordered(species),
#          trophic_level = ordered(trophic_level),
#          patch = ordered(patch),
#          X = ordered(X),
#          Y = ordered(Y)) %>%
#   ggplot(aes(time, m, color = species)) +
#   geom_line() +
#   scale_color_discrete_qualitative() +
#   facet_grid(Y ~ X, labeller = "label_both")

data %>% 
  filter(run == 1) %>% 
  filter(time > 1000) %>%
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(time, environment, color = species)) +
  geom_line(linetype = 2, size = 1.5, alpha = 0.5) +
  geom_line(aes(y = genotype_mean), linetype = 1, size = 1) +
  # geom_ribbon(aes(y = genotype_mean, ymin = genotype_mean - 1.96*sqrt(genotype_var), ymax = genotype_mean + 1.96*sqrt(genotype_var), fill = species, alpha = 0.01)) +
  scale_color_discrete_qualitative() +
  labs(y = "environment -- genotype") +
  facet_grid(Y ~ X, labeller = "label_both")

data %>% 
  filter(run == 1) %>% 
  filter(time > 1000) %>%
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(time, fitness_mean, color = species)) +
  geom_line(linetype = 1, size = 1) +
  scale_color_discrete_qualitative() +
  facet_grid(Y ~ X, labeller = "label_both")

data %>% 
  filter(run == 1) %>% 
  filter(time > 1000) %>%
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(time, N, color = species)) +
  geom_line(linetype = 1, size = 1) +
  scale_color_discrete_qualitative() +
  # scale_y_log10(breaks = c(1,2,3,5,10,20,30,50,100,200,300,500,1000,2000,3000,5000,10000)) +
  facet_grid(Y ~ X, labeller = "label_both")


data %>% 
  filter(run == 1) %>% 
  filter(time > 1000) %>%
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(m, genotype_mean, color = time)) +
  geom_point(size = 1) +
  scale_color_continuous_sequential(palette = "Plasma") +
  facet_grid(Y ~ X, labeller = "label_both")


data %>% 
  filter(run == 1) %>% 
  filter(time > 1000) %>%
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(m, fitness_mean, color = time)) +
  geom_point(size = 1) +
  scale_color_continuous_sequential(palette = "Plasma") +
  facet_grid(Y ~ X, labeller = "label_both")


data %>% 
  filter(run == 1) %>% 
  filter(time > 1000) %>%
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(m, N, color = time)) +
  geom_point(size = 1) +
  # scale_y_log10() +
  scale_color_continuous_sequential(palette = "Plasma") +
  facet_grid(Y ~ X, labeller = "label_both")

