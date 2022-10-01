

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


##### data #####

data <- 
  read_delim("output_X10_Y1_3tl.csv", delim = ";") %>% 
  mutate(origin = (species+2)%/%3,
         origin_X = ((origin-1)%%5)+1,
         origin_Y = ((origin-1)%/%5)+1,
         origin = ordered(origin),
         origin_X = ordered(origin_X),
         origin_Y = ordered(origin_Y),
         species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y))

data_origin <-
  data %>% 
  distinct(species, trophic_level, origin, origin_X, origin_Y)


##### plots #####


# data %>% 
#   filter(run == 1) %>% 
#   ggplot(aes(x = time)) +
#   geom_line(aes(y = genotype_mean, color = species), linetype = 1, size = 1) +
#   geom_line(aes(y = environment), linetype = 2, size = 1.5, alpha = 0.5) +
#   scale_color_discrete_qualitative() +
#   facet_grid(Y ~ X, labeller = "label_both") +
#   guides(color = FALSE)
# 
# 
# data %>% 
#   filter(run == 1) %>% 
#   ggplot(aes(time, fitness_mean, color = species)) +
#   geom_line(linetype = 1, size = 1) +
#   scale_color_discrete_qualitative() +
#   facet_grid(Y ~ X, labeller = "label_both") +
#   guides(color = FALSE)


data %>% 
  filter(run == 1,
         time > 100) %>% 
  ggplot(aes(time, environment)) +
  geom_line(aes(linetype = Y), size = 1) +
  scale_color_discrete_qualitative() +
  facet_wrap( ~ X, labeller = "label_both", nrow = 1) +
  guides(color = FALSE, linetype = FALSE)


data %>% 
  filter(run == 1,
         time > 100) %>% 
  ggplot(aes(time, N, group = species, color = origin)) +
  geom_line(aes(linetype = trophic_level), size = 1) +
  scale_color_discrete_qualitative() +
  # scale_y_log10() +
  facet_grid(Y + trophic_level ~ X, labeller = "label_both", scales = "free_y") +
  # facet_grid(trophic_level + Y ~ X, labeller = "label_both", scales = "free_y") +
  guides(color = FALSE, linetype = FALSE)


data %>% 
  filter(run == 1,
         time > 100) %>% 
  ggplot(aes(time, biomass, group = species, color = origin)) +
  geom_line(aes(linetype = trophic_level), size = 1) +
  scale_color_discrete_qualitative() +
  # scale_y_log10() +
  facet_grid(Y + trophic_level ~ X, labeller = "label_both", scales = "free_y") +
  guides(color = FALSE, linetype = FALSE)


data %>% 
  filter(run == 1,
         time > 100) %>% 
  filter(trophic_level == 1) %>% 
  ggplot(aes(time, N, group = species, color = origin_X)) +
  geom_line(aes(linetype = origin_Y), size = 1) +
  scale_color_discrete_qualitative() +
  # scale_y_log10() +
  facet_grid(Y + trophic_level ~ X, labeller = "label_both", scales = "free_y") +
  # facet_grid(trophic_level + Y ~ X, labeller = "label_both", scales = "free_y") +
  guides(color = FALSE, linetype = FALSE)


data %>% 
  filter(run == 1,
         time > 100) %>% 
  filter(trophic_level == 2) %>% 
  ggplot(aes(time, N, group = species, color = origin_X)) +
  geom_line(aes(linetype = origin_Y), size = 1) +
  scale_color_discrete_qualitative() +
  # scale_y_log10() +
  facet_grid(Y + trophic_level ~ X, labeller = "label_both", scales = "free_y") +
  # facet_grid(trophic_level + Y ~ X, labeller = "label_both", scales = "free_y") +
  guides(color = FALSE, linetype = FALSE)


data %>% 
  filter(run == 1,
         time > 100) %>% 
  filter(trophic_level == 3) %>% 
  ggplot(aes(time, N, group = species, color = origin_X)) +
  geom_line(aes(linetype = origin_Y), size = 1) +
  scale_color_discrete_qualitative() +
  # scale_y_log10() +
  facet_grid(Y + trophic_level ~ X, labeller = "label_both", scales = "free_y") +
  # facet_grid(trophic_level + Y ~ X, labeller = "label_both", scales = "free_y") +
  guides(color = FALSE, linetype = FALSE)


# # difference in mortality
# x <- 1:250
# y <- x
# d <- 0.1
# for(i in x[-1]) 
#   y[i] <- y[i-1]*(1-d)
# y2 <- x
# d2 <- 0.03
# for(i in x[-1]) 
#   y2[i] <- y2[i-1]*(1-d2)
# 
# ggplot(data = NULL, aes(x = x, y = y))+
#   geom_line() +
#   geom_line(aes(y = y2), color = "darkgrey")



