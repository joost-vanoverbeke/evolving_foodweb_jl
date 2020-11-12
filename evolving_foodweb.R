

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
  read_delim("output_evolving_foodweb_sexual_inv.csv",
             delim = ";")
data <- 
  read_delim("output_evolving_foodweb_rep.csv",
             delim = ";")
# data <- 
#   read_delim("output_evolving_foodweb_evo.csv",
#              delim = ";")


##### plots #####

data %>% 
  filter(run == 1) %>% 
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(time, environment)) +
  geom_line(aes(linetype = Y), size = 1) +
  # geom_line(aes(y = resource), linetype = 2, color = "darkgrey", size = 1) +
  scale_color_discrete_qualitative() +
  facet_wrap( ~ X, labeller = "label_both", ncol = 1)


data %>% 
  filter(run == 1,
         trophic_level == 3,
         time%%100 == 0) %>% 
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(time, N, group = interaction(Y, species), color = species)) +
  geom_line(aes(y = resource), linetype = 2, color = "darkgrey", size = 1) +
  # geom_line(aes(linetype = Y), size = 1) +
  geom_line(size = 1) +
  scale_color_discrete_qualitative() +
  scale_y_log10(breaks = c(1,3,10,30,100,300,1000,3000,10000,30000,100000)) +
  facet_grid(trophic_level + Y ~ X, labeller = "label_both")


data %>% 
  filter(run == 1) %>% 
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(time, N, group = interaction(Y, species), color = species)) +
  geom_line(aes(y = resource), linetype = 2, color = "darkgrey", size = 1) +
  # geom_line(aes(linetype = Y), size = 1) +
  geom_line(size = 1) +
  scale_color_discrete_qualitative() +
  scale_y_log10(breaks = c(1,3,10,30,100,300,1000,3000,10000,30000,100000)) +
  facet_grid(trophic_level + Y ~ X, labeller = "label_both")


data %>% 
  filter(run == 1) %>% 
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(time, biomass, color = species)) +
  # geom_line(aes(y = resource), linetype = 2, color = "darkgrey", size = 1) +
  geom_line(aes(linetype = Y), size = 1) +
  scale_color_discrete_qualitative() +
  scale_y_log10(breaks = c(1,3,10,30,100,300,1000,3000,10000,30000,100000)) +
  facet_grid(trophic_level ~ X, labeller = "label_both")


# difference in mortality
x <- 1:250
y <- x
d <- 0.1
for(i in x[-1]) 
  y[i] <- y[i-1]*(1-d)
y2 <- x
d2 <- 0.03
for(i in x[-1]) 
  y2[i] <- y2[i-1]*(1-d2)

ggplot(data = NULL, aes(x = x, y = y))+
  geom_line() +
  geom_line(aes(y = y2), color = "darkgrey")



