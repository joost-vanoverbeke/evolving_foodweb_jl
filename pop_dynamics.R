
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

# hcl_palettes(plot = TRUE)


##### data #####

data <- 
  read_delim("output_test_15p_3tl.csv",
             delim = ";")


##### figures #####


data %>% 
  filter(run == 1) %>% 
  # filter(time > 1000) %>%
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(time, N, color = species)) +
  geom_line(aes(linetype = trophic_level), size = 1) +
  scale_color_discrete_qualitative() +
  scale_y_log10(breaks = c(1,2,3,5,10,20,30,50,100,200,300,500,1000,2000,3000,5000,10000)) +
  facet_grid(Y ~ X, labeller = "label_both")

data %>% 
  filter(run == 1) %>% 
  # filter(time > 1000) %>%
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y)) %>% 
  ggplot(aes(time, biomass, color = species)) +
  geom_line(aes(linetype = trophic_level), size = 1) +
  scale_color_discrete_qualitative() +
  # scale_y_log10(breaks = c(1,2,3,5,10,20,30,50,100,200,300,500,1000,2000,3000,5000,10000)) +
  facet_grid(Y ~ X, labeller = "label_both")


data_spec <- 
  data %>% 
  filter(run == 1) %>% 
  # filter(time > 1000) %>%
  group_by(time, species, trophic_level) %>% 
  summarise(N = sum(N),
            biomass = sum(biomass)) %>% 
  ungroup() %>% 
  mutate(species = ordered(species),
         trophic_level = ordered(trophic_level))

data_spec %>% 
  ggplot(aes(time, N, color = species)) +
  geom_line(aes(linetype = trophic_level), size = 1) +
  scale_color_discrete_qualitative() +
  scale_y_log10(breaks = c(1,2,3,5,10,20,30,50,100,200,300,500,1000,2000,3000,5000,10000))


