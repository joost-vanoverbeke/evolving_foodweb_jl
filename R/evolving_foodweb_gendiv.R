

##### libraries #####

library(tidyverse)
library(magrittr)
# library(ggpubr)
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
  read_delim("../results/output_X1_Y1_tl1_loci10_oe1_d2e-1_dp5e-1_Yrange2e-1_stepLocal5e-3_sex_CC4e-4_1e3.csv", delim = ";")

data <- 
  read_delim("../results/output_X1_Y1_tl1_loci20_oe1_d2e-1_dp5e-1_Yrange2e-1_stepLocal1e-4_sex_noCC.csv", delim = ";") %>% 
  bind_rows(read_delim("../results/output_X1_Y1_tl1_loci20_oe1_d2e-1_dp5e-1_Yrange2e-1_stepLocal5e-4_sex_noCC.csv", delim = ";")) %>% 
  bind_rows(read_delim("../results/output_X1_Y1_tl1_loci20_oe1_d2e-1_dp5e-1_Yrange2e-1_stepLocal1e-3_sex_noCC.csv", delim = ";")) %>% 
  bind_rows(read_delim("../results/output_X1_Y1_tl1_loci20_oe1_d2e-1_dp5e-1_Yrange2e-1_stepLocal5e-3_sex_noCC.csv", delim = ";")) %>% 
  bind_rows(read_delim("../results/output_X1_Y1_tl1_loci10_oe1_d2e-1_dp5e-1_Yrange2e-1_stepLocal1e-4_sex_noCC.csv", delim = ";")) %>% 
  bind_rows(read_delim("../results/output_X1_Y1_tl1_loci10_oe1_d2e-1_dp5e-1_Yrange2e-1_stepLocal5e-4_sex_noCC.csv", delim = ";")) %>% 
  bind_rows(read_delim("../results/output_X1_Y1_tl1_loci10_oe1_d2e-1_dp5e-1_Yrange2e-1_stepLocal1e-3_sex_noCC.csv", delim = ";")) %>% 
  bind_rows(read_delim("../results/output_X1_Y1_tl1_loci10_oe1_d2e-1_dp5e-1_Yrange2e-1_stepLocal5e-3_sex_noCC.csv", delim = ";"))

tls <- max(data$trophic_level)  
nX <- max(data$X)
nY <- max(data$Y)

data <- 
  data %>% 
  mutate(run = ordered(run),
         origin = (species+tls-1)%/%tls,
         origin_X = ((origin-1)%%nX)+1,
         origin_Y = ((origin-1)%/%nX)+1,
         origin = ordered(origin),
         origin_X = ordered(origin_X),
         origin_Y = ordered(origin_Y),
         species = ordered(species),
         trophic_level = ordered(trophic_level),
         patch = ordered(patch),
         X = ordered(X),
         Y = ordered(Y))

# data_origin <-
#   data %>% 
#   distinct(species, trophic_level, origin, origin_X, origin_Y)

data_mean <- 
  data %>% 
  group_by(grid_X, grid_Y, m, rho, e_step_CC, e_step_local, time_CC, nbr_loci, sigma_z, mu, omega_e, d, rep_type, time) %>% 
  summarise(across(c(genotype_var, fitness_mean, fitness_var), ~mean(.))) %>% 
  ungroup()


##### plots #####

data %>% 
  # filter(time > 1000) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = genotype_mean, color = run), linetype = 1, size = 1) +
  geom_line(aes(y = environment, color = run), linetype = 2, size = 1.5, alpha = 0.5) +
  scale_color_discrete_qualitative() +
  facet_grid(e_step_local ~ nbr_loci, labeller = "label_both") +
  guides(color = "none")


data %>% 
  filter(time > 1000) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = sqrt(genotype_var), color = run), linetype = 1, size = 1) +
  scale_color_discrete_qualitative() +
  facet_grid(e_step_local ~ nbr_loci, labeller = "label_both") +
  guides(color = FALSE)


data %>% 
  filter(time > 1000) %>% 
  ggplot(aes(time, fitness_mean, color = run)) +
  geom_line(linetype = 1, size = 1) +
  scale_color_discrete_qualitative() +
  facet_grid(e_step_local ~ nbr_loci, labeller = "label_both") +
  guides(color = FALSE)


data_mean %>% 
  filter(time > 1000) %>% 
  mutate(e_step_local = ordered(e_step_local),
         nbr_loci = ordered(nbr_loci)) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = sqrt(genotype_var), color = e_step_local, linetype = nbr_loci), size = 1) +
  scale_color_discrete_qualitative()


data_mean %>% 
  filter(time > 1000) %>% 
  mutate(e_step_local = ordered(e_step_local),
         nbr_loci = ordered(nbr_loci)) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = fitness_mean, color = e_step_local, linetype = nbr_loci), size = 1) +
  scale_color_discrete_qualitative()


data_mean %>% 
  filter(time > 1000) %>% 
  mutate(e_step_local = ordered(e_step_local),
         nbr_loci = ordered(nbr_loci)) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = sqrt(fitness_var), color = e_step_local, linetype = nbr_loci), size = 1) +
  scale_color_discrete_qualitative()




