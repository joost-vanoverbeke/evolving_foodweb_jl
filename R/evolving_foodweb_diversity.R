

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
  read_delim("../results/output_X11_Y1_tl1_d2e-1_dp5e-1_noCC.csv", delim = ";") %>% 
  mutate(scenario = "noCC",
         evolution = "evolution") %>% 
  bind_rows(read_delim("../results/output_X11_Y1_tl1_d2e-1_dp5e-1_stCC2e-4_tCC1e+3.csv", delim = ";") %>% 
              mutate(scenario = "stCC2e-4_tCC1e+3",
                     evolution = "evolution")) %>% 
  bind_rows(read_delim("../results/output_X11_Y1_tl1_d2e-1_dp5e-1_stCC1e-4_tCC1e+3.csv", delim = ";") %>% 
              mutate(scenario = "stCC1e-4_tCC1e+3",
                     evolution = "evolution")) %>% 
  bind_rows(read_delim("../results/output_X11_Y1_tl1_d2e-1_dp5e-1_stCC1e-4_tCC2e+3.csv", delim = ";") %>% 
              mutate(scenario = "stCC1e-4_tCC2e+3",
                     evolution = "evolution")) %>% 
  bind_rows(read_delim("../results/output_X11_Y1_tl1_d2e-1_dp5e-1_noEvol_noCC.csv", delim = ";") %>% 
              mutate(scenario = "noCC",
                     evolution = "no evolution")) %>% 
  bind_rows(read_delim("../results/output_X11_Y1_tl1_d2e-1_dp5e-1_noEvol_stCC2e-4_tCC1e+3.csv", delim = ";") %>% 
              mutate(scenario = "stCC2e-4_tCC1e+3",
                     evolution = "no evolution")) %>% 
  bind_rows(read_delim("../results/output_X11_Y1_tl1_d2e-1_dp5e-1_noEvol_stCC1e-4_tCC1e+3.csv", delim = ";") %>% 
              mutate(scenario = "stCC1e-4_tCC1e+3",
                     evolution = "no evolution")) %>% 
  bind_rows(read_delim("../results/output_X11_Y1_tl1_d2e-1_dp5e-1_noEvol_stCC1e-4_tCC2e+3.csv", delim = ";") %>% 
              mutate(scenario = "stCC1e-4_tCC2e+3",
                     evolution = "no evolution")) 

tls <- max(data$trophic_level)  
nX <- max(data$X)
nY <- max(data$Y)

data <- 
  data %>% 
  mutate(origin = (species+tls-1)%/%tls,
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

data_origin <-
  data %>%
  distinct(species, trophic_level, origin, origin_X, origin_Y)


data_div_alpha <-
  data %>% 
  mutate(run = ordered(run)) %>% 
  filter(N > 0) %>% 
  group_by(scenario, evolution, run, time, patch, X, Y) %>% 
  mutate(p_N = N/sum(N),
         p_biomass = biomass/sum(biomass)) %>%
  summarise(alpha_H_N = -sum(p_N*log(p_N)),
            alpha_H_biomass = -sum(p_biomass*log(p_biomass))) %>% 
  ungroup() %>% 
  mutate(# alpha_div_N = exp(alpha_H_N),
         alpha_div_biomass = exp(alpha_H_biomass)) %>%
  group_by(scenario, evolution, run, time) %>% 
  summarise(alpha_div_N = exp((-1/n())*sum(-alpha_H_N)),
            alpha_div_biomass = exp((-1/n())*sum(-alpha_H_biomass)),
            alpha_H_N = mean(alpha_H_N),
            alpha_H_biomass = mean(alpha_H_biomass),
            # alpha_div_N = mean(alpha_div_N),
            # alpha_div_biomass = mean(alpha_div_biomass)
            ) %>% 
  ungroup()
  
data_div_gamma <-
  data %>% 
  mutate(run = ordered(run)) %>% 
  filter(N > 0) %>% 
  group_by(scenario, evolution, run, time, species) %>%
  summarise(N = sum(N),
            biomass = sum(biomass)) %>% 
  ungroup() %>% 
  group_by(scenario, evolution, run, time) %>%
  mutate(p_N = N/sum(N),
         p_biomass = biomass/sum(biomass)) %>%
  summarise(gamma_H_N = -sum(p_N*log(p_N)),
            gamma_H_biomass = -sum(p_biomass*log(p_biomass))) %>% 
  ungroup() %>% 
  mutate(gamma_div_N = exp(gamma_H_N),
         gamma_div_biomass = exp(gamma_H_biomass)) %>% 
  group_by(scenario, evolution, run, time) %>% 
  summarise(gamma_H_N = mean(gamma_H_N),
            gamma_H_biomass = mean(gamma_H_biomass),
            gamma_div_N = mean(gamma_div_N),
            gamma_div_biomass = mean(gamma_div_biomass)) %>% 
  ungroup()


data_div <-
  data_div_alpha %>% 
  left_join(data_div_gamma) %>% 
  mutate(beta_H_N = gamma_H_N-alpha_H_N,
         beta_H_biomass = gamma_H_biomass-alpha_H_biomass,
         beta_div_N = gamma_div_N/alpha_div_N,
         beta_div_biomass = gamma_div_biomass/alpha_div_biomass) %>% 
  select(-run) %>% 
  group_by(scenario, evolution, time) %>% 
  summarise(across(.fns = ~mean(.))) %>% 
  ungroup()

# analysis shift from origin


##### plots #####

data %>%
  filter(run == 1) %>%
  filter(evolution == "evolution",
         scenario == "noCC") %>% 
  ggplot(aes(time, environment)) +
  geom_line(aes(linetype = Y), size = 1) +
  scale_color_discrete_qualitative() +
  facet_wrap( ~ X, labeller = "label_both", nrow = 1) +
  guides(color = FALSE, linetype = FALSE)


data_div_alpha %>% 
  ggplot(aes(time, alpha_div_N)) +
  geom_line(aes(color = run)) +
  facet_grid(scenario ~ evolution)

data_div_gamma %>% 
  ggplot(aes(time, gamma_div_N)) +
  geom_line(aes(color = run)) +
  facet_grid(scenario ~ evolution)

data_div %>% 
  select(scenario, evolution, time, contains("div_N")) %>% 
  pivot_longer(cols = contains("div_N"),
               names_to = "level",
               values_to = "div_N") %>% 
  mutate(level = str_remove(level, "_div_N")) %>% 
  ggplot(aes(time, div_N, color = level)) +
  geom_line(aes(linetype = level), size = 1) +
  # geom_vline(aes(xintercept = 2500), size = 1, linetype = 2, color = "grey") +
  # geom_vline(aes(xintercept = 3500), size = 1, linetype = 2, color = "grey") +
  # geom_vline(aes(xintercept = 2000), size = 1, linetype = 2, color = "grey") +
  # geom_vline(aes(xintercept = 4000), size = 1, linetype = 2, color = "grey") +
  facet_grid(scenario ~ evolution)

data_div %>% 
  select(scenario, evolution, time, contains("div_N")) %>% 
  filter(time == max(time)) %>% 
  pivot_longer(cols = contains("div_N"),
               names_to = "level",
               values_to = "div_N") %>% 
  mutate(level = str_remove(level, "_div_N")) %>% 
  ggplot(aes(scenario, div_N, fill = scenario)) +
  geom_col() +
  facet_grid(evolution ~ level) +
  scale_x_discrete(labels = NULL)

data_div %>% 
  select(scenario, evolution, time, contains("div_N")) %>% 
  filter(time == max(time)) %>% 
  pivot_longer(cols = contains("div_N"),
               names_to = "level",
               values_to = "div_N") %>% 
  mutate(level = str_remove(level, "_div_N")) %>% 
  ggplot(aes(evolution, div_N, fill = evolution)) +
  geom_col() +
  facet_grid(level ~ scenario) +
  scale_x_discrete(labels = NULL)

