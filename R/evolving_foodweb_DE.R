

library(tidyverse)
library(deSolve)


trof_model_B <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    I1 <- ((a1 * R)/(1 + a * h1 * R))/bm1
    I2 <- ((a2 * B1)/(1 + a2 * h2 * B1))/bm2
    I3 <- ((a3 * B2)/(1 + a3 * h3 * B2))/bm3
    dR <- In - Out * R - I1 * B1
    dB1 <- I1 * c * B1 - m1 * B1 - I2 * B2
    dB2 <- I2 * g * B2 - m2 * B2 - I3 * B3
    dB3 <- (I3 * g - m3) * B3

    return(list(c(dR, dB1, dB2, dB3)))
  })
}

trof_model1_N <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    I1 <- (a1 * R)/(1 + a * h1 * R)
    I2 <- (a2 * bm1*N1)/(1 + a2 * h2 * bm1*N1)
    I3 <- (a3 * bm2*N2)/(1 + a3 * h3 * bm2*N2)
    dR <- In - Out * R - I1/bm1 * bm1*N1
    dN1 <- (I1/bm1 * c * bm1*N1 - m1 * bm1*N1 - I2/bm2 * bm2*N2)/bm1
    dN2 <- (I2/bm2 * g * bm2*N2 - m2 * bm2*N2 - I3/bm3 * bm3*N3)/bm2
    dN3 <- ((I3/bm3 * g - m3) * bm3*N3)/bm3
    
    return(list(c(dR, dN1, dN2, dN3)))
  })
}

trof_model2_N <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    I1 <- (a1 * R)/(1 + a1 * h1 * R)
    I2 <- (a2 * bm1*N1)/(1 + a2 * h2 * bm1*N1)
    I3 <- (a3 * bm2*N2)/(1 + a3 * h3 * bm2*N2)
    dR <- (In - Out * R - I1 * N1)
    dN1 <- (I1 / bm1 * c * N1 - m1 * N1 - I2 / bm1 * N2)
    dN2 <- (I2 / bm2 * g * N2 - m2 * N2 - I3 / bm2 * N3)
    dN3 <- (I3 / bm3 * g * N3 - m3 * N3)
    
    return(list(c(dR, dN1, dN2, dN3)))
  })
}

trof_model3_N <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    I1 <- (a1 * R)/(1 + a1 * h1 * R)
    I2 <- (a2 * N1)/(1 + a2 * h2 * N1)
    I3 <- (a3 * N2)/(1 + a3 * h3 * N2)
    dR <- (In - Out * R - I1 * N1)
    dN1 <- (I1 * g1 * N1 - m1 * N1 - I2 * N2)
    dN2 <- (I2 * g2 * N2 - m2 * N2 - I3 * N3)
    dN3 <- (I3 * g3 * N3 - m3 * N3)
    
    return(list(c(dR, dN1, dN2, dN3)))
  })
}


# body mass
bm1 <- 100
bm2 <- 10
bm3 <- 1

# mortality
x_0 = 0.1
pow_x = -0.25

# ingestion
a <- 1e-4
h <- 0.4
pw <- 0.75
c <- 10
g <- 0.7

In <- 100
Out <- 0.1
c_conv <- 100
scale <- 5
scale_ass <- 0


# parameters: a named vector
parameters_B <- c(In = In*scale, Out = Out,
                  bm1 = bm1, bm2 = bm2, bm3 = bm3,
                  m1 = x_0*bm1^pow_x, m2 = x_0*bm2^pow_x, m3 = x_0*bm3^pow_x,
                  a1 = a*bm1^(pw)/scale, h1 = h*bm1^(-pw),
                  a2 = a*bm2^(pw)/scale, h2 = h*bm2^(-pw),
                  a3 = a*bm3^(pw)/scale, h3 = h*bm3^(-pw),
                  c = c, g = g)

parameters1_N <- c(In = In*scale, Out = Out,
                  bm1 = bm1, bm2 = bm2, bm3 = bm3,
                  m1 = x_0*bm1^pow_x, m2 = x_0*bm2^pow_x, m3 = x_0*bm3^pow_x,
                  a1 = a*bm1^(pw)/scale, h1 = h*bm1^(-pw),
                  a2 = a*bm2^(pw)/scale, h2 = h*bm2^(-pw),
                  a3 = a*bm3^(pw)/scale, h3 = h*bm3^(-pw),
                  c = c, g = g)

parameters2_N <- c(In = In*scale, Out = Out,
                   bm1 = bm1, bm2 = bm2, bm3 = bm3,
                   m1 = x_0*bm1^pow_x, m2 = x_0*bm2^pow_x, m3 = x_0*bm3^pow_x,
                   a1 = a*bm1^(pw)/scale, h1 = h*bm1^(-pw),
                   a2 = a*bm2^(pw)/scale, h2 = h*(bm2^(-pw)),
                   a3 = a*bm3^(pw)/scale, h3 = h*(bm3^(-pw)),
                   c = c, g = g)

parameters3_N <- c(In = In*scale, Out = Out,
                   bm1 = bm1, bm2 = bm2, bm3 = bm3,
                   m1 = x_0*bm1^pow_x, m2 = x_0*bm2^pow_x, m3 = x_0*bm3^pow_x,
                   a1 = a*bm1^(pw)/scale, h1 = h*bm1^(-pw),
                   a2 = a*bm2^(pw)/scale, h2 = h*(bm2^(-pw))*bm1,
                   a3 = a*bm3^(pw)/scale, h3 = h*(bm3^(-pw))*bm2,
                   g1 = c/bm1, g2 = g*(1-scale_ass*bm2^(-pw))*bm1/bm2, g3 = g*(1-scale_ass*bm3^(-pw))*bm2/bm3)

# time sequence
time <- seq(0, 1000, by = 0.1)

# initial condition: a named vector
# state_B <- c(R = In*scale, B1 = 10000*bm1, B2 = 100*bm2, B3 = 100*bm3)
# state1_N <- c(R = In*scale, N1 = 10000, N2 = 100, N3 = 100)
# state2_N <- c(R = In*scale, N1 = 10000, N2 = 100, N3 = 100)
# state3_N <- c(R = In*scale, N1 = 10000, N2 = 100, N3 = 100)
state_B <- c(R = In*3*scale, B1 = 80*bm1*scale, B2 = 200*bm2*scale, B3 = 2000*bm3*scale)
state1_N <- c(R = In*3*scale, N1 = 80*scale, N2 = 200*scale, N3 = 2000*scale)
state2_N <- c(R = In*3*scale, N1 = 80*scale, N2 = 200*scale, N3 = 2000*scale)
state3_N <- c(R = In*3*scale, N1 = 80*scale, N2 = 200*scale, N3 = 2000*scale)


## Integration with 'ode'
out_B <- ode(y = state_B, times = time, func = trof_model_B, parms = parameters_B)
out1_N <- ode(y = state1_N, times = time, func = trof_model1_N, parms = parameters1_N)
out2_N <- ode(y = state2_N, times = time, func = trof_model2_N, parms = parameters2_N)
out3_N <- ode(y = state3_N, times = time, func = trof_model3_N, parms = parameters3_N)

# plot(out)

out_B %>% 
  data.frame() %>% 
  # filter(time < 50) %>%
  ggplot(aes(time)) +
  geom_line(aes(y = R, color = "R"), lwd = 1) +
  geom_line(aes(y = B1, color = "B1"), lwd = 1) +
  geom_line(aes(y = B2, color = "B2"), lwd = 1) +
  geom_line(aes(y = B3, color = "B3"), lwd = 1) +
  scale_y_log10() +
  labs(y = "B")

out_B %>% 
  data.frame() %>% 
  # filter(time < 50) %>%
  ggplot(aes(time)) +
  geom_line(aes(y = R, color = "R"), lwd = 1) +
  geom_line(aes(y = B1/bm1, color = "N1"), lwd = 1) +
  geom_line(aes(y = B2/bm2, color = "N2"), lwd = 1) +
  geom_line(aes(y = B3/bm3, color = "N3"), lwd = 1) +
  scale_y_log10() +
  labs(y = "N")


out1_N %>% 
  data.frame() %>% 
  # filter(time < 50) %>%
  ggplot(aes(time)) +
  geom_line(aes(y = R), linetype = 2, color = "darkgrey", lwd = 1) +
  geom_line(aes(y = N1, color = "N1"), lwd = 1) +
  geom_line(aes(y = N2, color = "N2"), lwd = 1) +
  geom_line(aes(y = N3, color = "N3"), lwd = 1) +
  scale_y_log10() +
  scale_color_discrete_qualitative() +
  labs(y = "N")


out2_N %>% 
  data.frame() %>% 
  # filter(time < 50) %>%
  ggplot(aes(time)) +
  geom_line(aes(y = R, color = "R"), lwd = 1) +
  geom_line(aes(y = N1, color = "N1"), lwd = 1) +
  geom_line(aes(y = N2, color = "N2"), lwd = 1) +
  geom_line(aes(y = N3, color = "N3"), lwd = 1) +
  scale_y_log10() +
  labs(y = "N")


out3_N %>% 
  data.frame() %>% 
  # filter(time < 50) %>%
  ggplot(aes(time)) +
  geom_line(aes(y = R, color = "R"), lwd = 1) +
  geom_line(aes(y = N1, color = "N1"), lwd = 1) +
  geom_line(aes(y = N2, color = "N2"), lwd = 1) +
  geom_line(aes(y = N3, color = "N3"), lwd = 1) +
  scale_y_log10() +
  labs(y = "N")


