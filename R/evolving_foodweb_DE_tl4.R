

library(tidyverse)
library(deSolve)


trof_model3_N <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    I1 <- (a1 * R)/(1 + a1 * h1 * R)
    I2 <- (a2 * N1)/(1 + a2 * h2 * N1)
    I3 <- (a3 * N2)/(1 + a3 * h3 * N2)
    I4 <- (a4 * N3)/(1 + a4 * h4 * N3)
    dR <- (In - Out * R - I1 * N1)
    dN1 <- (I1 * g1 * N1 - m1 * N1 - I2 * N2)
    dN2 <- (I2 * g2 * N2 - m2 * N2 - I3 * N3)
    dN3 <- (I3 * g3 * N3 - m3 * N3 - I4 * N4)
    dN4 <- (I4 * g4 * N4 - m4 * N4)
    
    return(list(c(dR, dN1, dN2, dN3, dN4)))
  })
}

# params 3 trophic levels
# body mass
bm_offset <- 1
bm_power <- 0.67
bm1 <- bm_offset*10^((1 - 1)*bm_power)
bm2 <- bm_offset*10^((2 - 1)*bm_power)
bm3 <- bm_offset*10^((3 - 1)*bm_power)
bm4 <- bm_offset*10^((4 - 1)*bm_power)

# mortality
x_0 = 0.2 #d
pow_x = -0.5 #d_power #standard -0.25

# ingestion
a <- 1e-4
h <- 0.4
pw <- 0.75 #i_power #standard 0.75
c <- 10 #resource_assimilation
g <- 0.7 #assimilation efficiency #standard 0.7

In <- 300 #in_rate
Out <- 0.1 #out_rate
c_conv <- 1 #resource conversion
scale <- 2 #scale_uptake
scale_ass <- 1 #scale_assim


parameters3_N <- c(In = In*scale, Out = Out,
                   bm1 = bm1, bm2 = bm2, bm3 = bm3, bm4 = bm4,
                   m1 = x_0*bm1^pow_x, m2 = x_0*bm2^pow_x, m3 = x_0*bm3^pow_x, m4 = x_0*bm4^pow_x,
                   a1 = a*bm1^(pw)/scale, h1 = h*bm1^(-pw)*c_conv,
                   a2 = a*bm2^(pw)/scale, h2 = h*(bm2^(-pw))*bm1,
                   a3 = a*bm3^(pw)/scale, h3 = h*(bm3^(-pw))*bm2,
                   a4 = a*bm4^(pw)/scale, h4 = h*(bm4^(-pw))*bm3,
                   g1 = c/(bm1^(1/scale_ass)), g2 = g*bm1/(bm2^(1/scale_ass)), g3 = g*bm2/(bm3^(1/scale_ass)), g4 = g*bm3/(bm4^(1/scale_ass)))


# time sequence
time <- seq(0, 1000, by = 0.1)

# initial condition: a named vector
state3_N <- c(R = In*scale, N1 = 5000, N2 = 1000, N3 = 500, N4 = 100)

## Integration with 'ode'
out3_N <- ode(y = state3_N, times = time, func = trof_model3_N, parms = parameters3_N)

# plot(out)


out3_N %>% 
  data.frame() %>% 
  # filter(time < 50) %>%
  ggplot(aes(time)) +
  geom_line(aes(y = R, color = "R"), lwd = 1) +
  geom_line(aes(y = N1, color = "N1"), lwd = 1) +
  geom_line(aes(y = N2, color = "N2"), lwd = 1) +
  geom_line(aes(y = N3, color = "N3"), lwd = 1) +
  geom_line(aes(y = N4, color = "N4"), lwd = 1) +
  scale_y_log10() +
  labs(y = "N")


