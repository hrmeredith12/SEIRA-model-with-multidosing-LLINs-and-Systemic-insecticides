## SEIRA_singleDose: simulates the intervention of LLIN and
# endectocide applied once.
# improved simulation of LLINs (Hill coef and better concentrations for LC50)
# by Hannah Meredith
# Updated: August 16, 2018

# libraries

# install.packages("pracma")
library(pracma)
# install.packages("deSolve")
library(deSolve)
# install.packages("ggplot2")
library(ggplot2)
# install.packages("readxl")
library(readxl)
library(grid)
library(gridExtra)
library(reshape2)

# baseline (no controls)
baseline <- function(y0, P, m, p_h){
  P = c(P,
        C_n = 0 ,
        C_l = 0 ,
        C_h = 0 ,
        D_death = 10,
        hl_d = .1,
        m = m,
        p_h =p_h)

  step = 1

  t=seq(0, P[13] + P[14], by = step) #duration of simulation = t_ss + treatment period

  base <- rk(y = y0,times = t, func = biting2, parms = P, method = "ode45", atol=1e-10, rtol=1e-10)

  total_infected_h <- base[,3] + base[,5]
  time <- base [,1]
  baseline <- list(time, total_infected_h)

}

# bednets alone
nets <- function(y0, P, m, p_h){
  P = c(P,
        C_n = 1*0.75 ,
        C_l = 0 ,
        C_h = 0 ,
        D_death = 10,
        hl_d = .1,
        dose = 1, # to visualize
        period_D = 365, # to visualize
        m = m,
        p_h = p_h)

  nets <- dosing(y0,P)
 }

# Dosing function

dosing <- function(y0, P) {
  # Initialize time window to achieve steady state and run ODEs
  # Establish overall matrix for collecting timecourses
  t_steadyState <- as.vector(P[13])
  treatment_pd <- as.vector(P[14])
  step <- 1
  times = seq(0, t_steadyState, by = step)
  out <- rk(y = y0, times = times,func = biting2, parms = P, method = "ode45", atol=1e-10, rtol=1e-10)
  overall <- out

  # After system reaches steady state, introduce control methods
  period_N <- P[16]
  net <- P[17]
  dose <- P[28]
  period_D <- P[29]

  # Define new time windows. These start and end times for nets and drug dosing
  # will determine the order of subsequent doses
  t_startD = as.vector(tail(overall[, 1], n = 1))
  t_endD = as.vector(tail(overall[, 1], n = 1) + period_D)
  t_startN = as.vector(tail(overall[, 1], n = 1))
  t_endN = as.vector(tail(overall[, 1], n = 1) + period_N)

  # Update initial conditions for next ODE solver. Add first nets and drug
  y0 <- c(
    E_h = as.vector(tail(overall[, 2], n = 1)),
    I_h = as.vector(tail(overall[, 3], n = 1)),
    R_h = as.vector(tail(overall[, 4], n = 1)),
    A_h = as.vector(tail(overall[, 5], n = 1)),
    E_m = as.vector(tail(overall[, 6], n = 1)),
    I_m = as.vector(tail(overall[, 7], n = 1)),
    D_l = as.vector(tail(overall[, 8], n = 1) + dose),
    D_h = as.vector(tail(overall[, 9], n = 1) + dose),
    N = as.vector(net)
  )
  out <- rk(y = y0,times = seq(t_startD, t_endD, by = step), func = biting2, parms = P, method = "ode45", atol=1e-10, rtol=1e-10)

  #Add new times to overall timecourse matrix
  overall <- rbind(overall, out)

  # Update timers
  t_startD = as.vector(tail(overall[, 1], n = 1))
  t_endD = as.vector(tail(overall[, 1], n = 1) + period_D)
  t_startN = as.vector(t_startN + period_N)
  t_endN = as.vector(t_startN + period_N)

  # Update counters
  doseNum = 1
  netNum = 1
  toggle = 0

  # Run simulation for set treatment period (here, 10 years starting after steady state)
  while (tail(overall[, 1], n = 1) < t_steadyState + treatment_pd & toggle == 0) {
    # Ensures simulation cuts off at 10 years
    if (t_endD >= t_steadyState + treatment_pd) {
      t_endD = t_steadyState + treatment_pd
      t_endN = t_steadyState + treatment_pd
      toggle = 1
    }

    # if Drugs and nets have same start time, add new doses to both control strategies
    if (abs(t_startD - t_startN) < 1){
      y0 <- c(
        E_h = as.vector(tail(overall[, 2], n = 1)),
        I_h = as.vector(tail(overall[, 3], n = 1)),
        R_h = as.vector(tail(overall[, 4], n = 1)),
        A_h = as.vector(tail(overall[, 5], n = 1)),
        E_m = as.vector(tail(overall[, 6], n = 1)),
        I_m = as.vector(tail(overall[, 7], n = 1)),
        D_l = as.vector(tail(overall[, 8], n = 1) + dose),
        D_h = as.vector(tail(overall[, 9], n = 1) + dose),
        N = as.vector(net)
      )
    out <-rk(y = y0, times = seq(t_startD, t_endD, by = step), func = biting2, parms = P, method = "ode45", atol=1e-10, rtol=1e-10)
    overall <- rbind(overall, out)
    t_startD = as.vector(tail(overall[, 1], n = 1))
    t_endD = as.vector(tail(overall[, 1], n = 1) + period_D)
    t_startN = as.vector(t_startN + period_N)
    t_endN = as.vector(t_startN + period_N)
    netNum = netNum + 1
    doseNum = doseNum + 1

    # print(c(t_startD, t_endD, t_startN, t_endN))

  }
  # if the drug's start and end time occur before the start of a new net,
  # then add another dose and update drug timers
  else if (t_startD < t_startN & t_endD < t_startN) {
    y0 <- c(
      E_h = as.vector(tail(overall[, 2], n = 1)),
      I_h = as.vector(tail(overall[, 3], n = 1)),
      R_h = as.vector(tail(overall[, 4], n = 1)),
      A_h = as.vector(tail(overall[, 5], n = 1)),
      E_m = as.vector(tail(overall[, 6], n = 1)),
      I_m = as.vector(tail(overall[, 7], n = 1)),
      D_l = as.vector(tail(overall[, 8], n = 1) + dose),
      D_h = as.vector(tail(overall[, 9], n = 1) + dose),
      N = as.vector(tail(overall[, 10], n = 1))
    )
    out <-rk(y = y0, times = seq(t_startD, t_endD, by = step), func = biting2, parms = P, method = "ode45", atol=1e-10, rtol=1e-10)
    overall <- rbind(overall, out)
    t_startD = as.vector(tail(overall[, 1], n = 1))
    t_endD = as.vector(tail(overall[, 1], n = 1) + period_D)
    doseNum = doseNum + 1

    # print(c(t_startD, t_endD, t_startN, t_endN))
    }
  # if the next drug is dosed before the next net, but the net needs to be dosed
  # during the next drug's time window, add a new dose of drug and pause simulation
  # when new net needs to be added
  else if (t_startD < t_startN & t_endD >= t_startN) {
    y0 <- c(
      E_h = as.vector(tail(overall[, 2], n = 1)),
      I_h = as.vector(tail(overall[, 3], n = 1)),
      R_h = as.vector(tail(overall[, 4], n = 1)),
      A_h = as.vector(tail(overall[, 5], n = 1)),
      E_m = as.vector(tail(overall[, 6], n = 1)),
      I_m = as.vector(tail(overall[, 7], n = 1)),
      D_l = as.vector(tail(overall[, 8], n = 1) + dose),
      D_h = as.vector(tail(overall[, 9], n = 1) + dose),
      N = as.vector(tail(overall[, 10], n = 1))
    )
    out <-rk(y = y0,times = seq(t_startD, t_startN, by = step), func = biting2, parms = P, method ="ode45", atol=1e-10, rtol=1e-10)
    overall <- rbind(overall, out)
    t_startD = as.vector(t_endD)
    t_endD = as.vector(t_endD + period_D)
    t_startN = as.vector(tail(overall[, 1], n = 1))
    t_endN = as.vector(tail(overall[, 1], n = 1) + period_N)
    doseNum = doseNum + 1

    # print(c(t_startD, t_endD, t_startN, t_endN))
  }
  # if the nets need to be replaced before a new dose is added, add a new net and run simulation
  # until new drug needs to be added. Then add new drug and let simulation run until end of that
  # dose

  else if (t_startN < t_startD){
    y0 <- c(
      E_h = as.vector(tail(overall[, 2], n = 1)),
      I_h = as.vector(tail(overall[, 3], n = 1)),
      R_h = as.vector(tail(overall[, 4], n = 1)),
      A_h = as.vector(tail(overall[, 5], n = 1)),
      E_m = as.vector(tail(overall[, 6], n = 1)),
      I_m = as.vector(tail(overall[, 7], n = 1)),
      D_l = as.vector(tail(overall[, 8], n = 1)),
      D_h = as.vector(tail(overall[, 9], n = 1)),
      N = as.vector(net)
    )
    out <- rk(y = y0, times = seq(t_startN, t_startD, by = step), func = biting2, parms = P, method = "ode45", atol=1e-10, rtol=1e-10)
    overall <- rbind(overall, out)
    netNum = netNum + 1
    t_startN = as.vector(t_startN + period_N)
    t_endN = as.vector(t_startN + period_N)

    # print(c(t_startD, t_endD, t_startN, t_endN))

    y0 <- c(
      E_h = as.vector(tail(overall[, 2], n = 1)),
      I_h = as.vector(tail(overall[, 3], n = 1)),
      R_h = as.vector(tail(overall[, 4], n = 1)),
      A_h = as.vector(tail(overall[, 5], n = 1)),
      E_m = as.vector(tail(overall[, 6], n = 1)),
      I_m = as.vector(tail(overall[, 7], n = 1)),
      D_l = as.vector(tail(overall[, 8], n = 1) + dose),
      D_h = as.vector(tail(overall[, 9], n = 1) + dose),
      N = as.vector(tail(overall[, 10], n = 1))
    )
    out <- rk(y = y0, times = seq(t_startD, t_endD, by = step), func = biting2, parms = P, method = "ode45", atol=1e-10, rtol=1e-10)
    overall <- rbind(overall, out)
    doseNum = doseNum + 1
    t_startD = as.vector(tail(overall[, 1], n = 1))
    t_endD = as.vector(tail(overall[, 1], n = 1) + period_D)

    # print(c(t_startD, t_endD, t_startN, t_endN))
  }

  }

  TI_h = overall[, 3] + overall[, 5]
  overall_TI_h <- cbind(overall, TI_h)
  overall_TI_h.df <- as.data.frame(overall_TI_h)
  #
  ggplot(overall_TI_h.df, aes(x = time)) +
    geom_line(aes(y = TI_h, colour = "Total infected")) +
    geom_line(aes(y = E_h, colour = "Exposed")) +
    geom_line(aes(y = I_h, colour = "Infected")) +
    geom_line(aes(y = R_h, colour = "Recovered")) +
    geom_line(aes(y = A_h, colour = "Asymptomatically Infected")) +
    ylab(label = "% Human hosts") +
    xlab(label = "Time") +
    coord_cartesian(xlim = c(0, 5000), ylim = c(0, 0.6))
  #
  ggplot(overall_TI_h.df, aes(x = time)) +
    geom_line(aes(y = D_l, colour = "Drug")) +
    # geom_line(aes(y = N, colour = "Net")) +
    ylab(label = "Drug") +
    xlab(label = "Time") +
    coord_cartesian(xlim = c(0, 5000), ylim = c(0, 10))

  total_infected_h <- overall[,3]+overall[,5]
  time <- overall [,1]
  dosing <- list(time, total_infected_h, overall[,10])
}

# ODEs for malaria transmission with control methods (LLINs and Systemic insecticides)
biting2 <- function(t, y, parms) {
  with(as.list(c(parms, y)), {

    k_n = log(2) / hl_n
    k_d = log(2) / hl_d


    # EQs and ODEs
    b_h <- a * b * (1 - C_n * N^hill_n / (N^hill_n + N_death^hill_n))
    b_m <- a * c * (1 - C_n * N^hill_n / (N^hill_n + N_death^hill_n))
    mu_2c <- mu_2 + 1 / 3 * (p_h * C_n * N^hill_n / (N^hill_n + N_death^hill_n) * mu_n + (1 - p_h) * C_l * D_l^hill_d / (D_l^hill_d + D_death^hill_d) *
          mu_d + (p_h) * C_h * D_h^hill_d / (D_h^hill_d + D_death^hill_d) * mu_d)

    S_h <- 1 - E_h - I_h - R_h - A_h
    S_m <- 1 - E_m - I_m

    dE_h <- m * b_h * p_h * S_h * I_m  - (1 / tau_h + r + mu_1) * E_h
    dI_h <- 1 / tau_h * E_h - (r + q1 + mu_1) * I_h
    dR_h <- q1 * (I_h + A_h) - (theta * b_h * m * p_h * I_m + q2 + mu_1) * R_h
    dA_h <- theta * b_h * m * p_h * I_m * R_h - (q1 + r + mu_1) * A_h

    dE_m <- b_m * p_h * (I_h + sigma * A_h) * S_m - (1 / tau_m + mu_2c) * E_m
    dI_m <- (1 / tau_m) * E_m - mu_2c * I_m
    dD_l <- -k_d * D_l
    dD_h <- -k_d * D_h
    dN <- -k_n * N

    # print(c(t,mu_2c))

    res <- c(dE_h, dI_h, dR_h, dA_h, dE_m, dI_m, dD_l, dD_h, dN)
    list(res)
  })
}

# define parameters and  ranges

# Parameters
P <- c(
  a = 0.2,
  b = 0.5,
  c = 0.5,
  r = 0.01,
  mu_1 = 1 / 21900,
  mu_2 = 0.12,
  tau_m = 10,
  tau_h = 21,
  q1 = 1 / 200,
  q2 = 1 / 1000,
  sigma = 0.25,
  theta = 0.5,
  t_ss = 1000,
  treatment_pd = 10*365,
  N_death = 0.73,
  period_N = 1095,
  net = 2,
  hl_n = 1906,
  mu_d = 0.6,
  mu_n = 0.1,
  hill_n = 2,
  hill_d = 2
)

cmax_range <- logspace(0, 4, n = 10) #n =10
period_D_range <- c(365 / 52, linspace(365 / 12, 365, n = 12)) #n =12
halflife_d_range <- logspace(-1, 2, n = 10) #n =10
D_death_range <- c(7.5, 20, 30, 65, 500, 1175 )#c(10, 25, 100, 500, 1100)
m_range <- c(10, 20, 40)
p_h_range <- c(0.35, 0.5, 0.8)

# define initial conditions
y0 <- c(
  E_h = 0.01,
  I_h = 0.01,
  R_h = 0.01,
  A_h = 0.01,
  E_m = 0.01,
  I_m = 0.01,
  D_l = 0,
  D_h = 0,
  N = 0
)

#target reduction in malaria prevalence (%)
target <- 10

# Simulate time courses, capture malaria prevalence at end of treatment window, and plot contour maps of
# dosing frequencies needed for each combination.
  x<-size(halflife_d_range, 2)
  y<-size(cmax_range, 2)
  z<-size(D_death_range,2)

someData <- rep(0, x*y*z)

periodUse <- array(someData, c( x, y, z))
rel_mal_H <- array(someData, c( x, y, z)) #matrix(,nrow = size(halflife_d_range, 2), ncol = size(cmax_range, 2))

# calculate number of infected cases without any control methods (baseline)
base <- baseline(y0,P,m=m_range[2], p_h=p_h_range[2])
time_base <- base[[1]]
infected_base <- base[[2]]
infected_base_steadyState <- tail(base[[2]],n=1)

# calculate number of infected cases with only bed nets
net <- nets(y0, P, m=m_range[2], p_h = p_h_range[2])
time_net <- net[[1]]
infected_net <- net[[2]]
net_dosing <- net[[3]]
net.df <- data.frame(time_net, infected_net)
infected_net_steadyState <- tail(net[[2]],n=1)
prevalence_ratio_net <- infected_net/infected_base_steadyState # used in plotting time courses

# x1 <- time_net
# y1 <- infected_net

# calculate prevalence with bednets and different systemic insecticide treatments
for (h in 1:z) {
   for (j in 1:x) {
     for (k in 1:y) {
      for (i in size(period_D_range, 2):1) {

        P = c(P[1:22],
          C_n = 1 * 0.75,
          C_l = 1 * 1,
          C_h = 0 * 1,
          D_death = D_death_range[h],
          hl_d = halflife_d_range[j], # 5
          dose = cmax_range[k], #30
          period_D = period_D_range[i], #365 / 12
          m = m_range[2],
          p_h = p_h_range[2])

        # calculate number of infected cases with bednets and systemic insecticides
        # compute the relative prevalence of malaria cases of avoided with addition of SI
        all <- dosing(y0, P)
        time_all <- all[[1]]
        infected_all <- all[[2]]
        net_dosing <- all[[3]]

        # x2 <- time_all
        # y2 <- infected_all
        # p <- plot(x1,y1,ylim=range(c(y1,y2)),xlim=range(c(x1,x2)), type="l",col="black")
        # p + lines(x2,y2,col="red")

        all.df <- data.frame(time_all, infected_all)

        prevalence_ratio_all <- infected_all/infected_base_steadyState # used in plotting time courses
        relative_malaria <- (infected_net_steadyState - tail(infected_all,n=1))/infected_net_steadyState * 100
        relative_malaria [relative_malaria < 0] = 0
        # relative_malaria.df <- as.data.frame(relative_malaria)

        print(c(i,k,j,h,relative_malaria))

        if (relative_malaria >= target){
          periodUse[j,k,h] = period_D_range[i]
          break
        }
      }

    }
  }


}

# eprinomectin and ivermectin share same contour map because they have similar LC50, so duplicate that data for plotting purposes:

test <- array(someData, c( x, y, z+1))
test[,,1:2] <- periodUse[,,1]
test[,,3:7] <- periodUse[,,2:6]

# Plotting contour maps and overlaying data

setwd("Users/Hannah/Desktop/LSHTM/Models/R_modeling/SEIRA")

dataset <- read_excel("20180810_cmax_halflife_LC50_forR.xlsx")
dataset$host <- as.factor(dataset$host)
dataset$drug <- as.factor(dataset$drug)
dataset$route <- as.factor(dataset$route)
dataset$loopID <- ifelse(is.na(dataset$cmaxMean),0,
                         ifelse(is.na(dataset$halflifeMean),0,
                                ifelse(dataset$drug == 'eprinomectin',1, 
                                       ifelse (dataset$drug == 'ivermectin',2,
                                               ifelse (dataset$drug == 'fluralaner',3,
                                                       ifelse (dataset$drug == 'doramectin',4,
                                                               ifelse (dataset$drug == 'afoxalaner',5,
                                                                       ifelse (dataset$drug == 'spinosad',6,
                                                                               ifelse (dataset$drug == 'moxidectin',7,0)))))))))

pltlist <- list()

for (i in 1:7){
  data <- round(test[,,i])#round(periodUse[,,i])
  colnames(data) <- cmax_range
  rownames(data) <- halflife_d_range
  longData <- melt(data)
  colnames(longData) <- c("halflife", "cmax", "period")
  par(pty = "s") # axis
  subdataset <- subset(dataset, dataset$loopID == i)
  
  p <- ggplot(longData, aes(halflife, cmax))
  p1 <- p +
    stat_contour(
      aes(z = period,
          colour = stat(level)),
      size = 1.5,
      breaks = c(
        round(365 / 52),
        round(365 / 12),
        round(365 / 4),
        round(365 / 2),
        round(365)
      )
    ) +
    scale_color_gradient(low = "black", high = "grey", guide = FALSE) +
    scale_x_log10() + scale_y_log10() +
    theme(aspect.ratio = 1) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    ) +
    xlab("Half life (days)") +
    theme(
      axis.text.x = element_text(colour = 'black', size = 10),
      axis.title.x = element_text(
        size = 12,
        hjust = 0.5,
        vjust = 0.5
      )
    ) +
    ylab("C_max (ng/mL)") +
    theme(
      axis.text.y = element_text(colour = 'black', size = 10),
      axis.title.y = element_text(
        size = 12,
        hjust = 0.5,
        vjust = 0.2
      )
    ) +
    # coord_cartesian(xlim = c(0.1, 100), ylim = c(1, 10000)) +
    geom_point(data = subdataset, aes(x=halflifeMean, y=cmaxMean, shape=as.factor(subdataset$route)), color = as.numeric(subdataset$host))
  
  
  pltlist[[i]] <- p1
  
}



grid.arrange(pltlist[[1]],pltlist[[2]],pltlist[[3]],pltlist[[4]],pltlist[[5]],pltlist[[6]],pltlist[[7]],nrow=2,ncol=4)

# order of plots: (1) Eprinomectin, (2) Ivermectin, (3) Fluralaner, (4) Doramectin, (5) Afoxalaner, (6) Spinosad, (7) Moxidectin 

