##################################################################################
#                                                                                #
#   PROJECT: Opportunities to interrupt transmission of enteropathogens of       #
#            poultry origin in Maputo, Mozambique: a transmission model analysis #                                                #
#                                                                                #
#      DATE: May 2023                                                            #
#  CODED BY: Kayoko Shioda, PhD, DVM, MPH                                        #
#                                                                                #
##################################################################################

#------------------------------------------------------------------------------#
# Description of this R script
#------------------------------------------------------------------------------#

# Objective:
# To conduct sampling-importance resampling for the transmission dynamic model 
# for Salmonella that has 6 pathways (person-to-person, water, soil, food, 
# live-animals, all other)

#------------------------------------------------------------------------------#
# Load packages and functions and other settings
#------------------------------------------------------------------------------#

# Packages
library(lhs)
library(deSolve) 
library(ggplot2) 
library(SobolSequence)

# Differential equations
SIWFCEOS_model <- function(t,       # Time 
                           x,       # Vector containing the initial condition
                           parms) { # Vector storing transmission rates, etc.
  
  # Initial condition
  S   <- x[1] 
  Ihh <- x[2]  # Number of individuals infected by other infected individuals
  Ihw <- x[3]  # Number of individuals infected by water
  Ihf <- x[4]  # Number of individuals infected by food
  Ihc <- x[5]  # Number of individuals infected by live chickens
  Ihe <- x[6]  # Number of individuals infected by soil
  Iho <- x[7]  # Number of individuals infected by other sources
  W   <- x[8]  # Pathogens in water
  F   <- x[9]  # Pathogens in food
  C   <- x[10] # Prevalence among live chickens
  E   <- x[11] # Pathogens in soil
  O   <- x[12] # Pathogens in "other:
  N   <- S + (Ihh + Ihw + Ihf + Ihc + Ihe + Iho) 
  
  # Define each rate moving from one compartment to the other
  gamma  <- parms[7]  # Recovery rate 
  alpha  <- parms[8]  # Shedding rate into water and soil
  xi     <- parms[9]  # Clearance rate from water and soil
  mu     <- parms[10] # Birth and death rate
  
  # Scale transmission rates by gamma+mu 
  beta_I <- parms[1]*(gamma+mu)/(Ihh+Ihw+Ihf+Ihc+Ihe+Iho) # Transmission rate from person to person
  beta_W <- parms[2]*(gamma+mu)/W # Transmission rate from water to person
  beta_F <- parms[3]*(gamma+mu)/F # Transmission rate from food to person
  beta_C <- parms[4]*(gamma+mu)/C # Transmission rate from live chicken to person 
  beta_E <- parms[5]*(gamma+mu)/E # Transmission rate from soil to person 
  beta_O <- parms[6]*(gamma+mu)/O # Transmission rate from all other sources
  
  dxdt <- numeric(length(x))
  dxdt[1] <- mu*N -beta_I*S*(Ihh+Ihw+Ihf+Ihc+Ihe+Iho) - beta_W*S*W - beta_F*S*F - beta_C*S*C - beta_E*S*E - beta_O*S*O - mu*S + gamma*(Ihh+Ihw+Ihf+Ihc+Ihe+Iho) # dS/dt
  dxdt[2] <-       beta_I*S*(Ihh+Ihw+Ihf+Ihc+Ihe+Iho) - gamma*Ihh - mu*Ihh # dIhh/dt (person-to-person transmission)
  dxdt[3] <-       beta_W*S*W                         - gamma*Ihw - mu*Ihw # dIhw/dt (water-to-person transmission)
  dxdt[4] <-       beta_F*S*F                         - gamma*Ihf - mu*Ihf # dIhf/dt (food-to-person transmission)
  dxdt[5] <-       beta_C*S*C                         - gamma*Ihc - mu*Ihc # dIhc/dt (live chicken-to-person transmission)
  dxdt[6] <-       beta_E*S*E                         - gamma*Ihe - mu*Ihe # dIhe/dt (soil-to-person transmission)
  dxdt[7] <-       beta_O*S*O                         - gamma*Iho - mu*Iho # dIho/dt (all other-to-person transmission) 
  dxdt[8] <-       alpha*(Ihh+Ihw+Ihf+Ihc+Ihe+Iho) - xi*W # dW/dt (Pathogen shedding into water and clearance from water)
  dxdt[9] <- 0 # dF/dt. Constant. Not time varying.
  dxdt[10]<- 0 # dC/dt. Constant. Not time varying.
  dxdt[11]<-       alpha*(Ihh+Ihw+Ihf+Ihc+Ihe+Iho) - xi*E # dE/dt (Pathogen shedding into soil and clearance from soil)
  dxdt[12]<- 0 # dO/dt. Constant. Not time varying.
  
  list(dxdt)
}



#------------------------------------------------------------------------------#
# Sobol sampling for the prevalence and contribution of each pathway
#------------------------------------------------------------------------------#

# Number of parameter sets you want to run
numsamples <- 2.5e+06 

# Number of parameters you want to estiamte
numparam <- 7 # 6 for the Campylobacter model and 7 for the Salmonella model

# Sobol sampling
p <- sobolSequence.points(numparam, count=numsamples)

# Lower and upper bounds for the Latin Hypercubic sampling
lhs_low  <-  rep(0, numparam) 
lhs_high <-  c(0.4, rep(1, numparam-1)) 

# Scale Sobol samples, using the lower and upper bounds defined above
rho <- t(lhs_low + (lhs_high - lhs_low)*t(p)) # Scaled parameter sets

# Need to have these scaled sampled attributable fractions sum to 1. 
for (i in 1:numsamples){
  rho[i, 2:numparam] <- rho[i, 2:numparam]/sum(rho[i,2:numparam]) 
}

#------------------------------------------------------------------------------#
# Set up parameters 
#------------------------------------------------------------------------------#

# Population size in the simulated population
N <- 1 # Should be 1, as we are modeling the prevalence

# Fixed parameters
gamma <- 1  # Recovery rate
alpha <- 1  # Shedding rate into water
xi    <- 1  # Clearance rate from water (Stays in the environment for 1/xi days)
b     <- 70/1000/365.25 * N  
mu    <- 70/1000/365.25      # Background death rate for children <5 yo (70 per 1,000 population during a year at midyear of 2019)

# Set the initial condition 
# (i.e., fraction of individuals and pathogen concentration in each status at Time 1)
S   <- N - 0.1  
Ihh <- 0
Ihw <- 0
Ihf <- 0
Ihc <- 0
Ihe <- 0
Iho <- 0.1 
W   <- 0.1 # Pathogen concentration in water, W_t 
F   <- 1   # Pathogen concentration in food, F_t is constant over time (set to 1 so that it is basically included in beta_F)
C   <- 1   # Prevalence of infected live chicken which is constant over time (set to 1 so that it is basically included in beta_C)
E   <- 0.1 # Pathogen concentration in soil, E_t 
O   <- 1   # Pathogen concentration in all other sources is constant over time (set to 1 so that it is basically included in beta_O)
set.ODEtime <- seq(from=0, to=365.25*2, by=1) # Run the model for two years

# Create a vector containing the initial condition 
set.x0        <- c(S[1], Ihh[1], Ihw[1], Ihf[1], Ihc[1], Ihe[1],  Iho[1], W[1], F[1], C[1], E[1], O[1])
names(set.x0) <- c('S', 'Ihh',  'Ihw',  'Ihf',  'Ihc',  'Ihe',   'Iho',  'W',  'F',  'C',  'E',  'O') 


#------------------------------------------------------------------------------#
# Run ODE with each parameter set and calculate nLL
#------------------------------------------------------------------------------#

# These are containers that we'll fill in
res.prev  <- matrix(NA, nrow=numsamples, ncol=(numparam)) # estimated prevalence (prev_I, prev_Ihh, prev_Ihw, prev_Ihf, prev_Ihc, prev_Iho) will be stored here
parsample <- matrix(NA, nrow=numsamples, ncol=(numparam)) # scaled Sobol sampled parameter sets (prevalence and AFs) will be stored here 
nllsample <- numeric(numsamples) # netative log-likelihood will be stored here
env.conc  <- matrix(NA, nrow=numsamples, ncol=2) # pathogen concentration in the environmental source (water in this case) in the endemic status

# Create a function to run ODE and calculate negative log-likelihood
ode_func <- function (lambda_gm, # prevalence (all pathways) scaled by gamma*mu 
                      partition) { # AFs
  out <- ode(y = set.x0,              # initial state values
             t = set.ODEtime,         # time step (t=0,1,2,3,...)
             func   = SIWFCEOS_model, # ODE function
             parms  = c(lambda_gm*partition, gamma, alpha, xi, mu), # Sobol samples and fixed parameters
             method = 'vode')         # Method that performs integration (lsode, ode45, vode)
  out <- as.data.frame(out)
  
  # Calculate the proportion of infected children (combining all transmission routes)
  out$I <- out$Ihh + out$Ihw + out$Ihf + out$Ihc + out$Ihe + out$Iho
  
  # Simulated prevalence at the stable stage (equilibrium)
  prev_Ihh <- tail(out[,"Ihh"],1)
  prev_Ihw <- tail(out[,"Ihw"],1)
  prev_Ihf <- tail(out[,"Ihf"],1)
  prev_Ihc <- tail(out[,"Ihc"],1)
  prev_Ihe <- tail(out[,"Ihe"],1)
  prev_Iho <- tail(out[,"Iho"],1)
  prev_I   <- tail(out[,"I"],  1)
  W_end    <- tail(out[,"W"],  1)
  E_end    <- tail(out[,"E"],  1)
  
  # Calculate the negative log-likelihood using the MapSan data and the WHO FERG 
  # attributable fractions (Multinomial distribution)
  NLL = -sum(759*c(0.21*c(0.18,0.10,0.46,0.15,0.01,0.10), 1-0.21) * log(c(prev_Ihh, prev_Ihw, prev_Ihf, prev_Ihc, prev_Ihe, prev_Iho, 1-prev_I)))
  
  # Output
  return(c(NLL, prev_I, prev_Ihh, prev_Ihw, prev_Ihf, prev_Ihc, prev_Ihe, prev_Iho, W_end, E_end))
}

# Run ODE and calculate negative log-likelihood for each sampled parameter set
for(i in 1:numsamples){
  if (i%%1000==0) {
    print(i)
  }
  
  # Scaled Sobol sampled prevalence of infection
  target_prev <- rho[i,1]
  
  # Scale the prevalence at steady state by (gamma + mu)
  lambda_gm <- target_prev/(1-target_prev)
  
  # Run ODE
  temp <- ode_func(lambda_gm, # Prevalence scaled by gamma+mu
                   rho[i,-1]) # Scaled Sobol sampled AFs
  
  # Save outputs
  nllsample[i]  <- temp[1]  # negative log-likelihood
  parsample[i,] <- rho[i,]  # scaled LH sampled parameter sets (prevalence and AFs)
  res.prev[i,]  <- temp[2:(numparam+1)] # estimated prevalence (prev_I, prev_Ihh, prev_Ihw, prev_Ihf, prev_Ihc, prev_Iho)
  env.conc[i,]  <- temp[c(numparam+2,numparam+3)]
}

parsample <- as.data.frame(parsample)
names(parsample) <- c("baseline.prev", "AF_h", "AF_w", "AF_f", "AF_c", "AF_e", "AF_o")
res.prev <- as.data.frame(res.prev)
names(res.prev) <- c("prev_I", "prev_Ihh", "prev_Ihw", "prev_Ihf", "prev_Ihc", "prev_Ihe", "prev_Iho")
env.conc <- as.data.frame(env.conc)
names(env.conc) <- c("W_end", "E_end")

#------------------------------------------------------------------------------#
# Resampling
#------------------------------------------------------------------------------#

# Calculate the relative negative log-likelihood
rel.nllsample <- nllsample - min(nllsample)

# Calculate the probability (weights) based on the relative nLL
prob.nllsample <- exp(-rel.nllsample)

# Set the number of samples you want to resample
num.resample <- 10000

# Resample the rows based on the sample importance, prob.nllsample (with replacement)
resample.ind <- sample(1:length(nllsample), num.resample, replace=T, prob=prob.nllsample) 
resampled.parsample <- parsample[resample.ind,]
resampled.res.prev  <- res.prev[resample.ind,]
resampled.env.conc  <- env.conc[resample.ind,]

