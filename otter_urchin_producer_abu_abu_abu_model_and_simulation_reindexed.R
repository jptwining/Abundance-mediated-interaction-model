library(nimble)
library(coda)
#########DATA SIMULATOR FOR OCCUPANCY - ABUNDANCE MODEL#########################

#########SIMULATING HIGH ABUNDANCE DOMINANT SPECIES, AND LOW OCCURENCE SUBORDINATE SPECIES ############## 

#Function definition with set of default values
data.fn <- function(I = 600, J = 4, mean.lambda.D = 0.5, beta1.D = 1, K = 7,
                    beta2.D = -1, beta3.D = 1, mean.individual.detection.D = 0.5, alpha1.D = -1, alpha2.D = 0.15,
                    mean.lambda.I = 0.5, beta1.I = 1,
                    beta2.I = -1, beta3.I = 1, mean.individual.detection.I = 0.5, alpha1.I = -1, alpha2.I = 0.15,
                    mean.lambda.s = 0.5, beta1.s = 1,
                    beta2.s = -1, beta3.s = 1, mean.individual.detection.s = 0.5, alpha1.s = -1, alpha2.s = 0.15, gamma0DI = -1, gamma0IS = -1){
  #
  # Function to simulate point counts replicated at I sites during J occasions.
  # Population closure is assumed for each site.
  # Expected abundance of dominant species may be affected by forest cover (forest),
  # road density (road) and agricultural land use (agr).
  # Expected detection probability may be affected by julian day and occassion of sampling
  # Expected abundance of subordinate species may be affected by forest cover (forest),
  # road density (road), agricultural land use (agr), the counts of the dominant species, and the
  #interaction between counts of dominant species and forest.
  # Function arguments:
  #     I: Number of spatial replicates (sites)
  #     J: Number of temporal replicates (occasions)
  #     mean.lambda.D: Mean abundance of dominant species at value 0 of abundance covariates
  #     beta1.D: Main effect of forest cover on dominant species abundance
  #     beta2.D: Main effect of road density on dominant species abundance
  #     beta3.D: Main effect of agriculture on domiant species abundance
  #     mean.detection.D: Mean detection prob. at value 0 of detection covariates of dominant species
  #     alpha1.D: Main effect of julian day on detection probability of dominant species
  #     alpha2.D: Main effect of occassion number on detection probability of dominant species
  #     mean.lambda.s: Mean abundance of subordinate species at value 0 of abundance covariates
  #     beta1.s: Main effect of forest cover on subordinate species abundance
  #     beta2.s: Main effect of road density on subordinate species abundance
  #     beta3.s: Main effect of agriculture on domiant species abundance
  #     mean.setection.s: Mean detection prob. at value 0 of detection covariates of subordinate species
  #     alpha1.s: Main effect of julian day on detection probability of subordinate species
  #     alpha2.s: Main effect of occassion number on detection probability of subordinate species
  #     gamma0: Main effect of dominant species on abundance of subordinate species
  #     gamma1: Interaction effect of dominant species and forest cover on subordinate species
  #     epsilon1: Main effect of dominant species on detection of subordinate species
  
  # Create covariates
  forest <- runif(n = I, -1, 1)                       #Scale forest
  road <- runif(n = I, -1, 1)                         #Scale road
  agr <- runif(n = I, -1, 1)                          #Scale agriculture
  day <- matrix(runif(I*J, -1, 1), ncol = J, nrow = I)  #scaled julian day
  
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  
  occ<- rep.row(1:4-mean(1:4), I)  
  str(occ)
  # Model for DOMINANT sp. abundance
  beta0.D <- log(mean.lambda.D)               # Mean abundance on link scale
  lambda.D <- exp(beta0.D + beta1.D*forest + beta2.D*road + beta3.D*agr)
  N.D <- rpois(n = I, lambda = lambda.D)      # Realised abundance
  Ntotal.D <- sum(N.D)                        # Total abundance (all sites)
  psi.true.D <- mean(N.D>0)                   # True occupancy in sample
  
  
  # Model for  DOMINANT sp. observations
  alpha0.D <- qlogis(mean.individual.detection.D)        # mean detection on link scale
  r.D <- plogis(alpha0.D + alpha1.D*day + alpha2.D*occ)
  p.D <- 1-((1-r.D) ^ N.D) 
  C.D <- matrix(NA, nrow = I, ncol = J)     # Prepare matrix for counts
  for (i in 1:J){                         # Generate counts by survey
    C.D[,i] <- rbinom(n = I, size = K, prob = p.D[,i])
  }
  summaxC.D <- sum(apply(C.D,1,max))          # Sum of max counts (all sites)
  psi.obs.D <- mean(apply(C.D,1,max)>0)       # Observed occupancy in sample

  
  # Model for INTERMEDIATE sp. abundance
  beta0.I <- log(mean.lambda.I)               # Mean abundance on link scale
  lambda.I <- exp(beta0.I + beta1.I*forest + beta2.I*road + beta3.I*agr + gamma0DI * N.D)
  N.I <- rpois(n = I, lambda = lambda.I)      # Realised abundance
  Ntotal.I <- sum(N.I)                        # Total abundance (all sites)
  psi.true.I <- mean(N.I>0)                   # True occupancy in sample
  
  # Model for  INTERMEDIATE sp. observations
  alpha0.I <- qlogis(mean.individual.detection.I)        # mean detection on link scale
  r.I <- plogis(alpha0.I + alpha1.I*day + alpha2.I*occ)
  p.I <- 1-((1-r.I) ^ N.I) 
  C.I <- matrix(NA, nrow = I, ncol = J)     # Prepare matrix for counts
  for (i in 1:J){                         # Generate counts by survey
    C.I[,i] <- rbinom(n = I, size = K, prob = p.I[,i])
  }
  summaxC.I <- sum(apply(C.I,1,max))          # Sum of max counts (all sites)
  psi.obs.I <- mean(apply(C.I,1,max)>0)       # Observed occupancy in sample
  
  
  #  Model for SUBORDINATE sp. abundance
  beta0.s <- log(mean.lambda.s)               # Mean abundance on link scale
  lambda.s <- exp(beta0.s + beta1.s*forest + beta2.s*road + beta3.s*agr + gamma0IS * N.I)
  N.s <- rpois(n = I, lambda = lambda.s)      # Realised abundance
  Ntotal.s <- sum(N.s)                        # Total abundance (all sites)
  psi.true.s <- mean(N.s>0)                       # Total number of occupied sites (all sites where z = 1)
  
  # Model for  SUBORDINATE sp. observations
  alpha0.s <- qlogis(mean.individual.detection.s)        # mean detection on link scale
  r.s <- plogis(alpha0.s + alpha1.s*day + alpha2.s*occ)
  p.s <- 1-((1-r.s) ^ N.s) 
  C.s <- matrix(NA, nrow = I, ncol = J)     # Prepare matrix for counts
  for (i in 1:J){                         # Generate counts by survey
    C.s[,i] <- rbinom(n = I, size = K, prob = p.s[,i])
  }
  summaxC.s <- sum(apply(C.s,1,max))          # Sum of max counts (all sites)
  psi.obs.s <- mean(apply(C.s,1,max)>0)       # Observed occupancy in sample
  
  # Output
  return(list(I = I, 
              J = J, 
              mean.lambda.D = mean.lambda.D, 
              beta0.D = beta0.D,
              beta1.D = beta1.D,
              beta2.D = beta2.D,
              beta3.D = beta3.D, 
              mean.individual.detection.D = mean.individual.detection.D,
              alpha0.D = alpha0.D,
              alpha1.D = alpha1.D,
              alpha2.D = alpha2.D,
              lambda.D = lambda.D, 
              N.D = N.D,
              p.D = p.D,
              C.D = C.D,
              Ntotal.D = Ntotal.D,
              psi.true.D = psi.true.D, 
              summaxC.D = summaxC.D,
              psi.obs.D = psi.obs.D, 
              mean.lambda.s = mean.lambda.s, 
              beta0.s = beta0.s,
              beta1.s = beta1.s,
              beta2.s = beta2.s,
              beta3.s = beta3.s, 
              mean.individual.detection.s = mean.individual.detection.s,
              alpha0.s = alpha0.s,
              alpha1.s = alpha1.s,
              alpha2.s = alpha2.s,
              lambda.s = lambda.s, 
              N.s = N.s,
              p.s = p.s,
              C.s = C.s,
              Ntotal.s = Ntotal.s,
              psi.true.s = psi.true.s, 
              summaxC.s = summaxC.s,
              psi.obs.s = psi.obs.s, 
              mean.lambda.D = mean.lambda.D, 
              beta0.D = beta0.D,
              beta1.D = beta1.D,
              beta2.D = beta2.D,
              beta3.D = beta3.D, 
              mean.individual.detection.I = mean.individual.detection.I,
              alpha0.I = alpha0.I,
              alpha1.I = alpha1.I,
              alpha2.I = alpha2.I,
              lambda.I = lambda.I, 
              N.I = N.I,
              p.I = p.I,
              C.I = C.I,
              Ntotal.I = Ntotal.I,
              psi.true.I = psi.true.I, 
              summaxC.I = summaxC.I,
              psi.obs.I = psi.obs.I, 
              gamma0DI = gamma0DI, gamma0IS = gamma0IS,
              road = road, forest = forest, day = day, occ = occ, agr = agr))
}

data <- data.fn()

#model run parameters
n.iter <- 50000
n.burn <- 20000

#sampling parameters
I = 600
J = 4

#pull out count data for dominant and subordinate species
yD <- data$C.D
yI <- data$C.I
yS <- data$C.s


#structure the data

yD=array(yD,dim=c(I,J)) #new data 
yI=array(yI, dim=c(I,J)) 
yS=array(yS,dim=c(I,J)) 
occ.new=array(data$occ,dim=c(I,J))
day.new=array(data$day,dim=c(I,J))
forest.new=array(data$forest,dim=c(I))
road.new=array(data$road,dim=c(I))
agr.new=array(data$agr,dim=c(I))

#parameters for multinomial obs function
K <- as.integer(7)

#setting I for each year



#Abu-abu-abu interaction model run
library(nimble)

NimModel <- nimbleCode({
  # --- Priors ---
  #Model for dominant species (D)
  #Priors for dominant species (D)
    beta0D ~ dnorm(0, 0.1) #use fixed effect
    alpha0D ~ dnorm(0, 0.1) #use fixed effects
  for (k in 1:2) {
    alphaD[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaD[k] ~ dnorm(0, 0.1)
  }
  # Likelihood
  # State model for true abundance of dominant species (D)
    for (i in 1:I) {
      nD[i] ~ dpois(lambdaD[i])
      log(lambdaD[i]) <-  beta0D + betaD[1] * forestD[i] + betaD[2] * roadD[i] + betaD[3] * agricultureD[i] 
    }
  #Observation model for replicated counts
    for (i in 1:I) {
      for (j in 1:nsurveys){
        logit(pD[i,j]) <- alpha0D + alphaD[1] * DATE[i,j] + alphaD[2] * occ[i,j]
        pD.cum[i,j] <- 1 - (1 - pD[i,j])^nD[i] #prob detecting at least 1 ind per day
        yD.new[i,j] ~ dbinom(pD.cum[i,j],size=K)
      }
    }
  # --- Priors ---
  #Model for dominant species (I)
  #Priors for intermediate species (I)
    beta0I ~ dnorm(0, 0.1) #use fixed effect
    alpha0I ~ dnorm(0, 0.1) #use fixed effects

  for (k in 1:2) {
    alphaI[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaI[k] ~ dnorm(0, 0.1)
  }
  gamma0DI ~ dnorm(0, 0.2)
  
  # Likelihood
  # State model for true abundance of dominant species (I)
  for (i in 1:I) {
    nI[i] ~ dpois(lambdaI[i])
    log(lambdaI[i]) <-  beta0I + betaI[1] * forestI[i] + betaI[2] * roadI[i] + betaI[3] * agricultureI[i] + gamma0DI * nD[i]
  }
  
  #Observation model for replicated counts
  for (i in 1:I) {
    for (j in 1:nsurveys){
      logit(pI[i,j]) <- alpha0I + alphaI[1] * DATE[i,j] + alphaI[2] * occ[i,j]
      pI.cum[i,j] <- 1 - (1 - pI[i,j])^nI[i] #prob detecting at least 1 ind per day
      yI.new[i,j] ~ dbinom(pI.cum[i,j],size=K)
    }
  }
  #Model for subordinate species (S)
  #Priors for subordinate species (s)
    beta0S ~ dlogis(0,1) #use fixed effects
    alpha0S ~ dlogis(0,1) #use fixed effects
  for (k in 1:2) {
    alphaS[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaS[k] ~ dnorm(0, 0.1)
  }
  gamma0IS ~ dnorm(0, 0.2)
  
  # # Likelihood for subordinate species (S)
  # State model for occupancy for subordinate species (S)
  for (i in 1:I) {
    nS[i] ~ dpois(lambdaS[i])
    log(lambdaS[i]) <-  beta0S + betaS[1] * forestS[i] + betaS[2] * roadS[i] + betaS[3] * agricultureS[i] + gamma0IS * nI[i]
  }
  
  #Observation model for replicated counts
  for (i in 1:I) {
    for (j in 1:nsurveys){
      logit(pS[i,j]) <- alpha0S + alphaS[1] * DATE[i,j] + alphaS[2] * occ[i,j]
      pS.cum[i,j] <- 1 - (1 - pS[i,j])^nS[i] #prob detecting at least 1 ind per day
      yS.new[i,j] ~ dbinom(pS.cum[i,j],size=K)
    }
  }
} )

# Parameters monitored
parameters <-  c("alpha0D", "alphaD", "beta0D", "betaD","alpha0I", "alphaI", "beta0I", "betaI", "alpha0S", "alphaS", "beta0S", "betaS", "gamma0DI", "gamma0IS")

#data and constants
Nimdata <- list(yD.new = as.matrix(yD), yS.new = as.matrix(yS), yI.new = as.matrix(yI))
constants=list(I =I, nsurveys = dim(yD)[2], 
               DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
               forestI = forest.new, roadI = road.new, agricultureI = agr.new,
               forestS = forest.new, roadS = road.new, agricultureS = agr.new, K = K)

#inits
nD.init=1*(apply(Nimdata$yD,c(1),sum)>0)+1
nI.init=1*(apply(Nimdata$yI,c(1),sum)>0)+1
nS.init=1*(apply(Nimdata$yS,c(1),sum)>0)+1


Niminits=list(beta0D=0,nD=nD.init,nS=nS.init,nI=nI.init,alphaD=c(0,0),betaD=c(0,0,0), beta0I=0,alphaI=c(0,0),betaI=c(0,0,0),alphaS=c(0,0),betaS=c(0,0,0),gamma0DI=0, gamma0IS=0)

#thinning rate
nt = 2

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = TRUE, enableWAIC=TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

# n.iter <- 100
# Run the model.
start.time2<-Sys.time()
Cmcmc$run(n.iter,reset=FALSE) #can keep running this line to get more sample

end.time<-Sys.time()
time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2=end.time-start.time2 # post-com

n.iter <- 25000
#get the chains
mvSamples = as.matrix(Cmcmc$mvSamples)
plot(mcmc(mvSamples[-c(1:n.burn),]))

