library(nimble)
library(coda)


#########DATA SIMULATOR FOR OCCUPANCY - ABUNDANCE MODEL#########################

#########SIMULATING HIGH ABUNDANCE DOMINANT SPECIES, AND LOW OCCURENCE SUBORDINATE SPECIES ############## #Function definition with set of default values
# data.fn <- function(M = 100, J = 4, mean.lambda.D = 1.5, beta1.D = 1, N.days = 7,
#                     beta2.D = -1, beta3.D = 1, mean.individual.detection.D = 0.05, alpha1.D = -1, alpha2.D = 0.15,
#                     mean.psi.s = 0.5, beta1.s = 1, beta2.s = -1, beta3.s = 1, mean.detection.s = 0.05,
#                     alpha1.s = -1, alpha2.s = 0.15, gamma0 = -1, gamma1 = 0){

data.fn <- function(I = 600, J = 4, mean.lambda.D = 0.5, beta1.D = 1, N.days = 7,
                    beta2.D = -1, beta3.D = 1, mean.individual.detection.D = 0.3, alpha1.D = -1, alpha2.D = 0.15,
                    mean.psi.s = 0.75, beta1.s = 1, beta2.s = -1, beta3.s = 1, mean.detection.s = 0.5,
                    alpha1.s = -1, alpha2.s = 0.15, gamma0 = -1, gamma1 = 0){
  
  #
  # Function to simulate point counts replicated at M sites during J occasions.
  # Population closure is assumed for each site.
  # Expected abundance of dominant species may be affected by forest cover (forest),
  # road density (road) and agricultural land use (agr).
  # Expected detection probability may be affected by julian day and occassion of sampling
  # Expected abundance of subordinate species may be affected by forest cover (forest),
  # road density (road), agricultural land use (agr), the counts of the dominant species, and the
  #interaction between counts of dominant species and forest.
  # Function arguments:
  #     M: Number of spatial replicates (sites)
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
  
  
  occ<- rep.row(1:J-mean(1:J), I)  
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
    C.D[,i] <- rbinom(n = I, size = N.days, prob = p.D[,i])
  }
  summaxC.D <- sum(apply(C.D,1,max))          # Sum of max counts (all sites)
  psi.obs.D <- mean(apply(C.D,1,max)>0)       # Observed occupancy in sample
  
  #  Model for SUBORDINATE sp. abundance
  beta0.s <- qlogis(mean.psi.s)               # Mean occupancy on logit scale
  psi.s <- plogis(beta0.s + beta1.s*forest + beta2.s*road + beta3.s*agr+gamma0*N.D+gamma1*forest*N.D)
  Z.s <- rbinom(I, 1, psi.s)
  Z.total.s <- sum(Z.s)                        # Total number of occupied sites (all sites where z = 1)
  
  # Model for  SUBORDINATE sp. observations
  alpha0.s <- qlogis(mean.detection.s)        # mean detection on link scale
  p.s <- plogis(alpha0.s + alpha1.s*day + alpha2.s*occ)
  C.s <- matrix(NA, nrow = I, ncol = J)     # Prepare matrix for counts
  for (i in 1:J){                         # Generate counts by survey
    C.s[,i] <- rbinom(n = I, size = Z.s, prob = p.s[,i])
  }
  sum.C.s <- sum(C.s)
  psi.obs.s <- mean(apply(C.s,1,max)>0)   
  
  # Output
  return(list(I = I, J = J, mean.lambda.D = mean.lambda.D, beta0.D = beta0.D, beta1.D = beta1.D,
              beta2.D = beta2.D, beta3.D = beta3.D, mean.individual.detection.D = mean.individual.detection.D, alpha0.D = alpha0.D,
              alpha1.D = alpha1.D, alpha2.D = alpha2.D, road = road, forest = forest, day = day, occ = occ,
              agr = agr, lambda.D = lambda.D, N.D = N.D, p.D = p.D, C.D = C.D, Ntotal.D = Ntotal.D,
              psi.true.D = psi.true.D, summaxC.D = summaxC.D, psi.obs.D = psi.obs.D, 
              mean.psi.s  = mean.psi.s, beta0.s = beta0.s, beta1.s = beta1.s,
              beta2.s = beta2.s, beta3.s = beta3.s, mean.detection.s = mean.detection.s, alpha0.s = alpha0.s,
              alpha1.s = alpha1.s, alpha2.s = alpha2.s, psi.s = psi.s, gamma0 = gamma0, gamma1 = gamma1,
              Z.s = Z.s, p.s = p.s, C.s = C.s, Z.total.s = Z.total.s))
}



data <- data.fn()


#pull out count data for dominant and subordinate species
yD <- data$C.D
yS <- data$C.s
yD.occ <- 1*(apply(data$C.D, c(1,2),sum)>0)

str(yS)
str
#sampling parameters
I = 600
J = 4

#structure the data

yD=array(yD,dim=c(I,J)) #new data 
yS=array(yS,dim=c(I,J)) 
occ.new=array(data$occ,dim=c(I,J))
day.new=array(data$day,dim=c(I,J))
forest.new=array(data$forest,dim=c(I))
road.new=array(data$road,dim=c(I))
agr.new=array(data$agr,dim=c(I))

#parameters for multinomial obs function
K <- as.integer(7)
# O <- as.integer(1)


#model run parameters
n.iter <- 2500
n.burn <- 2000

# Occupancy-abundance model with interaction parameters on p and multinomial observation model for cameras
#create dataset
str(fdata<- list(yD.new = yD, yS.new = yS, I =I, nsurveys = dim(yD)[2],
                 DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                 forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = K))



occuNimModel <- nimbleCode({
  # --- Priors ---
  #Model for dominant species (D) 
  #Priors for dominant species (D)
  
  
  #detection parameters
  alpha0D ~ dlogis(0,1)
  # alpha0D <- logit(mean.pD)
  # mean.pD ~ dunif(0, 1)
  for(k in 1:2){
    alphaD[k] ~ dnorm(0, 0.1)
  }
  
  #abundance parameters
  beta0D ~ dnorm(0, 0.1)
  for(k in 1:3){
    betaD[k] ~ dnorm(0, 0.1)
  }
  
  # Likelihood
  # State model for true abundance
  for(i in 1:I){
    nD[i] ~ dpois(lambda[i])
    zD[i] <- step(nD[i] - 0.1) #occupancy state for dominant
    log(lambda[i]) <- beta0D + betaD[1] * forestD[i] + betaD[2] * roadD[i] + betaD[3] * agricultureD[i] 
    
    #Observation model for replicated counts
    for(j in 1:nsurveys){
      yD.new[i,j] ~ dbinom(pD[i,j], N.days)
      logit(rD[i,j]) <- alpha0D + alphaD[1] * DATE[i,j] + alphaD[2] * occ[i,j] 
      pD[i,j] <- 1 - ((1-rD[i,j]) ^ nD[i])
    }
  }
  #Model for subordinate species (S)
  #Priors for subordinate species (s)
  
  #detection parameters
  alpha0S ~ dlogis(0,1)
  # alpha0S <- logit(mean.pS)
  # mean.pS ~ dunif(0, 1)
  for(k in 1:2){
    alphaS[k] ~ dnorm(0, 0.1)
  }
  
  #species occupancy parameters
  # beta0S ~ dnorm(0, 0.1)
  beta0S ~ dlogis(0,1)
  for(k in 1:3){
    betaS[k] ~ dnorm(0, 0.1)
  }
  
  #species interaction parameters
  gamma0 ~ dnorm(0, 0.2) 
  gamma1 ~ dnorm(0, 0.2)
  # epsilon0 ~ dnorm(0, 0.2)
  # Likelihood for subordinate species (S)
  # State model for occupancy for subordinate species (S)
  for(i in 1:I){
    zS[i] ~ dbern(psiB[i])
    logit(psiB[i]) <- beta0S + betaS[1] * forestS[i] + betaS[2] * roadS[i] + betaS[3] * agricultureS[i] + 
      gamma0 *zD[i] + gamma1*forestS[i]*zD[i]
    #Observation model for replicated cout
    for (j in 1:nsurveys){
      yS.new[i,j] ~ dbern(pS[i,j] * zS[i])
      logit(pS[i,j]) <- alpha0S + alphaS[1] * DATE[i,j] + alphaS[2] * occ[i,j]
    }
  }
  # Derived quantity: Occupancy across all surveyed sites
  # Derived variable
  totalZs <- sum(zS[1:I])
  # Derived quantity: Total abundance across all surveyed sites
  totalnD <- sum(nD[1:I])
}
)

# Parameters monitored
occuparameters <-  c("alpha0D", "alphaD", "beta0D", "betaD", "alpha0S", "alphaS",
                    "beta0S", "betaS", "totalZs","gamma0", "gamma1", "totalnD")

#data and constants
occuNimdata <- list(yD.new = as.matrix(yD), yS.new = as.matrix(yS))
occuconstants=list(I=I, nsurveys = dim(yD)[2],
                  DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                  forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = D)

#inits
occuNiminits=list(beta0D=0,nD=apply(yD,1,max),zS=apply(yS,1,max),alphaD=c(0,0),betaD=c(0,0,0),alphaS=c(0,0),betaS=c(0,0,0),gamma0=0, gamma1=0)

#thinning rate
nt = 5

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
occuRmodel <- nimbleModel(code=occuNimModel, constants=occuconstants, data=occuNimdata,check=FALSE,
                         inits=occuNiminits)
occuconf <- configureMCMC(occuRmodel,monitors=occuparameters, thin=nt,useConjugacy = TRUE, enableWAIC=TRUE)

#correlated
#alpha0D, beta0D. well sometimes.

# Build and compile
occuRmcmc <- buildMCMC(occuconf)
occuCmodel <- compileNimble(occuRmodel)
occuCmcmc <- compileNimble(occuRmcmc, project = occuRmodel)


# n.iter <- 25000
# Run the model.
start.time2<-Sys.time()
occuCmcmc$run(25000,reset=FALSE) #can keep running this line to get more sample
end.time<-Sys.time()
time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2=end.time-start.time2 # post-com

burnin=4000

#get the chains
occumvSamples = as.matrix(occuCmcmc$mvSamples)
# plot(mcmc(occumvSamples[-c(1:burnin),]))


abuNimModel <- nimbleCode({
  # --- Priors ---
  #Model for dominant species (D) 
  #Priors for dominant species (D)
  
  
  #detection parameters
  alpha0D ~ dlogis(0,1)
  # alpha0D <- logit(mean.pD)
  # mean.pD ~ dunif(0, 1)
  for(k in 1:2){
    alphaD[k] ~ dnorm(0, 0.1)
  }
  
  #abundance parameters
  beta0D ~ dnorm(0, 0.1)
  for(k in 1:3){
    betaD[k] ~ dnorm(0, 0.1)
  }
  
  # Likelihood
  # State model for true abundance
  for(i in 1:I){
    nD[i] ~ dpois(lambda[i])
    log(lambda[i]) <- beta0D + betaD[1] * forestD[i] + betaD[2] * roadD[i] + betaD[3] * agricultureD[i] 
    
    #Observation model for replicated counts
    for(j in 1:nsurveys){
      yD.new[i,j] ~ dbinom(pD[i,j], N.days)
      logit(rD[i,j]) <- alpha0D + alphaD[1] * DATE[i,j] + alphaD[2] * occ[i,j] 
      pD[i,j] <- 1 - ((1-rD[i,j]) ^ nD[i])
    }
  }
  #Model for subordinate species (S)
  #Priors for subordinate species (s)
  
  #detection parameters
  alpha0S ~ dlogis(0,1)
  # alpha0S <- logit(mean.pS)
  # mean.pS ~ dunif(0, 1)
  for(k in 1:2){
    alphaS[k] ~ dnorm(0, 0.1)
  }
  
  #species occupancy parameters
  # beta0S ~ dnorm(0, 0.1)
  beta0S ~ dlogis(0,1)
  for(k in 1:3){
    betaS[k] ~ dnorm(0, 0.1)
  }
  
  #species interaction parameters
  gamma0 ~ dnorm(0, 0.2) 
  gamma1 ~ dnorm(0, 0.2)
  # epsilon0 ~ dnorm(0, 0.2)
  # Likelihood for subordinate species (S)
  # State model for occupancy for subordinate species (S)
  for(i in 1:I){
    zS[i] ~ dbern(psiB[i])
    logit(psiB[i]) <- beta0S + betaS[1] * forestS[i] + betaS[2] * roadS[i] + betaS[3] * agricultureS[i] + 
      gamma0 *nD[i] + gamma1*forestS[i]*nD[i]
    #Observation model for replicated cout
    for (j in 1:nsurveys){
      yS.new[i,j] ~ dbern(pS[i,j] * zS[i])
      logit(pS[i,j]) <- alpha0S + alphaS[1] * DATE[i,j] + alphaS[2] * occ[i,j]
    }
  }
  # Derived quantity: Occupancy across all surveyed sites
  # Derived variable
  totalZs <- sum(zS[1:I])
  # Derived quantity: Total abundance across all surveyed sites
  totalnD <- sum(nD[1:I])
}
)

# Parameters monitored
abuparameters <-  c("alpha0D", "alphaD", "beta0D", "betaD", "alpha0S", "alphaS",
                     "beta0S", "betaS", "totalZs","gamma0", "gamma1", "totalnD")

#data and constants
abuNimdata <- list(yD.new = as.matrix(yD), yS.new = as.matrix(yS))
abuconstants=list(I=I, nsurveys = dim(yD)[2],
                   DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                   forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = D)

#inits
abuNiminits=list(beta0D=0,nD=apply(yD,1,max),zS=apply(yS,1,max),alphaD=c(0,0),betaD=c(0,0,0),alphaS=c(0,0),betaS=c(0,0,0),gamma0=0, gamma1=0)

#thinning rate
nt = 5

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
abuRmodel <- nimbleModel(code=abuNimModel, constants=abuconstants, data=abuNimdata,check=FALSE,
                          inits=abuNiminits)
abuconf <- configureMCMC(abuRmodel,monitors=abuparameters, thin=nt,useConjugacy = TRUE, enableWAIC=TRUE)

#correlated
#alpha0D, beta0D. well sometimes.

# Build and compile
abuRmcmc <- buildMCMC(abuconf)
abuCmodel <- compileNimble(abuRmodel)
abuCmcmc <- compileNimble(abuRmcmc, project = abuRmodel)


# n.iter <- 25000
# Run the model.
start.time2<-Sys.time()
abuCmcmc$run(25000,reset=FALSE) #can keep running this line to get more sample
end.time<-Sys.time()
time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2=end.time-start.time2 # post-com

burnin=4000

#get the chains
abumvSamples = as.matrix(abuCmcmc$mvSamples)
# plot(mcmc(abumvSamples[-c(1:burnin),]))


#set up directories 

# path<-getwd()
# dir.create(paste0(path,"/JoshSim_OccuAbu"))
# dir.create(paste0(path,"/JoshSim_OccuAbu/M=1000"))
# dir.create(paste0(path,"/JoshSim_OccuAbu/M=1000/g0=-1_g1=-1"))

i <- 250
library(snow)
library(doSNOW)
library(foreach)
cores=5
cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
registerDoSNOW(cl.tmp)

simrep <- 250   # Number of simulation replicates
estiabu <- array(NA, dim = c(3, (dim(abumvSamples)[2])))
estioccu <- array(NA, dim = c(3, (dim(occumvSamples)[2])))
colnames(estiabu) <- colnames(abumvSamples)
rownames(estiabu)<- c("mean", "lowerCI","upperCI")
colnames(estioccu) <- colnames(occumvSamples)
rownames(estioccu)<- c("mean", "lowerCI","upperCI")
true.vals <- array(NA, dim = c(18, 1))
rownames(true.vals) <- c("mean.individual.det_D",   "mean.det.s",   "alphaD[1]", "alphaD[2]", "alphaS[1]", "alphaS[2]", "beta0D",    "beta0S",    "betaD[1]",  "betaD[2]",  "betaD[3]",  "betaS[1]",  "betaS[2]",  "betaS[3]", 
                         "gamma0", "gamma1",   "totalZs",   "totalnD")

`%!in%` = Negate(`%in%`)
# Launch simulation
out = foreach(rep=1:simrep) %dopar% {
  # setwd(paste0(path,"/JoshSim_OccuAbu/M=1000/g0=-1_g1=-1"))
  # files = list.files(paste0(path,"/JoshSim_OccuAbu/M=1000/g0=-1_g1=-1"))
  setwd('C:/Users/jpt93/Documents/R/Occu-abu model/simulator_v3/occu-vs-abu/M=600/g0=-1,g0=1,k=4_lambdaD=2')
  files = list.files('C:/Users/jpt93/Documents/R/Occu-abu model/simulator_v3/occu-vs-abu/M=600/g0=-1,g0=1,k=4_lambdaD=2')
  
  if(paste("Occu_Abu_sim",rep, ".RData", sep="") %!in% files){
    library(nimble)
    library(coda)
    cat(paste("\n\n*** Simrep Nr.", i, "***\n\n"))
    
    
    data <- data.fn()
    true.vals[1,] <- data$mean.individual.detection.D
    true.vals[2,] <- data$mean.detection.s
    true.vals[3,] <- data$alpha1.D
    true.vals[4,] <- data$alpha2.D
    true.vals[5,] <- data$alpha1.s
    true.vals[6,] <- data$alpha2.s
    true.vals[7,] <- data$beta0.D
    true.vals[8,] <- data$beta0.s
    true.vals[9,] <- data$beta1.D
    true.vals[10,] <- data$beta2.D
    true.vals[11,] <- data$beta3.D
    true.vals[12,] <- data$beta1.s
    true.vals[13,] <- data$beta2.s
    true.vals[14,] <- data$beta3.s
    true.vals[15,] <- data$gamma0
    true.vals[16,] <- data$gamma1
    true.vals[17,] <- data$Z.total.s
    true.vals[18,] <- data$Ntotal.D
    
    #model run parameters
    n.iter <- 50000
    n.burn <- 25000
    
    #sampling parameters
    I = 600
    J = 4
    
    n.chains = 3
    occuchains = vector("list", n.chains)
    for(chain in 1:n.chains){
      
      
      #pull out count data for dominant and subordinate species
      yD <- data$C.D
      yS <- data$C.s
      yD.occ <- 1*(apply(data$C.D, c(1,2),sum)>0)
      
      
      #structure the data
      
      yD=array(yD,dim=c(I,J)) #new data 
      yS=array(yS,dim=c(I,J)) 
      yD.occ=array(yD.occ, dim=c(I,J))
      occ.new=array(data$occ,dim=c(I,J))
      day.new=array(data$day,dim=c(I,J))
      forest.new=array(data$forest,dim=c(I))
      road.new=array(data$road,dim=c(I))
      agr.new=array(data$agr,dim=c(I))
      
      #parameters for multinomial obs function
      K <- as.integer(7)
      
      
      #Occupancy-mediated interaction model run
      occuNimModel <- nimbleCode({
        # --- Priors ---
        #Model for dominant species (D) 
        #Priors for dominant species (D)
        #detection parameters
        alpha0D ~ dlogis(0,1)
        # alpha0D <- logit(mean.pD)
        # mean.pD ~ dunif(0, 1)
        for(k in 1:2){
          alphaD[k] ~ dnorm(0, 0.1)
        }
        
        #abundance parameters
        beta0D ~ dnorm(0,1)
        for(k in 1:3){
          betaD[k] ~ dnorm(0, 0.1)
        }
        
        # Likelihood
        # State model for true abundance
        for(i in 1:I){
          nD[i] ~ dpois(lambda[i])
          zD[i] <- step(nD[i] - 0.1) #occupancy state for dominant
          log(lambda[i]) <- beta0D + betaD[1] * forestD[i] + betaD[2] * roadD[i] + betaD[3] * agricultureD[i] 
          
          #Observation model for replicated counts
          for(j in 1:nsurveys){
            yD.new[i,j] ~ dbinom(pD[i,j], N.days)
            logit(rD[i,j]) <- alpha0D + alphaD[1] * DATE[i,j] + alphaD[2] * occ[i,j] 
            pD[i,j] <- 1 - ((1-rD[i,j]) ^ nD[i])
          }
        }
        #Model for subordinate species (S)
        #Priors for subordinate species (s)
        
        #detection parameters
        alpha0S ~ dlogis(0,1)
        # alpha0S <- logit(mean.pS)
        # mean.pS ~ dunif(0, 1)
        for(k in 1:2){
          alphaS[k] ~ dnorm(0, 0.1)
        }
        
        #species occupancy parameters
        # beta0S ~ dnorm(0, 0.1)
        beta0S ~ dlogis(0,1)
        for(k in 1:3){
          betaS[k] ~ dnorm(0, 0.1)
        }
        
        #species interaction parameters
        gamma0 ~ dnorm(0, 0.2) 
        gamma1 ~ dnorm(0, 0.2)
        # Likelihood for subordinate species (S)
        # State model for occupancy for subordinate species (S)
        for(i in 1:I){
          zS[i] ~ dbern(psiB[i])
          logit(psiB[i]) <- beta0S + betaS[1] * forestS[i] + betaS[2] * roadS[i] + betaS[3] * agricultureS[i] + 
            gamma0 *zD[i] + gamma1*forestS[i]*zD[i]
          #Observation model for replicated cout
          for (j in 1:nsurveys){
            yS.new[i,j] ~ dbern(pS[i,j] * zS[i])
            logit(pS[i,j]) <- alpha0S + alphaS[1] * DATE[i,j] + alphaS[2] * occ[i,j]
          }
        }
        # Derived quantity: Occupancy across all surveyed sites
        # Derived variable
        totalZs <- sum(zS[1:I])
        # Derived quantity: Total occupancy across all surveyed sites
        totalnD <- sum(nD[1:I])
      }
      )
      
      # Parameters monitored
      occuparameters <-  c("alpha0D", "alphaD", "beta0D", "betaD", "alpha0S", "alphaS",
                           "beta0S", "betaS", "totalZs","gamma0", "gamma1", "totalnD")
      
      #data and constants
      occuNimdata <- list(yD.new = as.matrix(yD), yS.new = as.matrix(yS))
      occuconstants=list(I=I, nsurveys = dim(yD)[2],
                         DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                         forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = K)
      
      #inits
      occuNiminits=list(beta0D=0,nD=apply(yD,1,max),zS=apply(yS,1,max),alphaD=c(0,0),betaD=c(0,0,0),alphaS=c(0,0),betaS=c(0,0,0),gamma0=0, gamma1=0)
      
      #thinning rate
      nt = 2
      
      # Build the model, configure the mcmc, and compile
      start.time<-Sys.time()
      occuRmodel <- nimbleModel(code=occuNimModel, constants=occuconstants, data=occuNimdata,check=FALSE,
                                inits=occuNiminits)
      occuconf <- configureMCMC(occuRmodel,monitors=occuparameters, thin=nt,useConjugacy = TRUE, enableWAIC=TRUE)
      
      
      #correlated
      #alpha0D, beta0D. well sometimes.
      
      
      # Build and compile
      occuRmcmc <- buildMCMC(occuconf)
      occuCmodel <- compileNimble(occuRmodel)
      occuCmcmc <- compileNimble(occuRmcmc, project = occuRmodel)
      
      
      # Run the model.
      start.time2<-Sys.time()
      occuCmcmc$run(n.iter,reset=FALSE) #can keep running this line to get more sample
      end.time<-Sys.time()
      time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
      time2=end.time-start.time2 # post-com
      
      
      #get the chains
      occumvSamples = as.matrix(occuCmcmc$mvSamples)
      occuchains[[chain]]=occumvSamples
    }
    
    
    n.chains = 3
    abuchains = vector("list", n.chains)
    for(chain in 1:n.chains){
      
      
      #pull out count data for dominant and subordinate species
 
      yD <- data$C.D
      yS <- data$C.s

      
      #structure the data
      
      yD=array(yD,dim=c(I,J)) #new data 
      yS=array(yS,dim=c(I,J)) 
      occ.new=array(data$occ,dim=c(I,J))
      day.new=array(data$day,dim=c(I,J))
      forest.new=array(data$forest,dim=c(I))
      road.new=array(data$road,dim=c(I))
      agr.new=array(data$agr,dim=c(I))
      
      #parameters for multinomial obs function
      K <- as.integer(7)
      
      
      abuNimModel <- nimbleCode({
        # --- Priors ---
        #Model for dominant species (D) 
        #Priors for dominant species (D)
        
        
        #detection parameters
        alpha0D ~ dlogis(0,1)
        # alpha0D <- logit(mean.pD)
        # mean.pD ~ dunif(0, 1)
        for(k in 1:2){
          alphaD[k] ~ dnorm(0, 0.1)
        }
        
        #abundance parameters
        beta0D ~ dnorm(0, 0.1)
        for(k in 1:3){
          betaD[k] ~ dnorm(0, 0.1)
        }
        
        # Likelihood
        # State model for true abundance
        for(i in 1:I){
          nD[i] ~ dpois(lambda[i])
          zD[i] <- step(nD[i] - 0.1) #occupancy state for dominant
          log(lambda[i]) <- beta0D + betaD[1] * forestD[i] + betaD[2] * roadD[i] + betaD[3] * agricultureD[i] 
          
          #Observation model for replicated counts
          for(j in 1:nsurveys){
            yD.new[i,j] ~ dbinom(pD[i,j], N.days)
            logit(rD[i,j]) <- alpha0D + alphaD[1] * DATE[i,j] + alphaD[2] * occ[i,j] 
            pD[i,j] <- 1 - ((1-rD[i,j]) ^ nD[i])
          }
        }
        #Model for subordinate species (S)
        #Priors for subordinate species (s)
        
        #detection parameters
        alpha0S ~ dlogis(0,1)
        # alpha0S <- logit(mean.pS)
        # mean.pS ~ dunif(0, 1)
        for(k in 1:2){
          alphaS[k] ~ dnorm(0, 0.1)
        }
        
        #species occupancy parameters
        # beta0S ~ dnorm(0, 0.1)
        beta0S ~ dlogis(0,1)
        for(k in 1:3){
          betaS[k] ~ dnorm(0, 0.1)
        }
        
        #species interaction parameters
        gamma0 ~ dnorm(0, 0.2) 
        gamma1 ~ dnorm(0, 0.2)
        # Likelihood for subordinate species (S)
        # State model for occupancy for subordinate species (S)
        for(i in 1:I){
          zS[i] ~ dbern(psiB[i])
          logit(psiB[i]) <- beta0S + betaS[1] * forestS[i] + betaS[2] * roadS[i] + betaS[3] * agricultureS[i] + 
            gamma0 *nD[i] + gamma1*forestS[i]*nD[i]
          #Observation model for replicated cout
          for (j in 1:nsurveys){
            yS.new[i,j] ~ dbern(pS[i,j] * zS[i])
            logit(pS[i,j]) <- alpha0S + alphaS[1] * DATE[i,j] + alphaS[2] * occ[i,j]
          }
        }
        # Derived quantity: Occupancy across all surveyed sites
        # Derived variable
        totalZs <- sum(zS[1:I])
        # Derived quantity: Total abundance across all surveyed sites
        totalnD <- sum(nD[1:I])
      }
      )
      
      # Parameters monitored
      abuparameters <-  c("alpha0D", "alphaD", "beta0D", "betaD", "alpha0S", "alphaS",
                           "beta0S", "betaS", "totalZs","gamma0", "gamma1", "totalnD")
      
      
      #data and constants
      abuNimdata <- list(yD.new = as.matrix(yD), yS.new = as.matrix(yS))
      abuconstants=list(I=I, nsurveys = dim(yD)[2],
                        DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                        forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = K)
      
      #inits
      abuNiminits=list(beta0D=0,nD=apply(yD,1,max),zS=apply(yS,1,max),alphaD=c(0,0),betaD=c(0,0,0),alphaS=c(0,0),betaS=c(0,0,0),gamma0=0, gamma1=0)
      
      #thinning rate
      nt = 2
      
      # Build the model, configure the mcmc, and compile
      start.time<-Sys.time()
      abuRmodel <- nimbleModel(code=abuNimModel, constants=abuconstants, data=abuNimdata,check=FALSE,
                               inits=abuNiminits)
      abuconf <- configureMCMC(abuRmodel,monitors=abuparameters, thin=nt,useConjugacy = TRUE, enableWAIC=TRUE)
      
      #correlated
      #alpha0D, beta0D. well sometimes.
      
      # Build and compile
      abuRmcmc <- buildMCMC(abuconf)
      abuCmodel <- compileNimble(abuRmodel)
      abuCmcmc <- compileNimble(abuRmcmc, project = abuRmodel)
      
      # Run the model.
      start.time2<-Sys.time()
      abuCmcmc$run(n.iter,reset=FALSE) #can keep running this line to get more sample
      end.time<-Sys.time()
      time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
      time2=end.time-start.time2 # post-com
      
      
      #get the chains
      abumvSamples = as.matrix(abuCmcmc$mvSamples)
      abuchains[[chain]]=abumvSamples
    }
    
    library(coda)
    # 
    
    str(abuchains)
    str(occuchains)
    
    n.iter = 25000
    n.burn = 20000
    #combine the chains and burn
    abua=mcmc.list(mcmc(abuchains[[1]][n.burn:n.iter,]),
                mcmc(abuchains[[2]][n.burn:n.iter,]),
                mcmc(abuchains[[3]][n.burn:n.iter,]))
    occua=mcmc.list(mcmc(occuchains[[1]][n.burn:n.iter,]),
                   mcmc(occuchains[[2]][n.burn:n.iter,]),
                   mcmc(occuchains[[3]][n.burn:n.iter,]))
    plot(occua)
    plot(abua)
    
    occugelman <- gelman.diag(occua)
    abugelman <- gelman.diag(abua)
    
    
    #combine chains
    occua=runjags::combine.mcmc(occua)
    abua=runjags::combine.mcmc(abua)
    
    # colnames(mvSamples)
    # [1] "alpha0D"   "alpha0S"   "alphaD[1]" "alphaD[2]" "alphaS[1]" "alphaS[2]" "beta0D"    "beta0S"    "betaD[1]"  "betaD[2]"  "betaD[3]"  "betaS[1]"  "betaS[2]"  "betaS[3]" 
    # [15] "gamma0"    "totalZs"   "totalnD" 
    
    #extract point estimates and credible intervals
    estioccu[1,1] <- plogis(mean(occua[,"alpha0D"]))
    estioccu[c(2,3),1]<- plogis(quantile(occua[,"alpha0D"], probs = c(2.5, 97.5)/100))
    
    estioccu[1,2] <- plogis(mean(occua[,"alpha0S"]))
    estioccu[c(2,3),2]<- plogis(quantile(occua[,"alpha0S"], probs = c(2.5, 97.5)/100))
    
    estioccu[1,3] <- mean(occua[,"alphaD[1]"])
    estioccu[c(2,3),3]<- quantile(occua[,"alphaD[1]"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,4] <- mean(occua[,"alphaD[2]"])
    estioccu[c(2,3),4]<- quantile(occua[,"alphaD[2]"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,5] <- mean(occua[,"alphaS[1]"])
    estioccu[c(2,3),5]<- quantile(occua[,"alphaD[1]"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,6] <- mean(occua[,"alphaS[2]"])
    estioccu[c(2,3),6]<- quantile(occua[,"alphaD[2]"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,7] <- mean(occua[,"beta0D"])
    estioccu[c(2,3),7]<- quantile(occua[,"beta0D"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,8] <- mean(occua[,"beta0S"])
    estioccu[c(2,3),8]<- quantile(occua[,"beta0S"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,9] <- mean(occua[,"betaD[1]"])
    estioccu[c(2,3),9]<- quantile(occua[,"betaD[1]"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,10] <- mean(occua[,"betaD[2]"])
    estioccu[c(2,3),10]<- quantile(occua[,"betaD[2]"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,11] <- mean(occua[,"betaD[3]"])
    estioccu[c(2,3),11]<- quantile(occua[,"betaD[3]"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,12] <- mean(occua[,"betaS[1]"])
    estioccu[c(2,3),12]<- quantile(occua[,"betaD[1]"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,13] <- mean(occua[,"betaS[2]"])
    estioccu[c(2,3),13]<- quantile(occua[,"betaD[2]"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,14] <- mean(occua[,"betaS[3]"])
    estioccu[c(2,3),14]<- quantile(occua[,"betaD[3]"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,15] <- mean(occua[,"gamma0"])
    estioccu[c(2,3),15]<- quantile(occua[,"gamma0"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,16] <- mean(occua[,"gamma1"])
    estioccu[c(2,3),16]<- quantile(occua[,"gamma1"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,17] <- mean(occua[,"totalZs"])
    estioccu[c(2,3),17]<- quantile(occua[,"totalZs"], probs = c(2.5, 97.5)/100)
    
    estioccu[1,18] <- mean(occua[,"totalnD"])
    estioccu[c(2,3),18]<- quantile(occua[,"totalnD"], probs = c(2.5, 97.5)/100)
    
    
    estiabu[1,1] <- plogis(mean(abua[,"alpha0D"]))
    estiabu[c(2,3),1]<- plogis(quantile(abua[,"alpha0D"], probs = c(2.5, 97.5)/100))
    
    estiabu[1,2] <- plogis(mean(abua[,"alpha0S"]))
    estiabu[c(2,3),2]<- plogis(quantile(abua[,"alpha0S"], probs = c(2.5, 97.5)/100))
    
    estiabu[1,3] <- mean(abua[,"alphaD[1]"])
    estiabu[c(2,3),3]<- quantile(abua[,"alphaD[1]"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,4] <- mean(abua[,"alphaD[2]"])
    estiabu[c(2,3),4]<- quantile(abua[,"alphaD[2]"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,5] <- mean(abua[,"alphaS[1]"])
    estiabu[c(2,3),5]<- quantile(abua[,"alphaD[1]"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,6] <- mean(abua[,"alphaS[2]"])
    estiabu[c(2,3),6]<- quantile(abua[,"alphaD[2]"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,7] <- mean(abua[,"beta0D"])
    estiabu[c(2,3),7]<- quantile(abua[,"beta0D"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,8] <- mean(abua[,"beta0S"])
    estiabu[c(2,3),8]<- quantile(abua[,"beta0S"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,9] <- mean(abua[,"betaD[1]"])
    estiabu[c(2,3),9]<- quantile(abua[,"betaD[1]"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,10] <- mean(abua[,"betaD[2]"])
    estiabu[c(2,3),10]<- quantile(abua[,"betaD[2]"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,11] <- mean(abua[,"betaD[3]"])
    estiabu[c(2,3),11]<- quantile(abua[,"betaD[3]"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,12] <- mean(abua[,"betaS[1]"])
    estiabu[c(2,3),12]<- quantile(abua[,"betaD[1]"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,13] <- mean(abua[,"betaS[2]"])
    estiabu[c(2,3),13]<- quantile(abua[,"betaD[2]"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,14] <- mean(abua[,"betaS[3]"])
    estiabu[c(2,3),14]<- quantile(abua[,"betaD[3]"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,15] <- mean(abua[,"gamma0"])
    estiabu[c(2,3),15]<- quantile(abua[,"gamma0"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,16] <- mean(abua[,"gamma1"])
    estiabu[c(2,3),16]<- quantile(abua[,"gamma1"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,17] <- mean(abua[,"totalZs"])
    estiabu[c(2,3),17]<- quantile(abua[,"totalZs"], probs = c(2.5, 97.5)/100)
    
    estiabu[1,18] <- mean(abua[,"totalnD"])
    estiabu[c(2,3),18]<- quantile(abua[,"totalnD"], probs = c(2.5, 97.5)/100)
    
    
    #put all my shit together and save it.
    setwd('C:/Users/jpt93/Documents/R/Occu-abu model/simulator_v3/occu-vs-abu/M=600/g0=-1,g0=1,k=4_lambdaD=2')
    allmyshit <- list(estioccu, estiabu, true.vals, occugelman, abugelman)
    save(allmyshit, file = paste("Occu_Abu_sim",rep, ".RData", sep=""))
    
    rm(chains)
    gc()
  }
  }
stopCluster(cl.tmp)

