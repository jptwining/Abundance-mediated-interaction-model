library(nimble)
library(coda)


#########DATA SIMULATOR FOR OCCUPANCY - ABUNDANCE MODEL#########################

#########SIMULATING HIGH ABUNDANCE DOMINANT SPECIES, AND LOW OCCURENCE SUBORDINATE SPECIES ############## #Function definition with set of default values
# data.fn <- function(M = 100, J = 4, mean.lambda.D = 1.5, beta1.D = 1, N.days = 7,
#                     beta2.D = -1, beta3.D = 1, mean.individual.detection.D = 0.05, alpha1.D = -1, alpha2.D = 0.15,
#                     mean.psi.s = 0.5, beta1.s = 1, beta2.s = -1, beta3.s = 1, mean.detection.s = 0.05,
#                     alpha1.s = -1, alpha2.s = 0.15, gamma0 = -1, gamma1 = 0){

data.fn <- function(M = 600, J = 10, mean.lambda.D = 1, beta1.D = 1, N.days = 7,
                    beta2.D = -1, beta3.D = 1, mean.individual.detection.D = 0.5, alpha1.D = -1, alpha2.D = 0.15,
                    mean.psi.s = 0.5, beta1.s = 1, beta2.s = -1, beta3.s = 1, mean.detection.s = 0.5,
                    alpha1.s = -1, alpha2.s = 0.15, gamma0 = -1, gamma1 = 1){
  
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
  forest <- runif(n = M, -1, 1)                       #Scale forest
  road <- runif(n = M, -1, 1)                         #Scale road
  agr <- runif(n = M, -1, 1)                          #Scale agriculture
  day <- matrix(runif(M*J, -1, 1), ncol = J, nrow = M)  #scaled julian day
  
  rep.row<-function(x,n){
    matrix(rep(x,each=n),nrow=n)
  }
  
  
  occ<- rep.row(1:J-mean(1:J), M)  
  str(occ)
  # Model for DOMINANT sp. abundance
  beta0.D <- log(mean.lambda.D)               # Mean abundance on link scale
  lambda.D <- exp(beta0.D + beta1.D*forest + beta2.D*road + beta3.D*agr)
  N.D <- rpois(n = M, lambda = lambda.D)      # Realised abundance
  Ntotal.D <- sum(N.D)                        # Total abundance (all sites)
  psi.true.D <- mean(N.D>0)                   # True occupancy in sample
  
  # Model for  DOMINANT sp. observations
  alpha0.D <- qlogis(mean.individual.detection.D)        # mean detection on link scale
  r.D <- plogis(alpha0.D + alpha1.D*day + alpha2.D*occ)
  p.D <- 1-((1-r.D) ^ N.D) 
  C.D <- matrix(NA, nrow = M, ncol = J)     # Prepare matrix for counts
  for (i in 1:J){                         # Generate counts by survey
    C.D[,i] <- rbinom(n = M, size = N.days, prob = p.D[,i])
  }
  summaxC.D <- sum(apply(C.D,1,max))          # Sum of max counts (all sites)
  psi.obs.D <- mean(apply(C.D,1,max)>0)       # Observed occupancy in sample
  
  #  Model for SUBORDINATE sp. abundance
  beta0.s <- qlogis(mean.psi.s)               # Mean occupancy on logit scale
  psi.s <- plogis(beta0.s + beta1.s*forest + beta2.s*road + beta3.s*agr+gamma0*N.D+gamma1*forest*N.D)
  Z.s <- rbinom(M, 1, psi.s)
  Z.total.s <- sum(Z.s)                        # Total number of occupied sites (all sites where z = 1)
  
  # Model for  SUBORDINATE sp. observations
  alpha0.s <- qlogis(mean.detection.s)        # mean detection on link scale
  p.s <- plogis(alpha0.s + alpha1.s*day + alpha2.s*occ)
  C.s <- matrix(NA, nrow = M, ncol = J)     # Prepare matrix for counts
  for (i in 1:J){                         # Generate counts by survey
    C.s[,i] <- rbinom(n = M, size = Z.s, prob = p.s[,i])
  }
  sum.C.s <- sum(C.s)
  psi.obs.s <- mean(apply(C.s,1,max)>0)   
  
  # Output
  return(list(M = M, J = J, mean.lambda.D = mean.lambda.D, beta0.D = beta0.D, beta1.D = beta1.D,
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

sum(data$Z.s)
sum(data$N.D)
sum(data$N.D>0)


#pull out count data for dominant and subordinate species
PA_highabundancesubordinate <- data$C.s
counts_highabundancedominant <- data$C.D

yD <- data$C.D
ys <- data$C.s


#sampling parameters
J = 600
K = 10

#structure the data

yD=array(yD,dim=c(J,K)) #new data 
yS=array(ys,dim=c(J,K)) 
occ.new=array(data$occ,dim=c(J,K))
day.new=array(data$day,dim=c(J,K))
forest.new=array(data$forest,dim=c(J))
road.new=array(data$road,dim=c(J))
agr.new=array(data$agr,dim=c(J))

#parameters for multinomial obs function
D <- as.integer(7)
# O <- as.integer(1)


#model run parameters
n.iter <- 2500
n.burn <- 2000

# Occupancy-abundance model with interaction parameters on p and multinomial observation model for cameras
#create dataset
str(fdata<- list(yD.new = yD, yS.new = yS, J =J, nsurveys = dim(yD)[2],
                 DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                 forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = D))


#Occu-abu interaction model run
NimModel <- nimbleCode({
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
  for(j in 1:J){
    nD[j] ~ dpois(lambda[j])
    zD[j] <- step(nD[j] - 0.1) #occupancy state for dominant
    log(lambda[j]) <- beta0D + betaD[1] * forestD[j] + betaD[2] * roadD[j] + betaD[3] * agricultureD[j] 
    
    #Observation model for replicated counts
    for(k in 1:nsurveys){
      yD.new[j,k] ~ dbinom(pD[j,k], N.days)
      logit(rD[j,k]) <- alpha0D + alphaD[1] * DATE[j,k] + alphaD[2] * occ[j,k] 
      pD[j,k] <- 1 - ((1-rD[j,k]) ^ nD[j])
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
  for(j in 1:J){
    zS[j] ~ dbern(psiB[j])
    logit(psiB[j]) <- beta0S + betaS[1] * forestS[j] + betaS[2] * roadS[j] + betaS[3] * agricultureS[j] + 
      gamma0 *nD[j] + gamma1*forestS[j]*nD[j]
    #Observation model for replicated cout
    for (k in 1:nsurveys){
      yS.new[j,k] ~ dbern(pS[j,k] * zS[j])
      logit(pS[j,k]) <- alpha0S + alphaS[1] * DATE[j,k] + alphaS[2] * occ[j,k]
    }
  }
  # Derived quantity: Occupancy across all surveyed sites
  # Derived variable
  totalZs <- sum(zS[1:J])
  # Derived quantity: Total abundance across all surveyed sites
  totalnD <- sum(nD[1:J])
}
)



# Parameters monitored
parameters <-  c("alpha0D", "alphaD", "beta0D", "betaD", "alpha0S", "alphaS",
                 "beta0S", "betaS", "totalZs","gamma0", "gamma1", "totalnD")

#data and constants
Nimdata <- list(yD.new = as.matrix(yD), yS.new = as.matrix(yS))
constants=list(J =J, nsurveys = dim(yD)[2],
               DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
               forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = D)

#inits
Niminits=list(beta0D=0, beta0S=0, nD=apply(yD,1,max),zS=apply(yS,1,max), alpha0D = 0, alpha0S = 0, alphaD=c(0,0),betaD=c(0,0,0),alphaS=c(0,0),betaS=c(0,0,0),gamma0=0, gamma1=0)

#thinning rate
nt = 5

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = TRUE, enableWAIC=TRUE)

#correlated
#alpha0D, beta0D. well sometimes.


# Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


# n.iter <- 25000
# Run the model.
start.time2<-Sys.time()
Cmcmc$run(25000,reset=FALSE) #can keep running this line to get more sample
end.time<-Sys.time()
time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2=end.time-start.time2 # post-com


burnin=4000

#get the chains
mvSamples = as.matrix(Cmcmc$mvSamples)
# plot(mcmc(mvSamples[-c(1:burnin),]))

sum(data$Z.s)
sum(data$N.D)
sum(data$N.D>0)


#set up directories 

# path<-getwd()
# dir.create(paste0(path,"/JoshSim_OccuAbu_v3"))
# dir.create(paste0(path,"/JoshSim_OccuAbu_v3/M=1000"))
# dir.create(paste0(path,"/JoshSim_OccuAbu_v3/M=1000/g0=-1_g1=1_K=4"))
# dir.create(paste0(path,"/JoshSim_OccuAbu_v3/M=1000/g0=-1_g1=1_K=4/rerun"))
# 

i <- 250
library(snow)
library(doSNOW)
library(foreach)
cores=5
cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
registerDoSNOW(cl.tmp)

simrep <- 250   # Number of simulation replicates
esti <- array(NA, dim = c(3, (dim(mvSamples)[2])))
colnames(esti) <- colnames(mvSamples)
rownames(esti)<- c("mean", "lowerCI","upperCI")
true.vals <- array(NA, dim = c(18, 1))
rownames(true.vals) <- c("mean.individual.det_D",   "mean.det.s",   "alphaD[1]", "alphaD[2]", "alphaS[1]", "alphaS[2]", "beta0D",    "beta0S",    "betaD[1]",  "betaD[2]",  "betaD[3]",  "betaS[1]",  "betaS[2]",  "betaS[3]", 
                         "gamma0", "gamma1",   "totalZs",   "totalnD")

`%!in%` = Negate(`%in%`)
# Launch simulation
out = foreach(rep=1:simrep) %dopar% {
  # setwd(paste0(path,"/JoshSim_OccuAbu_v3/M=600/g0=-1,g1=1,K=4/rerun"))
  # files = list.files(paste0(path,"/JoshSim_OccuAbu_v3/M=600/g0=-1,g1=1,K=4/rerun"))
  setwd('C:/Users/jpt93/Documents/R/Occu-abu model/simulator_v3/two_species//M=600/g0=-1,g1=1,K=10/rerun')
  files = list.files('C:/Users/jpt93/Documents/R/Occu-abu model/simulator_v3/two_species/M=600/g0=-1,g1=1,K=10/rerun')
  
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
    n.burn <- 45000
    
    
    #sampling parameters
    J = 600
    K = 10
    
    n.chains = 3
    chains = vector("list", n.chains)
    for(chain in 1:n.chains){
      
      
      #pull out count data for dominant and subordinate species
      PA_highabundancesubordinate <- data$C.s
      counts_highabundancedominant <- data$C.D
      
      
      yD <- data$C.D
      ys <- data$C.s
      
      
      #structure the data
      
      yD=array(yD,dim=c(J,K)) #new data 
      yS=array(ys,dim=c(J,K)) 
      occ.new=array(data$occ,dim=c(J,K))
      day.new=array(data$day,dim=c(J,K))
      forest.new=array(data$forest,dim=c(J))
      road.new=array(data$road,dim=c(J))
      agr.new=array(data$agr,dim=c(J))
      
      #parameters for multinomial obs function
      D <- as.integer(7)
      
      
      #setting J for each year
      
      
      # Occupancy-abundance model with interaction parameters on p and multinomial observation model for cameras
      #create dataset
      str(fdata<- list(yD.new = yD, yS.new = yS, J =J, nsurveys = dim(yD)[2],
                       DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                       forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = D))
      
      
      
      #Occu-abu interaction model run
      
      NimModel <- nimbleCode({
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
        for(j in 1:J){
          nD[j] ~ dpois(lambda[j])
          zD[j] <- step(nD[j] - 0.1) #occupancy state for dominant
          log(lambda[j]) <- beta0D + betaD[1] * forestD[j] + betaD[2] * roadD[j] + betaD[3] * agricultureD[j] 
          
          #Observation model for replicated counts
          for(k in 1:nsurveys){
            yD.new[j,k] ~ dbinom(pD[j,k], N.days)
            logit(rD[j,k]) <- alpha0D + alphaD[1] * DATE[j,k] + alphaD[2] * occ[j,k] 
            pD[j,k] <- 1 - ((1-rD[j,k]) ^ nD[j])
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
        for(j in 1:J){
          zS[j] ~ dbern(psiB[j])
          logit(psiB[j]) <- beta0S + betaS[1] * forestS[j] + betaS[2] * roadS[j] + betaS[3] * agricultureS[j] + 
            gamma0 *nD[j] + gamma1*forestS[j]*nD[j]
          #Observation model for replicated cout
          for (k in 1:nsurveys){
            yS.new[j,k] ~ dbern(pS[j,k] * zS[j])
            logit(pS[j,k]) <- alpha0S + alphaS[1] * DATE[j,k] + alphaS[2] * occ[j,k]
          }
        }
        # Derived quantity: Occupancy across all surveyed sites
        # Derived variable
        totalZs <- sum(zS[1:J])
        # Derived quantity: Total abundance across all surveyed sites
        totalnD <- sum(nD[1:J])
      }
      )
      
      #nimble options
      nimbleOptions(clearNimbleFunctionsAfterCompiling = TRUE) 
      
      
      # Parameters monitored
      parameters <-  c("alpha0D", "alphaD", "beta0D", "betaD", "alpha0S", "alphaS",
                       "beta0S", "betaS", "totalZs","gamma0", "gamma1", "totalnD")
      
      #data and constants
      Nimdata <- list(yD.new = as.matrix(yD), yS.new = as.matrix(yS))
      constants=list(J =J, nsurveys = dim(yD)[2],
                     DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                     forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = D)
      #inits
      Niminits=list(beta0D=0,beta0S=0,nD=apply(yD,1,max),zS=apply(yS,1,max), alpha0D = 0, alpha0S = 0, alphaD=c(0,0),betaD=c(0,0,0),alphaS=c(0,0),betaS=c(0,0,0),gamma0=0, gamma1=0)
      
      #thinning rate
      nt = 1
      
      # Build the model, configure the mcmc, and compile
      start.time<-Sys.time()
      Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                            inits=Niminits)
      conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = TRUE, enableWAIC=TRUE)
      
      # Build and compile
      Rmcmc <- buildMCMC(conf)
      Cmodel <- compileNimble(Rmodel)
      Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
      
      # Run the model.
      start.time2<-Sys.time()
      Cmcmc$run(n.iter,reset=FALSE) #can keep running this line to get more sample
      end.time<-Sys.time()
      time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
      time2=end.time-start.time2 # post-com
      
      
      #get the chains
      mvSamples = as.matrix(Cmcmc$mvSamples)
      chains[[chain]]=mvSamples
    }
    
    # setwd(paste0(path,"/JoshSim_OccuAbu_v3/M=1000/g0=-1_g1=1_K=4/rerun"))
    setwd('C:/Users/jpt93/Documents/R/Occu-abu model/simulator_v3/two_species/M=600/g0=-1,g1=1,K=10/rerun')
    
    
    #save the chains
    save(chains, file=paste("Occu_abu_chains_sim",rep,".RData", sep=""))
    
    library(coda)
    # plot(mcmc(mvSamples[n.burn:nrow(mvSamples),]))
    # 
    #combine the chains and burn
    a=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
                mcmc(chains[[2]][n.burn:n.iter,]),
                mcmc(chains[[3]][n.burn:n.iter,]))
    # plot(a)
    
    gelman <- gelman.diag(a)
    
    # 
    a=runjags::combine.mcmc(a)
    
    
    # colnames(mvSamples)
    # [1] "alpha0D"   "alpha0S"   "alphaD[1]" "alphaD[2]" "alphaS[1]" "alphaS[2]" "beta0D"    "beta0S"    "betaD[1]"  "betaD[2]"  "betaD[3]"  "betaS[1]"  "betaS[2]"  "betaS[3]" 
    # [15] "gamma0"    "totalZs"   "totalnD" 
    
    #extract point estimates and credible intervals
    esti[1,1] <- plogis(mean(a[,"alpha0D"]))
    esti[c(2,3),1]<- plogis(quantile(a[,"alpha0D"], probs = c(2.5, 97.5)/100))
    
    esti[1,2] <- plogis(mean(a[,"alpha0S"]))
    esti[c(2,3),2]<- plogis(quantile(a[,"alpha0S"], probs = c(2.5, 97.5)/100))
    
    esti[1,3] <- mean(a[,"alphaD[1]"])
    esti[c(2,3),3]<- quantile(a[,"alphaD[1]"], probs = c(2.5, 97.5)/100)
    
    esti[1,4] <- mean(a[,"alphaD[2]"])
    esti[c(2,3),4]<- quantile(a[,"alphaD[2]"], probs = c(2.5, 97.5)/100)
    
    esti[1,5] <- mean(a[,"alphaS[1]"])
    esti[c(2,3),5]<- quantile(a[,"alphaD[1]"], probs = c(2.5, 97.5)/100)
    
    esti[1,6] <- mean(a[,"alphaS[2]"])
    esti[c(2,3),6]<- quantile(a[,"alphaD[2]"], probs = c(2.5, 97.5)/100)
    
    esti[1,7] <- mean(a[,"beta0D"])
    esti[c(2,3),7]<- quantile(a[,"beta0D"], probs = c(2.5, 97.5)/100)
    
    esti[1,8] <- mean(a[,"beta0S"])
    esti[c(2,3),8]<- quantile(a[,"beta0S"], probs = c(2.5, 97.5)/100)
    
    esti[1,9] <- mean(a[,"betaD[1]"])
    esti[c(2,3),9]<- quantile(a[,"betaD[1]"], probs = c(2.5, 97.5)/100)
    
    esti[1,10] <- mean(a[,"betaD[2]"])
    esti[c(2,3),10]<- quantile(a[,"betaD[2]"], probs = c(2.5, 97.5)/100)
    
    esti[1,11] <- mean(a[,"betaD[3]"])
    esti[c(2,3),11]<- quantile(a[,"betaD[3]"], probs = c(2.5, 97.5)/100)
    
    esti[1,12] <- mean(a[,"betaS[1]"])
    esti[c(2,3),12]<- quantile(a[,"betaD[1]"], probs = c(2.5, 97.5)/100)
    
    esti[1,13] <- mean(a[,"betaS[2]"])
    esti[c(2,3),13]<- quantile(a[,"betaD[2]"], probs = c(2.5, 97.5)/100)
    
    esti[1,14] <- mean(a[,"betaS[3]"])
    esti[c(2,3),14]<- quantile(a[,"betaD[3]"], probs = c(2.5, 97.5)/100)
    
    esti[1,15] <- mean(a[,"gamma0"])
    esti[c(2,3),15]<- quantile(a[,"gamma0"], probs = c(2.5, 97.5)/100)
    
    esti[1,16] <- mean(a[,"gamma1"])
    esti[c(2,3),16]<- quantile(a[,"gamma1"], probs = c(2.5, 97.5)/100)
    
    esti[1,17] <- mean(a[,"totalZs"])
    esti[c(2,3),17]<- quantile(a[,"totalZs"], probs = c(2.5, 97.5)/100)
    
    esti[1,18] <- mean(a[,"totalnD"])
    esti[c(2,3),18]<- quantile(a[,"totalnD"], probs = c(2.5, 97.5)/100)
    
    #put all my shit together and save it.
    allmyshit <- list(esti, true.vals, gelman)
    save(allmyshit, file = paste("Occu_Abu_sim",rep, ".RData", sep=""))
    
    rm(chains)
    gc()
  }}
stopCluster(cl.tmp)
