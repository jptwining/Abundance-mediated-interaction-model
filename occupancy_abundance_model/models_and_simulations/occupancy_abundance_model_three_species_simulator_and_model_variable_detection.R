library(nimble)
library(coda)
#########DATA SIMULATOR FOR OCCUPANCY - ABUNDANCE MODEL#########################

#########SIMULATING HIGH ABUNDANCE DOMINANT SPECIES, AND LOW OCCURENCE SUBORDINATE SPECIES ############## 

#Function definition with set of default values
data.fn <- function(I = 600, J = 4, mean.lambda.D = 0.5, beta1.D = 1, N.days = 7,
                    beta2.D = -1, beta3.D = 1, mean.individual.detection.D = 0.1, alpha1.D = -1, alpha2.D = 0.15,
                    mean.lambda.I = 1, beta1.I = 1, 
                    beta2.I = -1, beta3.I = 1, mean.individual.detection.I = 0.5, alpha1.I = -1, alpha2.I = 0.15,
                    mean.psi.s = 0.25, beta1.s = 1, beta2.s = -1, beta3.s = 1, mean.detection.s = 0.5,
                    alpha1.s = -1, alpha2.s = 0.15, gamma0DI = -1, gamma0IS = -1, gamma0DS = 1){
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
    matrix(rep(x,each=n),nrow=I)
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
    C.D[,i] <- rbinom(n = I, size = N.days, prob = p.D[,i])
  }
  summaxC.D <- sum(apply(C.D,1,max))          # Sum of max counts (all sites)
  psi.obs.D <- mean(apply(C.D,1,max)>0)       # Observed occupancy in sample
  
  
  # Model for INTERMEDIATE sp. abundance
  beta0.I <- log(mean.lambda.I)               # Mean abundance on link scale
  lambda.I <- exp(beta0.I + beta1.I*forest + beta2.I*road + beta3.I*agr + gamma0DI*N.D)
  N.I <- rpois(n = I, lambda = lambda.I)      # Realised abundance
  Ntotal.I <- sum(N.I)                        # Total abundance (all sites)
  psi.true.I <- mean(N.I>0)                   # True occupancy in sample
  
  
  # Model for  INTERMEDIATE sp. observations
  alpha0.I <- qlogis(mean.individual.detection.I)        # mean detection on link scale
  r.I <- plogis(alpha0.I + alpha1.I*day + alpha2.I*occ)
  p.I <- 1-((1-r.I) ^ N.I) 
  C.I <- matrix(NA, nrow = I, ncol = J)     # Prepare matrix for counts
  for (i in 1:J){                         # Generate counts by survey
    C.I[,i] <- rbinom(n = I, size = N.days, prob = p.I[,i])
  }
  summaxC.I <- sum(apply(C.I,1,max))          # Sum of max counts (all sites)
  psi.obs.I <- mean(apply(C.I,1,max)>0)       # Observed occupancy in sample
  
  #  Model for SUBORDINATE sp. abundance
  beta0.s <- qlogis(mean.psi.s)               # Mean occupancy on logit scale
  psi.s <- plogis(beta0.s + beta1.s*forest + beta2.s*road + beta3.s*agr+gamma0IS*N.I +gamma0DS*N.D)
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
              alpha1.D = alpha1.D, alpha2.D = alpha2.D, lambda.D = lambda.D, N.D = N.D, p.D = p.D, C.D = C.D, Ntotal.D = Ntotal.D,
              psi.true.D = psi.true.D, summaxC.D = summaxC.D, psi.obs.D = psi.obs.D, 
              mean.psi.s  = mean.psi.s, beta0.s = beta0.s, beta1.s = beta1.s,
              beta2.s = beta2.s, beta3.s = beta3.s, mean.detection.s = mean.detection.s, alpha0.s = alpha0.s,
              alpha1.s = alpha1.s, alpha2.s = alpha2.s, psi.s = psi.s, 
              mean.lambda.I = mean.lambda.I, beta0.I = beta0.I, beta1.I = beta1.I,
              beta2.I = beta2.I, beta3.I = beta3.I, mean.individual.detection.I = mean.individual.detection.I, alpha0.I = alpha0.I,
              alpha1.I = alpha1.I, alpha2.I = alpha2.I, lambda.I = lambda.I, N.I = N.I, p.I = p.I, C.I = C.I, Ntotal.I = Ntotal.I,
              psi.true.I = psi.true.I, summaxC.I = summaxC.I, psi.obs.I = psi.obs.I,gamma0DI = gamma0DI, gamma0IS = gamma0IS, gamma0DS = gamma0DS,
              Z.s = Z.s, p.s = p.s, C.s = C.s, Z.total.s = Z.total.s, road = road, forest = forest, day = day, occ = occ,
              agr = agr))
}

data <- data.fn()

#model run parameters
n.iter <- 50000
n.burn <- 20000

#sampling parameters
I =  600
J = 4

#pull out count data for dominant and subordinate species
PA_highabundancesubordinate <- data$C.s
counts_highabundancedominant <- data$C.D


yD <- data$C.D
yI <- data$C.I
ys <- data$C.s


#structure the data

yD=array(yD,dim=c(I,J)) #new data 
yI=array(yI, dim=c(I,J)) 
yS=array(ys,dim=c(I,J)) 
occ.new=array(data$occ,dim=c(I,J))
day.new=array(data$day,dim=c(I,J))
forest.new=array(data$forest,dim=c(I))
road.new=array(data$road,dim=c(I))
agr.new=array(data$agr,dim=c(I))

#parameters for multinomial obs function
D <- as.integer(7)
O <- as.integer(1)


#setting J for each year


#Occupancy-abundance model with interaction parameters on p and multinomial observation model for cameras
#create dataset
str(fdata<- list(yD.new = yD, yI.new = yI, yS.new = yS, I = I, nsurveys = dim(yD)[2],
                 DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                 forestI = forest.new, roadI = road.new, agricultureI = agr.new,
                 forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = D, O = O))


#Occu-abu interaction model run

NimModel <- nimbleCode({
  # --- Priors ---
  #Model for dominant species (D) 
  #Priors for dominant species (D)
  alpha0D ~ dlogis(0,1)
  for (k in 1:2){
    alphaD[k] ~ dnorm(0, 0.1)
  }
  beta0D ~ dnorm(0, 0.1)
  for (k in 1:3){
    betaD[k] ~ dnorm(0, 0.1)
  }
  # Likelihood
  # State model for true abundance
  for (i in 1:I){
    nD[i] ~ dpois(lambdaD[i])
    log(lambdaD[i]) <- beta0D + betaD[1] * forestD[i] + betaD[2] * roadD[i] + betaD[3] * agricultureD[i] 
    
    #Observation model for replicated counts
    for (j in 1:nsurveys){
      yD.new[i,j] ~ dbinom(pD[i,j], N.days)
      logit(rD[i,j]) <- alpha0D + alphaD[1] * DATE[i,j] + alphaD[2] * occ[i,j] 
      pD[i,j] <- 1 - ((1-rD[i,j]) ^ nD[i])
    }
  }
  
  #Model for intermediate species (I)
  #Priors for intermediate species (I)
  alpha0I ~ dlogis(0,1)
  for (k in 1:2){
    alphaI[k] ~ dnorm(0, 0.1)
  }
  beta0I ~ dnorm(0, 0.1)
  for (k in 1:3){
    betaI[k] ~ dnorm(0, 0.1)
  }
  gamma0DI ~ dnorm(0, 0.2)
  # Likelihood for intermediate species (I)
  # State model for true abundance for intermediate species (I)
  for (i in 1:I){
    nI[i] ~ dpois(lambdaI[i])
    log(lambdaI[i]) <- beta0I + betaI[1] * forestI[i] + betaI[2] * roadI[i] + betaI[3] * agricultureI[i] + gamma0DI*nD[i]
    
    #Observation model for replicated counts
    for (j in 1:nsurveys){
      yI.new[i,j] ~ dbinom(pI[i,j], N.days)
      logit(rI[i,j]) <- alpha0I + alphaI[1] * DATE[i,j] + alphaI[2] * occ[i,j] 
      pI[i,j] <- 1 - ((1-rI[i,j]) ^ nI[i])
    }
  }
  #Model for subordinate species (S)
  #Priors for subordinate species (s)
  alpha0S ~ dlogis(0,1)
  
  for (k in 1:2){
    alphaS[k] ~ dnorm(0, 0.1)
  }
  beta0S ~ dlogis(0,1)
  for (k in 1:3){
    betaS[k] ~ dnorm(0, 0.1)
  }
  gamma0IS ~ dnorm(0, 0.2)
  gamma0DS ~ dnorm(0, 0.2)
  # gamma1 ~ dnorm(0, 0.2)
  # epsilon0 ~ dnorm(0, 0.2)
  # Likelihood for subordinate species (S)
  # State model for occupancy for subordinate species (S)
  for (i in 1:I){
    zS[i] ~ dbern(psiB[i])
    logit(psiB[i]) <- beta0S + betaS[1] * forestS[i] + betaS[2] * roadS[i] + betaS[3] * agricultureS[i] +gamma0IS *nI[i] +gamma0DS*nD[i]
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
  # Derived quantity: Total abundance across all surveyed sites
  totalnI <- sum(nI[1:I])
  
}
)


# Parameters monitored
parameters <-  c("alpha0D", "alphaD", "beta0D", "betaD","alpha0I", "alphaI", "beta0I", "betaI", "alpha0S", "alphaS", "beta0S", "betaS", "totalZs","totalnI", "totalnD", "gamma0DI", "gamma0IS", "gamma0DS")

#data and constants
Nimdata <- list(yD.new = as.matrix(yD), yS.new = as.matrix(yS), yI.new = as.matrix(yI))
constants=list(I = I, nsurveys = dim(yD)[2],
               DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
               forestI = forest.new, roadI = road.new, agricultureI = agr.new,
               forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = D)

#inits
Niminits=list(beta0D=0,beta0I=0, beta0S = 0, nD=apply(yD,1,max),nI=apply(yI,1,max),zS=apply(yS,1,max),alpha0D=0, alpha0I=0, alpha0S=0, alphaD=c(0,0),alphaI=c(0,0),alphaS=c(0,0), betaD=c(0,0,0), betaI=c(0,0,0),betaS=c(0,0,0),gamma0DI=0, gamma0IS=0, gamma0DS=0)

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

path<-getwd()
dir.create(paste0(path,"/three_species"))
dir.create(paste0(path,"/three_species/M=600"))
dir.create(paste0(path,"/three_species/M=600/rD=0.1,rI=0.5,pS=0.5"))

`%!in%` = Negate(`%in%`)

i <- 5
library(snow)
library(doSNOW)
library(foreach)
cores=5
cl.tmp = makeCluster(rep("localhost",cores), type="SOCK")
registerDoSNOW(cl.tmp)


simrep <- 5  # Number of simulation replicates
esti <- array(NA, dim = c(3, (dim(mvSamples)[2])))
colnames(esti) <- colnames(mvSamples)
rownames(esti)<- c("mean", "lowerCI","upperCI")
true.vals <- array(NA, dim = c(27, 1))
rownames(true.vals) <- c("mean.individual.det_D",  "mean.individual.D_I", "mean.det.s",   "alphaD[1]", "alphaD[2]","alphaI[1]", "alphaI[2]", "alphaS[1]", "alphaS[2]",
                         "beta0D", "beta0I",  "beta0S",    "betaD[1]",  "betaD[2]",  "betaD[3]", "betaI[1]",  "betaI[2]",  "betaI[3]",  "betaS[1]",  "betaS[2]",  "betaS[3]", 
                         "gamma0DI", "gamma0DS", "gamma0IS", "totalZs",   "totalnD", "totalnI")

files = list.files(paste0(path,"/three_species/M=600/rD=0.1,rI=0.5,pS=0.5"))

# Launch simulation
out = foreach(rep=1:simrep) %dopar% {
  setwd(paste0(path,"/three_species/M=600/rD=0.1,rI=0.5,pS=0.5"))
  if(paste("Occu_Abu_sim",rep, ".RData", sep="") %!in% files){
    library(nimble)
    library(coda)
    
    
    data <- data.fn()
    true.vals[1,] <- data$mean.individual.detection.D
    true.vals[2,] <- data$mean.individual.detection.I
    true.vals[3,] <- data$mean.detection.s
    true.vals[4,] <- data$alpha1.D
    true.vals[5,] <- data$alpha2.D
    true.vals[6,] <- data$alpha1.I
    true.vals[7,] <- data$alpha2.I
    true.vals[8,] <- data$alpha1.s
    true.vals[9,] <- data$alpha2.s
    true.vals[10,] <- data$beta0.D
    true.vals[11,] <- data$beta0.I
    true.vals[12,] <- data$beta0.s
    true.vals[13,] <- data$beta1.D
    true.vals[14,] <- data$beta2.D
    true.vals[15,] <- data$beta3.D
    true.vals[16,] <- data$beta1.I
    true.vals[17,] <- data$beta2.I
    true.vals[18,] <- data$beta3.I
    true.vals[19,] <- data$beta1.s
    true.vals[20,] <- data$beta2.s
    true.vals[21,] <- data$beta3.s
    true.vals[22,] <- data$gamma0DI
    true.vals[23,] <- data$gamma0DS
    true.vals[24,] <- data$gamma0IS
    true.vals[25,] <- data$Z.total.s
    true.vals[26,] <- data$Ntotal.D
    true.vals[27,] <- data$Ntotal.I
    
    
    n.chains = 3
    chains = vector("list", n.chains)
    for(chain in 1:n.chains){
      # (M = 600, J = 4, mean.lambda.D = 0.5, beta1.D = 1, N.days = 7,
      #   beta2.D = -1, beta3.D = 1, mean.individual.detection.D = 0.5, alpha1.D = -1, alpha2.D = 0.15,
      #   mean.psi.s = 0.2, beta1.s = 1, beta2.s = -1, beta3.s = 1, mean.detection.s = 0.5,
      #   alpha1.s = -1, alpha2.s = 0.15, gamma0 = -1)
      # set gamma0 and gamma1 and save them
      
      #model run parameters
      n.iter <- 50000
      n.burn <- 20000
      
      #sampling parameters
      I =  600
      J = 4
      
      #pull out count data for dominant and subordinate species
      PA_highabundancesubordinate <- data$C.s
      counts_highabundancedominant <- data$C.D
      
      
      yD <- data$C.D
      yI <- data$C.I
      ys <- data$C.s
      
      
      #structure the data
      
      yD=array(yD,dim=c(I,J)) #new data 
      yI=array(yI, dim=c(I,J)) 
      yS=array(ys,dim=c(I,J)) 
      occ.new=array(data$occ,dim=c(I,J))
      day.new=array(data$day,dim=c(I,J))
      forest.new=array(data$forest,dim=c(I))
      road.new=array(data$road,dim=c(I))
      agr.new=array(data$agr,dim=c(I))
      
      #parameters for multinomial obs function
      D <- as.integer(7)
      O <- as.integer(1)
      
      
      #setting J for each year
      
      
      #Occupancy-abundance model with interaction parameters on p and multinomial observation model for cameras
      #create dataset
      str(fdata<- list(yD.new = yD, yI.new = yI, yS.new = yS, I = I, nsurveys = dim(yD)[2],
                       DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                       forestI = forest.new, roadI = road.new, agricultureI = agr.new,
                       forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = D, O = O))
      
      
      #Occu-abu interaction model run
      
      NimModel <- nimbleCode({
        # --- Priors ---
        #Model for dominant species (D) 
        #Priors for dominant species (D)
        alpha0D ~ dlogis(0,1)
        for (k in 1:2){
          alphaD[k] ~ dnorm(0, 0.1)
        }
        beta0D ~ dnorm(0, 0.1)
        for (k in 1:3){
          betaD[k] ~ dnorm(0, 0.1)
        }
        # Likelihood
        # State model for true abundance
        for (i in 1:I){
          nD[i] ~ dpois(lambdaD[i])
          log(lambdaD[i]) <- beta0D + betaD[1] * forestD[i] + betaD[2] * roadD[i] + betaD[3] * agricultureD[i] 
          
          #Observation model for replicated counts
          for (j in 1:nsurveys){
            yD.new[i,j] ~ dbinom(pD[i,j], N.days)
            logit(rD[i,j]) <- alpha0D + alphaD[1] * DATE[i,j] + alphaD[2] * occ[i,j] 
            pD[i,j] <- 1 - ((1-rD[i,j]) ^ nD[i])
          }
        }
        
        #Model for intermediate species (I)
        #Priors for intermediate species (I)
        alpha0I ~ dlogis(0,1)
        for (k in 1:2){
          alphaI[k] ~ dnorm(0, 0.1)
        }
        beta0I ~ dnorm(0, 0.1)
        for (k in 1:3){
          betaI[k] ~ dnorm(0, 0.1)
        }
        gamma0DI ~ dnorm(0, 0.2)
        # Likelihood for intermediate species (I)
        # State model for true abundance for intermediate species (I)
        for (i in 1:I){
          nI[i] ~ dpois(lambdaI[i])
          log(lambdaI[i]) <- beta0I + betaI[1] * forestI[i] + betaI[2] * roadI[i] + betaI[3] * agricultureI[i] + gamma0DI*nD[i]
          
          #Observation model for replicated counts
          for (j in 1:nsurveys){
            yI.new[i,j] ~ dbinom(pI[i,j], N.days)
            logit(rI[i,j]) <- alpha0I + alphaI[1] * DATE[i,j] + alphaI[2] * occ[i,j] 
            pI[i,j] <- 1 - ((1-rI[i,j]) ^ nI[i])
          }
        }
        #Model for subordinate species (S)
        #Priors for subordinate species (s)
        alpha0S ~ dlogis(0,1)
        
        for (k in 1:2){
          alphaS[k] ~ dnorm(0, 0.1)
        }
        beta0S ~ dlogis(0,1)
        for (k in 1:3){
          betaS[k] ~ dnorm(0, 0.1)
        }
        gamma0IS ~ dnorm(0, 0.2)
        gamma0DS ~ dnorm(0, 0.2)
        # gamma1 ~ dnorm(0, 0.2)
        # epsilon0 ~ dnorm(0, 0.2)
        # Likelihood for subordinate species (S)
        # State model for occupancy for subordinate species (S)
        for (i in 1:I){
          zS[i] ~ dbern(psiB[i])
          logit(psiB[i]) <- beta0S + betaS[1] * forestS[i] + betaS[2] * roadS[i] + betaS[3] * agricultureS[i] +gamma0IS *nI[i] +gamma0DS*nD[i]
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
        # Derived quantity: Total abundance across all surveyed sites
        totalnI <- sum(nI[1:I])
        
      }
      )
      
      
      # Parameters monitored
      parameters <-  c("alpha0D", "alphaD", "beta0D", "betaD","alpha0I", "alphaI", "beta0I", "betaI", "alpha0S", "alphaS", "beta0S", "betaS", "totalZs","totalnI", "totalnD", "gamma0DI", "gamma0IS", "gamma0DS")
      
      #data and constants
      Nimdata <- list(yD.new = as.matrix(yD), yS.new = as.matrix(yS), yI.new = as.matrix(yI))
      constants=list(I = I, nsurveys = dim(yD)[2],
                     DATE = day.new, occ = occ.new, forestD = forest.new, roadD = road.new, agricultureD = agr.new,
                     forestI = forest.new, roadI = road.new, agricultureI = agr.new,
                     forestS = forest.new, roadS = road.new, agricultureS = agr.new, N.days = D)
      
      #inits
      Niminits=list(beta0D=0,beta0I=0, beta0S = 0, nD=apply(yD,1,max),nI=apply(yI,1,max),zS=apply(yS,1,max),alpha0D=0, alpha0I=0, alpha0S=0, alphaD=c(0,0),alphaI=c(0,0),alphaS=c(0,0), betaD=c(0,0,0), betaI=c(0,0,0),betaS=c(0,0,0),gamma0DI=0, gamma0IS=0, gamma0DS=0)
      
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
      # plot(mcmc(mvSamples[-c(1:n.burn),]))
      chains[[chain]]=mvSamples
    }
    
    setwd(paste0(path,"/three_species/M=600/rD=0.1,rI=0.5,pS=0.5"))
    
    #save the chains
    # save(chains, file=paste("Occu_abu_chains_sim",rep,".RData", sep=""))
    n.iter <- 25000
    n.burn <- 20000
    #combine the chains and burn
    a=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
                mcmc(chains[[2]][n.burn:n.iter,]),
                mcmc(chains[[3]][n.burn:n.iter,]))
    
    gelman = gelman.diag(a)
    
    a=runjags::combine.mcmc(a)
    
    plot(a)
    # colnames(mvSamples)
    # [1] "alpha0D"   "alpha0I"   "alpha0S"   "alphaD[1]" "alphaD[2]" "alphaI[1]" "alphaI[2]" "alphaS[1]" "alphaS[2]" "beta0D"    "beta0I"    "beta0S"    "betaD[1]" 
    # [14] "betaD[2]"  "betaD[3]"  "betaI[1]"  "betaI[2]"  "betaI[3]"  "betaS[1]"  "betaS[2]"  "betaS[3]"  "gamma0DI"  "gamma0DS"  "gamma0IS"  "totalZs"   "totalnD"  
    # [27] "totalnI"
    
    #extract point estimates and credible intervals
    esti[1,1] <- plogis(mean(a[,"alpha0D"]))
    esti[c(2,3),1]<- plogis(quantile(a[,"alpha0D"], probs = c(2.5, 97.5)/100))
    
    esti[1,2] <- plogis(mean(a[,"alpha0I"]))
    esti[c(2,3),2]<- plogis(quantile(a[,"alpha0I"], probs = c(2.5, 97.5)/100))
    
    esti[1,3] <- plogis(mean(a[,"alpha0S"]))
    esti[c(2,3),3]<- plogis(quantile(a[,"alpha0S"], probs = c(2.5, 97.5)/100))
    
    esti[1,4] <- mean(a[,"alphaD[1]"])
    esti[c(2,3),4]<- quantile(a[,"alphaD[1]"], probs = c(2.5, 97.5)/100)
    
    esti[1,5] <- mean(a[,"alphaD[2]"])
    esti[c(2,3),5]<- quantile(a[,"alphaD[2]"], probs = c(2.5, 97.5)/100)
    
    esti[1,6] <- mean(a[,"alphaI[1]"])
    esti[c(2,3),6]<- quantile(a[,"alphaI[1]"], probs = c(2.5, 97.5)/100)
    
    esti[1,7] <- mean(a[,'alphaI[2]'])
    esti[c(2,3),7]<- quantile(a[,"alphaI[2]"], probs = c(2.5, 97.5)/100)
    
    esti[1,8] <- mean(a[,"alphaS[1]"])
    esti[c(2,3),8]<- quantile(a[,"alphaS[1]"], probs = c(2.5, 97.5)/100)
    
    esti[1,9] <- mean(a[,"alphaS[2]"])
    esti[c(2,3),9]<- quantile(a[,"alphaS[2]"], probs = c(2.5, 97.5)/100)
    
    esti[1,10] <- mean(a[,"beta0D"])
    esti[c(2,3),10]<- quantile(a[,"beta0D"], probs = c(2.5, 97.5)/100)
    
    esti[1,11] <- mean(a[,"beta0I"])
    esti[c(2,3),11]<- quantile(a[,"beta0I"], probs = c(2.5, 97.5)/100)
    
    esti[1,12] <- mean(a[,"beta0S"])
    esti[c(2,3),12]<- quantile(a[,"beta0S"], probs = c(2.5, 97.5)/100)
    
    esti[1,13] <- mean(a[,"betaD[1]"])
    esti[c(2,3),13]<- quantile(a[,"betaD[1]"], probs = c(2.5, 97.5)/100)
    
    esti[1,14] <- mean(a[,"betaD[2]"])
    esti[c(2,3),14]<- quantile(a[,"betaD[2]"], probs = c(2.5, 97.5)/100)
    
    esti[1,15] <- mean(a[,"betaD[3]"])
    esti[c(2,3),15]<- quantile(a[,"betaD[3]"], probs = c(2.5, 97.5)/100)
    
    esti[1,16] <- mean(a[,"betaI[1]"])
    esti[c(2,3),16]<- quantile(a[,"betaI[1]"], probs = c(2.5, 97.5)/100)
    
    esti[1,17] <- mean(a[,"betaI[2]"])
    esti[c(2,3),17]<- quantile(a[,"betaI[2]"], probs = c(2.5, 97.5)/100)
    
    esti[1,18] <- mean(a[,"betaI[3]"])
    esti[c(2,3),18]<- quantile(a[,"betaI[3]"], probs = c(2.5, 97.5)/100)
    
    esti[1,19] <- mean(a[,"betaS[1]"])
    esti[c(2,3),19]<- quantile(a[,"betaS[1]"], probs = c(2.5, 97.5)/100)
    
    esti[1,20] <- mean(a[,"betaS[2]"])
    esti[c(2,3),20]<- quantile(a[,"betaS[2]"], probs = c(2.5, 97.5)/100)
    
    esti[1,21] <- mean(a[,"betaS[3]"])
    esti[c(2,3),21]<- quantile(a[,"betaS[3]"], probs = c(2.5, 97.5)/100)
    
    esti[1,22] <- mean(a[,"gamma0DI"])
    esti[c(2,3),22]<- quantile(a[,"gamma0DI"], probs = c(2.5, 97.5)/100)
    
    esti[1,23] <- mean(a[,"gamma0DS"])
    esti[c(2,3),23]<- quantile(a[,"gamma0DS"], probs = c(2.5, 97.5)/100)
    
    esti[1,24] <- mean(a[,"gamma0IS"])
    esti[c(2,3),24]<- quantile(a[,"gamma0IS"], probs = c(2.5, 97.5)/100)
    
    esti[1,25] <- mean(a[,"totalZs"])
    esti[c(2,3),25]<- quantile(a[,"totalZs"], probs = c(2.5, 97.5)/100)
    
    esti[1,26] <- mean(a[,'totalnD'])
    esti[c(2,3),26]<- quantile(a[,"totalnD"], probs = c(2.5, 97.5)/100)
    
    esti[1,27] <- mean(a[,"totalnI"])
    esti[c(2,3),27]<- quantile(a[,"totalnI"], probs = c(2.5, 97.5)/100)
    
    
    #put all my shit together and save it.
    allmyshit <- list(esti, true.vals, gelman)
    save(allmyshit, file = paste("Occu_Abu_sim",rep, ".RData", sep=""))
    
    rm(chains)
    
    gc()
  }
}
stopCluster(cl.tmp)


