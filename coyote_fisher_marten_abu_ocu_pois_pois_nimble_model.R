#occu-abu model for fisher and marten

#Occu-abu interaction model run
library(nimble)

NimModel <- nimbleCode({
  # --- Priors ---
  #Model for dominant species (D)
  #Priors for dominant species (D)
  alpha0D ~ dlogis(0,1)
  beta0D ~ dnorm(0, 0.1)
  for (k in 1:3) {
    alphaD[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaD[k] ~ dnorm(0, 0.1)
  }
  
  tau.pD <- pow(sd.pD, -2)
  sd.pD ~ dunif(0, 2)
  tau.lamD <- pow(sd.lamD, -2)
  sd.lamD ~ dunif(0, 2)
  # Likelihood
  # State model for true abundance of dominant species (D)
  #year random intercepts
  
  for (t in 1:nyears) {
    beta0D.it[t] ~ dnorm(beta0D, tau.lamD)
  }
  
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      nD[j,t] ~ dpois(lambda[j,t])
      log(lambda[j,t]) <-  beta0D.it[t] + betaD[1] * forestD[j,t] + betaD[2] * deerD[j,t] + betaD[3] * edgeD[j,t] 
    }
  }
  for (t in 1:nyears) {
    #site random intercepts nested within year
    alpha0D.it[t] ~ dnorm(alpha0D, tau.pD)
  }
  #Observation model for replicated counts
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      for (k in 1:nsurveys){
        yD.new[j,k,t] ~ dpois(pD[j,k,t] * nD[j,t])
        log(pD[j,k,t]) <- alpha0D.it[t] + alphaD[1] * DATE[j,k,t] + alphaD[2] * pow(DATE[j,k,t], 2) + alphaD[3] * occ[j,k,t]
      }
    }
  }
  
  # --- Priors ---
  #Model for dominant species (I)
  #Priors for intermediate species (I)
  alpha0I<- logit(mean.pI)
  mean.pI ~ dunif(0, 1)
  beta0I ~ dnorm(0, 0.1)
  for (k in 1:3) {
    alphaI[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaI[k] ~ dnorm(0, 0.1)
  }
  
  tau.pI <- pow(sd.pI, -2)
  sd.pI ~ dunif(0, 2)
  tau.lamI <- pow(sd.lamI, -2)
  sd.lamI ~ dunif(0, 2)
  gamma0DI ~ dnorm(0, 0.2)
  # Likelihood
  # State model for true abundance of dominant species (I)
  #year random intercepts
  
  for (t in 1:nyears) {
    beta0I.it[t] ~ dnorm(beta0I, tau.lamI)
  }
  
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      nI[j,t] ~ dpois(lambdaI[j,t])
      log(lambdaI[j,t]) <-  beta0I.it[t] + betaI[1] * conifmixedI[j,t] + betaI[2] * snowdepthI[j,t] + betaI[3] * deciduousI[j,t] + gamma0DI * nD[j,t]
    }
  }
  for (t in 1:nyears) {
    #site random intercepts nested within year
    alpha0I.it[t] ~ dnorm(alpha0I, tau.pI)
  }
  #Observation model for replicated counts
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      for (k in 1:nsurveys){
        yI.new[j,k,t] ~ dpois(pI[j,k,t] * nI[j,t])
        log(pI[j,k,t]) <- alpha0I.it[t] + alphaI[1] * DATE[j,k,t] + alphaI[2] * pow(DATE[j,k,t], 2) + alphaI[3] * occ[j,k,t]
      }
    }
  }
  #Model for subordinate species (S)
  #Priors for subordinate species (s)
  alpha0S <- logit(mean.pS)
  mean.pS ~ dunif(0, 1)
  beta0S ~ dnorm(0, 0.1)
  tau.pS <- pow(sd.pS, -2)
  sd.pS ~ dunif(0, 3)
  tau.psiS <- pow(sd.psiS, -2)
  sd.psiS ~ dunif(0, 3)
  gamma0DS ~ dnorm(0, 0.2)
  gamma0IS ~ dnorm(0, 0.2)
  for (k in 1:3) {
    alphaS[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaS[k] ~ dnorm(0, 0.1)
  }
  # # Likelihood for subordinate species (S)
  # #year level random effect on psi
  for (t in 1:nyears) {
    beta0S.it[t] ~ dnorm(beta0S, tau.psiS)
  }
  # State model for occupancy for subordinate species (S)
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      zS[j, t] ~ dbern(psi[j, t])
      logit(psi[j, t]) <-  beta0S.it[t] + betaS[1] * conifmixedS[j, t] + betaS[2] * deciduousS[j, t] + betaS[3] * snowdepthS[j, t] + gamma0DS * nD[j, t] + gamma0IS * nI[j, t] 
    }
  }
  #year level random effect on p
  for (t in 1:nyears) {
    alpha0S.it[t] ~ dnorm(alpha0S, tau.pS)
    
  }
  #Observation model for replicated cout
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      for (k in 1:nsurveys) {
        yS.new[j, k, t] ~ dbern((pS[j, k, t] * K2D[j,k,t]) * zS[j, t])
        logit(pS[j, k, t]) <-alpha0S.it[t] + alphaS[1] * DATE[j, k, t] +  alphaS[2] * pow(DATE[j,k,t], 2) + alphaS[3] *  occ[j, k, t]
        
      }
    }
  }
  
  # Derived quantity: Occupancy across all surveyed sites
  # Derived variable
  for (t in 1:nyears) {
    totalZs[t] <- sum(zS[1:J[t],t])
  }
  
  # Derived quantity: Total abundance across all surveyed sites
  for (t in 1:nyears) {
    totalnD[t] <- sum(nD[1:J[t],t])
  }
  
  # Derived quantity: Total abundance across all surveyed sites
  for (t in 1:nyears) {
    totalnI[t] <- sum(nI[1:J[t],t])
  }
  
  # Derived quantity: Total abundance across all surveyed sites
  for (t in 1:nyears) {
    meannD[t] <- mean(nD[1:J[t],t])
  }
  
  # Derived quantity: Total abundance across all surveyed sites
  for (t in 1:nyears) {
    maxnD[t] <- max(nD[1:J[t],t])
  }
  #components for fit diagnostics #DOMINANT
  for (t in 1:nyears) {
    for (j in 1:J[t]){
      for (k in 1:nsurveys) {
        YsimD[j,k,t] ~  dpois(pD[j,k,t] * nD[j,t])
        YexpD[j,k,t] <- pD[j,k,t] * nD[j,t]

        #components for test 1
        err1obsD[j,k,t] <- (sqrt(yD.new[j,k,t]) - sqrt(YexpD[j,k,t]))^2
        err1simD[j,k,t] <- (sqrt(YsimD[j,k,t]) - sqrt(YexpD[j,k,t]))^2
      }
    }
  }
  for (t in 1:nyears){
    T1obstempD[t] <- sum(err1obsD[1:J[t], 1:nsurveys, t])
    T1simtempD[t] <- sum(err1simD[1:J[t], 1:nsurveys, t])
  }

  #fit diagnostic totals
  T1obsD <- sum(T1obstempD[1:nyears])
  T1simD <- sum(T1simtempD[1:nyears])

  #components for fit diagnostics #INTERMEDIATE
  for (t in 1:nyears) {
    for (j in 1:J[t]){
      for (k in 1:nsurveys) {
        YsimI[j,k,t] ~  dpois(pI[j,k,t] * nI[j,t])
        YexpI[j,k,t] <- pI[j,k,t] * nI[j,t]

        #components for test 1
        err1obsI[j,k,t] <- (sqrt(yI.new[j,k,t]) - sqrt(YexpI[j,k,t]))^2
        err1simI[j,k,t] <- (sqrt(YsimI[j,k,t]) - sqrt(YexpI[j,k,t]))^2
      }
    }
  }
  for (t in 1:nyears){
    T1obstempI[t] <- sum(err1obsI[1:J[t], 1:nsurveys, t])
    T1simtempI[t] <- sum(err1simI[1:J[t], 1:nsurveys, t])
  }

  #fit diagnostic totals
  T1obsI <- sum(T1obstempI[1:nyears])
  T1simI <- sum(T1simtempI[1:nyears])

  #components for fit diagnostics #SUBORDINATE
  for (t in 1:nyears) {
    for (j in 1:J[t]){
      for (k in 1:nsurveys) {
        YsimS[j,k,t] ~ dbern((pS[j,k,t] * K2D[j,k,t]) * zS[j, t])
        YexpS[j,k,t] <- (pS[j,k,t] * K2D[j,k,t]) * zS[j, t]

        #components for test 1
        err1obsS[j,k,t] <- (sqrt(yS.new[j,k,t]) - sqrt(YexpS[j,k,t]))^2
        err1simS[j,k,t] <- (sqrt(YsimS[j,k,t]) - sqrt(YexpS[j,k,t]))^2
      }
    }
  }
  for (t in 1:nyears){
    T1obstempS[t] <- sum(err1obsS[1:J[t], 1:nsurveys, t])
    T1simtempS[t] <- sum(err1simS[1:J[t], 1:nsurveys, t])
  }
  #fit diagnostic totals
  T1obsS <- sum(T1obstempS[1:nyears])
  T1simS <- sum(T1simtempS[1:nyears])
} )