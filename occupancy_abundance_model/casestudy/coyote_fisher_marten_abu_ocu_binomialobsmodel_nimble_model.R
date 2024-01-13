#occu-abu model for fisher and marten

#Occu-abu interaction model run
library(nimble)

NimModel <- nimbleCode({
  # --- Priors ---
  #Model for dominant species (D)
  #Priors for dominant species (D)
  for (t in 1:nyears) {
    beta0D[t] ~ dnorm(0, 0.1) #use fixed effects
  }
  for (t in 1:nyears) {
    alpha0D[t] ~ dnorm(0, 0.1) #use fixed effects
  }
    for (k in 1:3) {
    alphaD[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaD[k] ~ dnorm(0, 0.1)
  }

  # Likelihood
  # State model for true abundance of dominant species (D)
  
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      nD[j,t] ~ dpois(lambda[j,t])
      log(lambda[j,t]) <-  beta0D[t] + betaD[1] * forestD[j,t] + betaD[2] * deerD[j,t] + betaD[3] * edgeD[j,t] 
    }
  }

  #Observation model for replicated counts
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      for (k in 1:nsurveys){
        logit(pD[j,k,t]) <- alpha0D[t] + alphaD[1] * DATE[j,k,t] + alphaD[2] * pow(DATE[j,k,t], 2) + alphaD[3] * occ[j,k,t]
        pD.cum[j,k,t] <- 1 - (1 - pD[j,k,t])^nD[j,t] #prob detecting at least 1 ind per day
        yD.new[j,k,t] ~ dbinom(pD.cum[j,k,t],size=D)
      }
    }
  }
  
  # --- Priors ---
  #Model for dominant species (I)
  #Priors for intermediate species (I)
  
  for (t in 1:nyears) {
    beta0I[t] ~ dnorm(0, 0.1) #use fixed effects
  }
  for (t in 1:nyears) {
    alpha0I[t] ~ dnorm(0, 0.1) #use fixed effects
  }
  for (k in 1:3) {
    alphaI[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaI[k] ~ dnorm(0, 0.1)
  }
  gamma0DI ~ dnorm(0, 0.2)
  
  # Likelihood
  # State model for true abundance of dominant species (I)
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      nI[j,t] ~ dpois(lambdaI[j,t])
      log(lambdaI[j,t]) <- beta0I[t] + betaI[1] * conifmixedI[j,t] + betaI[2] * snowdepthI[j,t] + betaI[3] * deciduousI[j,t] + gamma0DI * nD[j,t]
    }
  }

  #Observation model for replicated counts
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      for (k in 1:nsurveys){
        logit(pI[j,k,t]) <- alpha0I[t] + alphaI[1] * DATE[j,k,t] + alphaI[2] * pow(DATE[j,k,t], 2) + alphaI[3] * occ[j,k,t]
        pI.cum[j,k,t] <- 1 - (1 - pI[j,k,t])^nI[j,t] #prob detecting at least 1 ind per day
        yI.new[j,k,t] ~ dbinom(pI.cum[j,k,t],size=D)
      }
    }
  }
  #Model for subordinate species (S)
  #Priors for subordinate species (s)
  for (t in 1:nyears) {
    beta0S[t] ~ dlogis(0,1) #use fixed effects
  }
  for (t in 1:nyears) {
    alpha0S[t] ~ dlogis(0,1) #use fixed effects
  }
  for (k in 1:3) {
    alphaS[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaS[k] ~ dnorm(0, 0.1)
  }
  gamma0DS ~ dnorm(0, 0.2)
  gamma0IS ~ dnorm(0, 0.2)
  
    # # Likelihood for subordinate species (S)
  # State model for occupancy for subordinate species (S)
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      zS[j, t] ~ dbern(psi[j, t])
      logit(psi[j, t]) <-  beta0S[t] + betaS[1] * conifmixedS[j, t] + betaS[2] * deciduousS[j, t] + betaS[3] * snowdepthS[j, t] + gamma0DS * nD[j, t] + gamma0IS * nI[j, t] 
    }
  }

  #Observation model for replicated count
  for (t in 1:nyears) {
    for (j in 1:J[t]) {
      for (k in 1:nsurveys) {
        yS.new[j, k, t] ~ dbern((pS[j, k, t] * K2D[j,k,t]) * zS[j, t])
        logit(pS[j, k, t]) <-alpha0S[t] + alphaS[1] * DATE[j, k, t] +  alphaS[2] * pow(DATE[j,k,t], 2) + alphaS[3] *  occ[j, k, t]
        
      }
    }
  }
  
  # # Derived quantity: Occupancy across all surveyed sites
  # # Derived variable
  # for (t in 1:nyears) {
  #   totalZs[t] <- sum(zS[1:J[t],t])
  # }
  # 
  # # Derived quantity: Total abundance across all surveyed sites
  # for (t in 1:nyears) {
  #   totalnD[t] <- sum(nD[1:J[t],t])
  # }
  # 
  # # Derived quantity: Total abundance across all surveyed sites
  # for (t in 1:nyears) {
  #   totalnI[t] <- sum(nI[1:J[t],t])
  # }
  # 
  # # Derived quantity: Total abundance across all surveyed sites
  # for (t in 1:nyears) {
  #   meannD[t] <- mean(nD[1:J[t],t])
  # }
  # 
  # # Derived quantity: Total abundance across all surveyed sites
  # for (t in 1:nyears) {
  #   maxnD[t] <- max(nD[1:J[t],t])
  # }
  # #components for fit diagnostics #DOMINANT
  # for (t in 1:nyears) {
  #   for (j in 1:J[t]){
  #     for (k in 1:nsurveys) {
  #       YsimD[j,k,t] ~  dpois(pD[j,k,t] * nD[j,t])
  #       YexpD[j,k,t] <- pD[j,k,t] * nD[j,t]
  # 
  #       #components for test 1
  #       err1obsD[j,k,t] <- (sqrt(yD.new[j,k,t]) - sqrt(YexpD[j,k,t]))^2
  #       err1simD[j,k,t] <- (sqrt(YsimD[j,k,t]) - sqrt(YexpD[j,k,t]))^2
  #     }
  #   }
  # }
  # for (t in 1:nyears){
  #   T1obstempD[t] <- sum(err1obsD[1:J[t], 1:nsurveys, t])
  #   T1simtempD[t] <- sum(err1simD[1:J[t], 1:nsurveys, t])
  # }
  # 
  # #fit diagnostic totals
  # T1obsD <- sum(T1obstempD[1:nyears])
  # T1simD <- sum(T1simtempD[1:nyears])
  # 
  # #components for fit diagnostics #INTERMEDIATE
  # for (t in 1:nyears) {
  #   for (j in 1:J[t]){
  #     for (k in 1:nsurveys) {
  #       YsimI[j,k,t] ~  dpois(pI[j,k,t] * nI[j,t])
  #       YexpI[j,k,t] <- pI[j,k,t] * nI[j,t]
  # 
  #       #components for test 1
  #       err1obsI[j,k,t] <- (sqrt(yI.new[j,k,t]) - sqrt(YexpI[j,k,t]))^2
  #       err1simI[j,k,t] <- (sqrt(YsimI[j,k,t]) - sqrt(YexpI[j,k,t]))^2
  #     }
  #   }
  # }
  # for (t in 1:nyears){
  #   T1obstempI[t] <- sum(err1obsI[1:J[t], 1:nsurveys, t])
  #   T1simtempI[t] <- sum(err1simI[1:J[t], 1:nsurveys, t])
  # }
  # 
  # #fit diagnostic totals
  # T1obsI <- sum(T1obstempI[1:nyears])
  # T1simI <- sum(T1simtempI[1:nyears])
  # 
  # #components for fit diagnostics #SUBORDINATE
  # for (t in 1:nyears) {
  #   for (j in 1:J[t]){
  #     for (k in 1:nsurveys) {
  #       YsimS[j,k,t] ~ dbern((pS[j,k,t] * K2D[j,k,t]) * zS[j, t])
  #       YexpS[j,k,t] <- (pS[j,k,t] * K2D[j,k,t]) * zS[j, t]
  # 
  #       #components for test 1
  #       err1obsS[j,k,t] <- (sqrt(yS.new[j,k,t]) - sqrt(YexpS[j,k,t]))^2
  #       err1simS[j,k,t] <- (sqrt(YsimS[j,k,t]) - sqrt(YexpS[j,k,t]))^2
  #     }
  #   }
  # }
  # for (t in 1:nyears){
  #   T1obstempS[t] <- sum(err1obsS[1:J[t], 1:nsurveys, t])
  #   T1simtempS[t] <- sum(err1simS[1:J[t], 1:nsurveys, t])
  # }
  # #fit diagnostic totals
  # T1obsS <- sum(T1obstempS[1:nyears])
  # T1simS <- sum(T1simtempS[1:nyears])
} )
