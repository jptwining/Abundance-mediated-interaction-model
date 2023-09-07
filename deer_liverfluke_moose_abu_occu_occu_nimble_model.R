NimModel <- nimbleCode({
  # --- Priors ---
  #Model for dominant species (D)
  #Priors for dominant species (D)
  beta0D ~ dnorm(0, 0.1) #use fixed effects
  alpha0D ~ dnorm(0, 0.1) #use fixed effects
  
  for (k in 1:2) {
    alphaD[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaD[k] ~ dnorm(0, 0.1)
  }
  
  # Likelihood
  # State model for true abundance of dominant species (D)
  for (j in 1:J) {
    nD[j] ~ dpois(lambda[j])
    log(lambda[j]) <-  beta0D + betaD[1] * forestD[j] + betaD[2] * roadD[j] + betaD[3] * agricultureD[j] 
  }
  
  #Observation model for replicated counts
  for (j in 1:J) {
    for (k in 1:nsurveys){
      logit(pD[j,k]) <- alpha0D + alphaD[1] * DATE[j,k] + alphaD[2] * occ[j,k]
      pD.cum[j,k] <- 1 - (1 - pD[j,k])^nD[j] #prob detecting at least 1 ind per day
      yD.new[j,k] ~ dbinom(pD.cum[j,k],size=D)
    }
  }
  
  
  # --- Priors ---
  #Model for dominant species (I)
  #Priors for intermediate species (I)
  beta0I ~ dlogis(0,1) #use fixed effects
  alpha0I ~ dlogis(0,1) #use fixed effects
  for (k in 1:2) {
    alphaI[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:3) {
    betaI[k] ~ dnorm(0, 0.1)
  }
  gamma0DI ~ dnorm(0, 0.2)
  
  # Likelihood
  # State model for true abundance of dominant species (I)
  for (j in 1:J) {
    zI[j] ~ dbern(psiI[j])
    logit(psiI[j]) <- beta0I + betaI[1] * forestI[j] + betaI[2] * roadI[j] + betaI[3] * agricultureI[j] + gamma0DI * nD[j]
  }
  
  #Observation model for replicated counts
  for (j in 1:J) {
    for (k in 1:nsurveys){
      yI.new[j, k] ~ dbern(pI[j, k] * zI[j])
      logit(pI[j, k]) <-alpha0I + alphaI[1] * DATE[j, k] + alphaI[2] *  occ[j, k]
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
  for (j in 1:J) {
    zS[j] ~ dbern(psiS[j])
    logit(psiS[j]) <-  beta0S + betaS[1] * forestS[j] + betaS[2] * roadS[j] + betaS[3] * agricultureS[j] + gamma0IS * zI[j] 
  }
  
  #Observation model for replicated count
  for (j in 1:J) {
    for (k in 1:nsurveys) {
      yS.new[j, k] ~ dbern(pS[j, k] * zS[j])
      logit(pS[j, k]) <-alpha0S + alphaS[1] * DATE[j, k] + alphaS[2] *  occ[j, k]
      
    }
  }
} )
