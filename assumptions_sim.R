rm(list=ls())
library(np); library(sn)

estimator <- function(DATA,est='probit',theta=1) {
  # Estimating P(D = 1|Z,X)
  if (est == 'probit') {
    fit <- glm(D~1+Z+X,data=DATA,family=binomial('probit'))
    G <- predict(fit,type='link'); G <- (G-mean(G))/sqrt(var(G)) # Standardise variance 
  }
  if (est == 'np') {
    invisible(capture.output(bw <- npindexbw(formula=D~1+Z+X,data=DATA,method='kleinspady')))
    invisible(capture.output(fit <- npindex(formula=D~1+Z+X,bws=bw,data=DATA,method='kleinspady')))
    G <- fit$beta[1]*DATA$Z + fit$beta[2]*DATA$X; G <- (G-mean(G))/sqrt(var(G)) # Standardise variance 
  }
  if (est == 'logit') {
    fit <- glm(D~1+Z+X,data=DATA,family=binomial)
    G <- predict(fit,type='link'); G <- (G-mean(G))/sqrt(var(G)) # Standardise variance 
  }

  # Compute estimate
  beta <- (cov(Z,Y)/cov(Z,G))/sqrt(theta)
  return(beta)
}

# Simulation parameters
alpha0 <- 0; alphaz <- 1; alphax <- -1
beta0 <- 0; betal <- 1; betax <- 0.2; betav <- 0.2

# Exponential distribution parameters
lambda <- 0.5; varX <- 1/lambda^2

# Skew normal parameters
loc <- 0; scale <- 1

# Theta^2: share of variation in total liability which is genetic
theta <- 0.1

# Share of variance in G from instrument
s <- 0.4

# Sample size of N, Number of simulation samples of nruns
N <- 1000; nruns <- 1000

a <- c(0,1,2,3,4,5)
b <- 0 # c(0,0.1,0.25,0.5,1)
grid <- as.matrix(expand.grid(a,b))
ng <- NROW(grid)

simout <- as.data.frame(matrix(nrow=ng,ncol=8))
colnames(simout) <- c('a','b','p','p_se', 'n', 'n_se', 'l', 'l_se')
simout[,c('a','b')] <- grid

for (gval in 1:ng) {
  sim_track <- matrix(nrow=nruns,ncol=3)
  colnames(sim_track) <- c('p','n','l')
  
  # Simulation parameters for this run
  a <- simout[gval,'a']
  b <- simout[gval,'b']
  
  varV <- scale^2*(1-2*(a^2/(1+a^2))/pi)
  
  for (k in 1:nruns) {
    set.seed(k)
    # Components of latent liability
    Z <- rnorm(N,0,sqrt(s*theta)); X <- (rexp(N,lambda)-1/lambda)*sqrt(((1-s)*theta)/varX)
    V <- (rsn(N,loc,scale,a) - loc - sqrt(2/pi)*scale*a/sqrt(1+a^2))*sqrt((1-theta)/varV)
    
    # True genetic liability
    Gtrue <- (alpha0 + alphaz*Z + alphax*X) #; Gtrue <- Gtrue/sqrt(var(Gtrue)/gamma)
    
    # Latent liability (Genetic minus non-genetic)0.0
    L <- Gtrue - V
    
    # Observed dichotomisation (threshold crossing)
    D <- ifelse(L >= b*X,1,0)
    
    # Constructing the outcome (linear function of CDF of latent liability)
    Y <- beta0 + betal*L + betax*X + betav*V + rnorm(N,0,1)
    
    DATA <- data.frame(Z,X,D,Y)
    sim_track[k,'p'] <- estimator(DATA,est='probit',theta=theta)
    sim_track[k,'n'] <- estimator(DATA,est='np',theta=theta)
    sim_track[k,'l'] <- estimator(DATA,est='logit',theta=theta)
    
    if (k%%1 == 0) cat(paste("Parameter ", gval, ", Iteration ",k, "\n   Probit: ", round(mean(sim_track[,'p'], na.rm=T)/betal,3),
                         "\n   Non-parametric: ", round(mean(sim_track[,'n'], na.rm=T)/betal,3),
                         "\n   Logit: ", round(mean(sim_track[,'l'], na.rm=T)/betal,3), "\n", sep=""))

  }
  simout[gval,c('p','p_se','n','n_se','l','l_se')] <- 
    c(mean(sim_track[,'p'])/betal, sqrt(var(sim_track[,'p'])/nruns),
      mean(sim_track[,'n'])/betal, sqrt(var(sim_track[,'p'])/nruns),
      mean(sim_track[,'l'])/betal, sqrt(var(sim_track[,'p'])/nruns))
  
  cat(paste(
    "Probit 95% CI: [", round(simout$p[gval] - 1.96*simout$p_se[gval],3), ", ", round(simout$p[gval] + 1.96*simout$p_se[gval],3), "] \n",
    "Non-parametric 95% CI: [", round(simout$n[gval] - 1.96*simout$n_se[gval],3), ", ", round(simout$n[gval] + 1.96*simout$n_se[gval],3), "] \n",
    "Logit 95% CI: [", round(simout$l[gval] - 1.96*simout$l_se[gval],3), ", ", round(simout$l[gval] + 1.96*simout$l_se[gval],3), "] \n", sep="")
  )
}
