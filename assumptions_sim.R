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
beta0 <- 0; betal <- 1; betax <- 1; betav <- 1

# Exponential distribution parameters
lambda <- 0.5; varX <- 1/lambda^2

# Skew normal parameters
loc <- 0; scale <- 1

# Theta^2: share of variation in total liability which is genetic
theta <- 0.05

# Share of variance in G from instrument
s <- 0.2

# Sample size of N, Number of simulation samples of nruns
N <- 5000; nruns <- 1000

a <- 0
b <- c(0,0.1,0.25,0.5,1)
grid <- as.matrix(expand.grid(a,b))
ng <- NROW(grid)

sim_output <- matrix(nrow=ng,ncol=5)
colnames(sim_output) <- c('a','b','probit','np','logit')
sim_output[,c('a','b')] <- grid

for (gridval in 1:ng) {
  sim_track <- matrix(nrow=nruns,ncol=3)
  colnames(sim_track) <- c('probit','np','logit')
  
  # Simulation parameters for this run
  a <- sim_output[gridval,'a']
  b <- sim_output[gridval,'b']
  
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
    sim_track[k,'probit'] <- estimator(DATA,est='probit',theta=theta)
    sim_track[k,'np'] <- estimator(DATA,est='np',theta=theta)
    sim_track[k,'logit'] <- estimator(DATA,est='logit',theta=theta)
    
    if (k%%1 == 0) print(paste("Iteration: ",k, " --- Probit: ", round(median(sim_track[,'probit'], na.rm=T),3),
                                " --- Non-parametric: ", round(median(sim_track[,'np'], na.rm=T),3),
                                " --- Logit: ", round(median(sim_track[,'logit'], na.rm=T),3), sep=""))
  }
  sim_output[gridval,c('probit','np','logit')] <- c(100*abs((median(sim_track[,'probit'])-betal)/betal),
                         100*abs((median(sim_track[,'np'])-betal)/betal),
                         100*abs((median(sim_track[,'logit'])-betal)/betal))
  
  print(paste("Probit 95% CI: [", round(quantile(sim_track[,'probit'],0.025),3), ", ", 
              round(quantile(sim_track[,'probit'],0.975),3), "] -- Logit 95% CI: [",
              round(quantile(sim_track[,'logit'],0.025),3), ", ", 
              round(quantile(sim_track[,'logit'],0.975),3), "]", sep=""))
}
