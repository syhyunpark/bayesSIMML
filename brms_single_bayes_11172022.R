# use brms package; This function fit a simple Bayesian linear model for the outcome
#no main effect of A
library("brms")
library("rstan")
bayes_single_model <- function(X, p,y,A){
  #fit_sng_data <- cbind(A,X,y)
  #colnames(fit_sng_data) <- c("A", paste0("X",1:p),"y")
  A.mean <- mean(A)
  Ac <- A - A.mean
  Ac.unique <- unique(sort(Ac))
  
  fit_sng_data <- data.frame(cbind(Ac,X,y))
  fit_sng_data$Ac <- factor(fit_sng_data$Ac)
  #setting contrasts in logistic model using brms: https://vasishth.github.io/bayescogsci/book/sec-contrast-covariate.html
  contrasts(fit_sng_data$Ac) <-  Ac.unique
  
  colnames(fit_sng_data) <- c("Ac", paste0("X",1:p),"y")
  
  if (p == 5){
    fit_single <-brm(formula = y ~  X1+X2+X3+X4+X5 + Ac*(X1+X2+X3+X4+X5),
                     family = bernoulli(link = "logit"),
                     data=fit_sng_data,
                     iter = 3000, warmup = 500, chains = 4,
                     prior = c(prior(normal(0,5),class=b),
                               prior(student_t(3,0,8), class=Intercept)))
  }
  if (p == 10){
    fit_single <-brm(formula = y ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10 + 
                       Ac*(X1+X2+X3+X4+X5+X6+X7+X8+X9+X10), 
                     family = bernoulli(link = "logit"),
                     data=fit_sng_data,
                     iter = 3000, warmup = 500, chains = 4,
                     prior = c(prior(normal(0,5),class=b),
                               prior(student_t(3,0,8), class=Intercept)))
  }
  
  #head(as.matrix(fit_single))
  # Get posterior draws of all parameters
  sparams <- get_sampler_params(fit_single$fit)
  div_single <- sum(sapply(sparams, function(x) sum(x[, 'divergent__'])))
  
  beta_0<- as.matrix(
    fit_single,
    variable = "b_Ac1",
  ) # estimation of treatment main effect
  
  beta_1 <- as.matrix(
    fit_single,
    variable = paste0("b_X",1:p,":Ac1"),
  ) #estimation of coefficients for treatment by covariate interaction
  m <- as.matrix(
    fit_single,
    variable = paste0("b_X",1:p),
  ) #estimation of covariates' main effects 
  tau <- as.matrix(
    fit_single,
    variable = "b_Intercept",
  ) # intercept
  return(data.table(beta_0,beta_1,
                    m, tau, div_single)) 
}

pred_sng_bayes <- function(xnew=NULL, xmnew=NULL, newA=NULL, A.mean = NULL, coefs = sgbayes.coefs,thresh = 0){
  ##Given patients characteristics, derive optimal treatment for each patients using simple Bayesian linear model
  
  n <- nrow(xnew)
  p <- ncol(xnew)
  
  beta_0 <- as.matrix(coefs[,c("b_Ac1")])
  beta_1 <- as.matrix(coefs%>%select(paste0("b_X",1:p,":Ac1"))) 
  m <- as.matrix(coefs%>%select(paste0("b_X",1:p))) #estimation of covariates' main effects 
  tau <- as.matrix(coefs[,c("b_Intercept")]) # int
  
  
  nsample <-  nrow(beta_0)#number of MCMC samples
  
  eta.distr = contrast.distr = tbi.tmp  <- matrix(rep(0, n*nsample), n, nsample)
  
  newA.mean <- newA-A.mean
  
  for(i in 1:n){
    eta.distr[i,] <- tau + m%*%xnew[i,] + beta_0*newA.mean[i] + beta_1%*%xnew[i,]*newA.mean[i]#eta.distr (n.test*n.samples): a canonical parameter matrix with a row for each observation(test data point) and a column for each posterior sample
    contrast.distr[i,] <- beta_0 + beta_1%*%xnew[i,]# samples of distribution of tbi (canonical parameter of A=2 - canonical parameter of A=1)
    tbi.tmp[i,] <- (contrast.distr[i,] < thresh)#each row: ? tbi < 0 for the ith patient
  }
  
  tbi <- apply(tbi.tmp, 1, mean)#each row: Pr(tbi < 0) for the ith patient
  
  return(list(eta.distr=eta.distr, contrast.distr =contrast.distr, tbi=tbi))
}

