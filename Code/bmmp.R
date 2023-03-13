##### Bayesian Multivariate Matched Proportions #####
# R functions to run Bayesian Multivariate Matched Proportions		
#	from Meyer, Cheng, and Knutson (2023) 									
#																	
# Created:  07/25/2018												
# Modified: 11/09/2022												
#																	
# By: Mark J Meyer													
#																	



#### Require libraries ####
suppressPackageStartupMessages(require(MASS))
suppressPackageStartupMessages(require(MCMCpack))
suppressPackageStartupMessages(require(mvnfast))
suppressPackageStartupMessages(require(mvtnorm)) # drop?
suppressPackageStartupMessages(require(numbers))
suppressPackageStartupMessages(require(truncnorm))
suppressPackageStartupMessages(require(geepack))
suppressPackageStartupMessages(require(statmod))
# suppressPackageStartupMessages(require(MCMCglmm))



#### General Use Functions ####
getGelman	<- function(mcmcOut, chains = 4){
  if(chains != 4 & chains != 2){
    stop('Split chains in half or quarters only')
  }
  K		<- ncol(mcmcOut)
  B		<- nrow(mcmcOut)
  GRR		<- matrix(0, ncol = 2, nrow = K)
  cuts		<- split(1:B, rep(1:chains, each = B/chains))
  
  for(k in 1:K){
    if(chains == 4){
      mcmc1	<- mcmc(mcmcOut[cuts[[1]],k])
      mcmc2	<- mcmc(mcmcOut[cuts[[2]],k])
      mcmc3	<- mcmc(mcmcOut[cuts[[3]],k])
      mcmc4	<- mcmc(mcmcOut[cuts[[4]],k])
      gdb		<- gelman.diag(mcmc.list(mcmc1, mcmc2, mcmc3, mcmc4))
      
      GRR[k,]	<- gdb$psrf
    } else {
      mcmc1	<- mcmc(mcmcOut[cuts[[1]],k])
      mcmc2	<- mcmc(mcmcOut[cuts[[2]],k])
      gdb		<- gelman.diag(mcmc.list(mcmc1, mcmc2))
      
      GRR[k,]	<- gdb$psrf		
    }
  }
  colnames(GRR)	<- c('Point est.', 'Upper C.I.')
  rownames(GRR)	<- 1:K
  
  return(GRR)
  
}

getCounts	<- function(X1, X2, K){
  n21		<- vector('numeric', length = K)
  n12		<- vector('numeric', length = K)
  n11		<- vector('numeric', length = K)
  n22		<- vector('numeric', length = K)
  
  for(j in 1:K){
    rawtab	<- table(X2[,j], X1[,j])
    temp	<- matrix(0, nrow = 2, ncol = 2)
    temp[1:dim(rawtab)[1], 1:dim(rawtab)[2]]	<- rawtab
    
    patab	<- matrix(c(temp[2,2], temp[1,2], temp[2,1], temp[1,1]), nrow = 2)
    
    n21[j]	<- patab[2,1]
    n12[j]	<- patab[1,2]
    
    n11[j]	<- patab[1,1]
    n22[j]	<- patab[2,2]
  }
  
  l	<- list(n21 = n21, n12 = n12, n11 = n11, n22 = n22)
  return(l)
}



#### Bayesian Probit Model ####
bpm  <- function(X1, X2, B = 10000, burnin = B,
				 verbose = TRUE, up = 2000, dots = 200){
  
  ### data ###
  n		    <- nrow(X1)
  K		    <- ncol(X1)
  ymat		<- as.matrix(cbind(X1, X2))
  bDat	  <- data.frame(Y	= c(t(ymat)), s = rep(1:(2*K), n), id = rep(1:n, each = 2*K))
  form    <- model.frame(Y ~ factor(s) - 1, data = bDat)
  Y       <- model.response(form)
  X       <- model.matrix(form, data = bDat)
  n1		  <- sum(Y)
  n0		  <- length(Y) - n1
  
  model   <- list(Y = Y, X = X, bdat = bDat)
  mcmcs   <- list(B = B, burnin = burnin)
  
  # Vb      <- solve(t(X)%*%X)
  XtX     <- t(X)%*%X
  Pb          <- XtX
  Vb          <- solve(Pb)
  
  # nL      <- ncol(L)
  nX      <- ncol(X)
  n       <- nrow(X)
  
  ### build matrices ###
  betam   <- matrix(0, nrow = B + burnin, ncol = nX)
  Z       <- rep(0, n)
  
  ### sampler ###
  for(b in 2:(B + burnin)){
    
    ### sample latent Zs ###
    XB			<- X%*%betam[b-1,]
    
    # yi = 0 #
    Z[which(Y == 0)]  <- rtruncnorm(n0, a = -Inf, b = 0, mean = XB[which(Y == 0)], sd = 1)

    # yi = 1 #
    Z[which(Y == 1)]  <- rtruncnorm(n1, a = 0, b = Inf, mean = XB[which(Y == 1)], sd = 1)
    
    ### sample betas ###
    mub         <- Vb%*%crossprod(X,Z)
    betam[b,]   <- c(rmvn(1, mu = mub, sigma = Vb))
    
    
    if(verbose){
      if(mod(b, dots) == 0){
        cat('.')
      }
      
      if(mod(b, up) == 0){
        cat(paste("\n",b,"samples completed\n"))
      }
    }
  }
  
  betap             <- as.matrix(betam[-c(1:burnin),])
  colnames(betap)   <- colnames(X)
  theta             <- pnorm(t(betap))
  L	                <- cbind(diag(K), -diag(K))
  rhop              <- t(L%*%theta)
  colnames(rhop)    <- paste('rho', 1:K)
  geweke  <- c(geweke.diag(betap)$z)
  out <- list(rho = rhop, betas = betap,
              Z = Z, geweke = geweke, model = model, mcmc = mcmcs)
  
  ### set class ###
  class(out) <- 'bpm'
  
  return(out)
} 

print.bpm		<- function(mod){
  K		<- ncol(mod$model$X)/2
  n		<- nrow(mod$model$X)/(2*K)
  cat("Penalized Bayesian Probit Model for Multivariate Matched Proportions\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probabilities ##
  dmp				<- matrix(apply(mod$rho, 2, median), ncol = K, nrow = 1)
  colnames(dmp)	<- 1:K
  rownames(dmp)	<- ""
  cat("\nMedian difference in marginal probability\n")
  print(round(dmp, digits = 3))
  
}

summary.bpm	<- function(mod, ci.alpha = 0.05, chains = 4){
  K		<- ncol(mod$model$X)/2
  n		<- nrow(mod$model$X)/(2*K)
  B		<- mod$mcmc$B
  burnin	<- mod$mcmc$burnin
  cat("Penalized Bayesian Probit Model for Multivariate Matched Proportions\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probability ##
  rho		<- mod$rho
  low		<- ci.alpha/2
  high	<- 1-(ci.alpha/2)
  tab		<- t(apply(rho, 2, quantile, probs = c(0.5, low, high)))
  prob	<- apply(rho > 0, 2, mean)
  GR		<- getGelman(rho, chains)
  tab1	<- cbind(tab, prob, GR)
  colnames(tab1)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""), "P(d > 0)", 'GR Est.', 'Upper GR')
  rownames(tab1)	<- 1:K
  cat("\nDifference in marginal probability\n\n")
  print(round(tab1, 3))
  cat('Notes: P(r > 0) is posterior probability difference > 0.\nGR denotes Gelman-Rubin, < 1.1 suggests convergence.\n')
  
  cat(paste("\nPosterior samples based on", B, "samples after a burn-in of\n"))
  cat(paste(burnin, "(discarded). "))
}

coef.bpm		<- function(mod, CI = FALSE, ci.alpha = 0.05){
  dmpTable	<- apply(mod$rho, 2, quantile, probs = c(0.5, ci.alpha/2, (1 - ci.alpha/2)))
  
  if(CI){
    dmpt	<- dmpTable
  } else {
    dmpt	<- dmpTable[1, ]
  }
  
  return(dmpt)
}

pprob.bpm	<- function(mod, rho0 = 0){
  probs	<- apply(mod$rho > rho0, 2, mean)
  
  return(probs)
}



#### Penalized Bayesian Probit Model ####
pbpm  <- function(X1, X2, random = FALSE, B = 10000, burnin = B, 
                  Al = 1, Ae = 1, hc = TRUE, al = 0.01, bl = 0.01,
                  verbose = TRUE, up = 2000, dots = 200){
  
  ### data ###
  n		    <- nrow(X1)
  K		    <- ncol(X1)
  ymat		<- as.matrix(cbind(X1, X2))
  bDat	  <- data.frame(Y	= c(t(ymat)), s = rep(1:(2*K), n), id = rep(1:n, each = 2*K))
  form    <- model.frame(Y ~ factor(s) - 1, data = bDat)
  Y       <- model.response(form)
  X       <- model.matrix(form, data = bDat)
  n1		  <- sum(Y)
  n0		  <- length(Y) - n1
  
  model   <- list(Y = Y, X = X, random = random, bdat = bDat)
  mcmcs   <- list(B = B, burnin = burnin)
  if(random){
    U         <- model.matrix(~ as.factor(id) - 1, data = bDat)
    model$U   <- U
  }
  
  # Vb      <- solve(t(X)%*%X)
  XtX     <- t(X)%*%X
  # LtL     <- t(L)%*%L
  if(random){
    UtU     <- t(U)%*%U
  }
  
  # nL      <- ncol(L)
  nX      <- ncol(X)
  if(random){
    nU      <- ncol(U)
  }
  n       <- nrow(X)
  
  ### build matrices ###
  betam   <- matrix(0, nrow = B + burnin, ncol = nX)
  # alpham  <- matrix(0, nrow = B + burnin, ncol = nL)
  if(random){
    omegam  <- matrix(0, nrow = B + burnin, ncol = nU)
    eta     <- rep(1, B + burnin)
    xi      <- rep(1, B + burnin)
  }
  lambda  <- rep(1, B + burnin)
  mu      <- rep(1, B + burnin)
  # gammas  <- matrix(cut.start, nrow = B + burnin, ncol = J-1, byrow = TRUE)
  Z       <- rep(0, n)
  
  ### sampler ###
  for(b in 2:(B + burnin)){
    
    ### sample latent Zs ###
    XB			<- X%*%betam[b-1,] + ifelse(random, U%*%omegam[b-1,], 0)
    
    # yi = 0 #
    Z[which(Y == 0)]  <- rtruncnorm(n0, a = -Inf, b = 0, mean = XB[which(Y == 0)], sd = 1)

    # yi = 1 #
    Z[which(Y == 1)]  <- rtruncnorm(n1, a = 0, b = Inf, mean = XB[which(Y == 1)], sd = 1)
    
    ### sample betas ###
    Pb          <- XtX + (1/lambda[b-1])*diag(2*K)
    Vb          <- solve(Pb)
    mub         <- Vb%*%t(X)%*%(Z - ifelse(random, U%*%omegam[b-1,], 0))
    betam[b,]   <- c(rmvn(1, mu = mub, sigma = Vb))
    
    ### sample lambda, mu ###
    btb         <- t(betam[b-1,])%*%betam[b-1,]
    if(hc == TRUE){
      lambda[b]   <- rinvgamma(1, K + 1/2, 1/mu[b-1] + (1/2)*btb)
      mu[b]       <- rinvgamma(1, 1, 1/(Al^2) + 1/lambda[b-1])
    }
    else{
      lambda[b]   <- rgamma(1, K + al, btb/2 + bl)
    }
      
    ### sample omegas ###
    if(random){
      Po          <- UtU + (1/xi[b-1])*diag(nU)
      Vo          <- solve(Po)
      muo         <- Vo%*%t(U)%*%(Z - X%*%betam[b-1,])
      omegam[b,]  <- c(rmvn(1, mu = muo, sigma = Vo))
    }
    
    ### sample sigma, lambda, and mu ###
    # ata         <- t(alpham[b-1,])%*%alpham[b-1,]
    # if(hc == TRUE){
    #   lambda[b]   <- rinvgamma(1, nL/2 + 1/2, 1/mu[b-1] + (1/2)*ata)
    #   mu[b]       <- rinvgamma(1, 1, 1/(Al^2) + 1/lambda[b-1])
    # } else{
    #   lambda[b]   <- rgamma(1, nL/2 + al, ata/2 + bl)
    # }
    
    ### sample xi and eta ###
    if(random){
      oto         <- t(omegam[b-1,])%*%omegam[b-1,]
      # xi[b]       <- rgamma(1, nU/2 + ax, oto/2 + bx)
      xi[b]       <- rinvgamma(1, nU/2 + 1/2, oto/2 + 1/eta[b-1])
      eta[b]      <- rinvgamma(1, 1, 1/(Ae^2) + 1/xi[b-1])
    }
    
    ### sample gammas ###
    # cuts        <- c(gammas[b-1,], Inf)
    # for(j in 2:length(gammas[b-1,])){
    #   au  <- max(max(Z[Y == Yl[j]]), cuts[j-1])
    #   bu  <- min(min(Z[Y == Yl[j+1]]), cuts[j+1])
    #   gammas[b,j] <- runif(1, au, bu)
    # }
    
    if(verbose){
      if(mod(b, dots) == 0){
        cat('.')
      }
      
      if(mod(b, up) == 0){
        cat(paste("\n",b,"samples completed\n"))
      }
    }
  }
  
  betap             <- as.matrix(betam[-c(1:burnin),])
  colnames(betap)   <- colnames(X)
  theta             <- pnorm(t(betap))
  L	                <- cbind(diag(K), -diag(K))
  rhop              <- t(L%*%theta)
  colnames(rhop)    <- paste('rho', 1:K)
  lamv              <- matrix(lambda[-c(1:burnin)], ncol = 1)
  colnames(lamv)    <- 'lambda'
  if(hc == TRUE){
    muv             <- matrix(mu[-c(1:burnin)], ncol = 1)
    colnames(muv)   <- 'mu'
  } else{
    
  }
  if(random){
    xiv     <- matrix(xi[-c(1:burnin)], ncol = 1)
    colnames(xiv) <- 'xi'
    etav              <- matrix(eta[-c(1:burnin)], ncol = 1)
    colnames(etav)    <- 'eta'
    omegap            <- omegam[-c(1:burnin),]
    colnames(omegap)  <- paste('id', sort(unique(bDat$id)), sep = '')
    geweke  <- c(geweke.diag(betap)$z, 
                 geweke.diag(xiv)$z, geweke.diag(etav)$z)
    out <- list(rho = rhop, betas = betap, lambda = lamv, xi = xiv, eta = etav,
                Z = Z, geweke = geweke, model = model, mcmc = mcmcs)
  } else{
    geweke  <- c(geweke.diag(betap)$z)
    out <- list(rho = rhop, betas = betap, lambda = lamv,
                Z = Z, geweke = geweke, model = model, mcmc = mcmcs)
  }
  
  ### find WAIC ###
  # out$fit     <- WAIC_odlm(out)
  
  ### set class ###
  class(out) <- 'pbpm'
  
  return(out)
} 

print.pbpm		<- function(mod){
  K		<- ncol(mod$model$X)/2
  n		<- nrow(mod$model$X)/(2*K)
  cat("Penalized Bayesian Probit Model for Multivariate Matched Proportions\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probabilities ##
  dmp				<- matrix(apply(mod$rho, 2, median), ncol = K, nrow = 1)
  colnames(dmp)	<- 1:K
  rownames(dmp)	<- ""
  cat("\nMedian difference in marginal probability\n")
  print(round(dmp, digits = 3))
  
}

summary.pbpm	<- function(mod, ci.alpha = 0.05, chains = 4){
  K		<- ncol(mod$model$X)/2
  n		<- nrow(mod$model$X)/(2*K)
  B		<- mod$mcmc$B
  burnin	<- mod$mcmc$burnin
  cat("Penalized Bayesian Probit Model for Multivariate Matched Proportions\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probability ##
  rho		<- mod$rho
  low		<- ci.alpha/2
  high	<- 1-(ci.alpha/2)
  tab		<- t(apply(rho, 2, quantile, probs = c(0.5, low, high)))
  prob	<- apply(rho > 0, 2, mean)
  GR		<- getGelman(rho, chains)
  tab1	<- cbind(tab, prob, GR)
  colnames(tab1)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""), "P(d > 0)", 'GR Est.', 'Upper GR')
  rownames(tab1)	<- 1:K
  cat("\nDifference in marginal probability\n\n")
  print(round(tab1, 3))
  cat('Notes: P(r > 0) is posterior probability difference > 0.\nGR denotes Gelman-Rubin, < 1.1 suggests convergence.\n')
  
  cat(paste("\nPosterior samples based on", B, "samples after a burn-in of\n"))
  cat(paste(burnin, "(discarded). "))
}

coef.pbpm		<- function(mod, CI = FALSE, ci.alpha = 0.05){
  dmpTable	<- apply(mod$rho, 2, quantile, probs = c(0.5, ci.alpha/2, (1 - ci.alpha/2)))
  
  if(CI){
    dmpt	<- dmpTable
  } else {
    dmpt	<- dmpTable[1, ]
  }
  
  return(dmpt)
}

pprob.pbpm	<- function(mod, rho0 = 0){
  probs	<- apply(mod$rho > rho0, 2, mean)
  
  return(probs)
}



#### Bayesian Multivariate Probit Model ####
library(splines)

bmvp  <- function(X1, X2, B = 10000, burnin = B, pen = TRUE, cholesky = FALSE,
                  Al = 1, Ae = 1, hc = FALSE, al = 0.01, bl = 0.01,
                  Kp = 2, alpha = 0.001, As = 1, Bs = 1,
                  verbose = TRUE, up = 2000, dots = 200){
  
  ### data ###
  n		    <- nrow(X1)
  K		    <- ncol(X1)
  ymat		<- as.matrix(cbind(X1, X2))
  bDat	  <- data.frame(Y	= c(t(ymat)), s = rep(1:(2*K), n), id = rep(1:n, each = 2*K))
  form    <- model.frame(Y ~ factor(s) - 1, data = bDat)
  Y       <- model.response(form)
  X       <- model.matrix(form, data = bDat)
  n1		  <- sum(Y)
  n0		  <- length(Y) - n1
  
  model   <- list(Y = Y, X = X, bdat = bDat)
  mcmcs   <- list(B = B, burnin = burnin)
  
  bhat    <- coef(lm(Y ~ X - 1))
  # S       <- (t(ymat)-bhat)%*%t(t(ymat)-bhat)
  S       <- diag(2*K)
  Cf      <- cov2cor(S)
  Cfc     <- chol(Cf)
  Cfk     <- kronecker(diag(n), Cf)
  Ccf     <- chol(Cfk)
  Cci     <- solve(Ccf)
  
  # nL      <- ncol(L)
  nX      <- ncol(X)
  N       <- nrow(X)

  ### build matrices ###
  xi        <- matrix(1, nrow = B + burnin, ncol = nX)
  xi[1,]    <- betas <- qnorm(bhat + 0.0001)
  Zc        <- rep(0, N)
  if(pen == TRUE){
    lambda      <- rep(1, B + burnin)
    mu          <- rep(1, B + burnin)
  }

  sig2.me <- 0.01
  Kt      <- 2*K
  Theta   <- bs(1:(2*K), df = Kt, intercept = TRUE, degree = 3)
  diff0   <- diag(1, 2*K, 2*K)
  diff2   <- matrix(rep(c(1, -2, 1, rep(0, 2*K - 2)), 2*K - 2)[1:((2*K - 2) * 2*K)], 2*K - 2, 2*K, byrow = TRUE)
  P0      <- t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2      <- t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat   <- alpha * P0 + (1 - alpha) * P2
  
  # BPSI = array(NA, c(Kt, Kp, B + burnin))
  bpsi    <- matrix(0, Kt, Kp)
  c.mat   <- matrix(rnorm(n * Kp, 0, 0.01), n, Kp)
  
  lambda.psi <- rep(1, Kp)
  
  psi.cur     <- t(bpsi) %*% t(Theta)
  pcaef.cur   <- c.mat %*% psi.cur
  
  ### sampler ###
  for(b in 2:(B + burnin)){
    
    ### project X ###
    if(cholesky){
      Xc      <- Cci%*%X
      XcB			<- Xc%*%betas
    } else{
      # XcB			<- X%*%betas
      XcB			<- as.vector(X%*%betas) + as.vector(pcaef.cur)
    }
    
    ### sample latent Zs ###
    
    # yi = 0 #
    # Zc[which(Y == 0)]  <- rtnorm(n = n0, lower = -Inf, upper = 0, mean = XcB[which(Y == 0)], sd = 1) #rtmvnorm
    Zc[which(Y == 0)]  <- rtruncnorm(n0, a = -Inf, b = 0, mean = XcB[which(Y == 0)], sd = 1) #rtmvnorm
    # Z[which(Y == 0)]  <- rtmvnorm(n = 1, mean = XB[which(Y == 0)], sigma = Cc[which(Y == 0), which(Y == 0)], 
    #                               lower = rep(-Inf, n0), upper = rep(0, n0), algorithm = "gibbs")
    
    # yi = 1 #
    # Zc[which(Y == 1)]  <- rtnorm(n = n1, lower = 0, upper = Inf, mean = XcB[which(Y == 1)], sd = 1) #rtmvnorm
    Zc[which(Y == 1)]  <- rtruncnorm(n1, a = 0, b = Inf, mean = XcB[which(Y == 1)], sd = 1) #rtmvnorm
    # Z[which(Y == 1)]  <- rtmvnorm(n = 1, mean = XB[which(Y == 1)], sigma = Vc[which(Y == 1), which(Y == 1)], 
    #                               lower = rep(0, n1), upper = rep(Inf, n1), algorithm = "gibbs")
    # Z           <- Ccf%*%Zc
    
    ### sample betas ###
    # if(pen){
    #   XVX         <- t(X)%*%Cfki%*%X + (1/lambda[b-1])*diag(2*K)
    # } else{
    #   XVX         <- t(X)%*%Cfki%*%X
    # }
    # XVXi        <- solve(XVX)
    # XVt         <- crossprod(X,Cfki)
    # mub         <- XVXi%*%XVt%*%Z
    # betam[b,]   <- c(rmvn(1, mu = mub, sigma = XVXi))
    if(cholesky == TRUE){
      if(pen){
        XVX         <- t(Xc)%*%Xc + (1/lambda[b-1])*diag(2*K)
        XVXi        <- solve(XVX)
      } else{
        XVX         <- t(Xc)%*%Xc
        XVXi        <- solve(XVX)
      }
      XtZ         <- crossprod(Xc,Zc)
      mub         <- XVXi%*%XtZ
      betas       <- c(rmvn(1, mu = mub, sigma = XVXi))
      xi[b,]   <- c(Cfc%*%betas)
    } else{
      mean.cur  <- as.vector(pcaef.cur)
      if(pen){
        XVX         <- (1/sig2.me) * t(X)%*%X + (1/lambda[b-1])*diag(2*K)
        XVXi        <- solve(XVX)
      } else{
        XVX         <- (1/sig2.me) * t(X)%*%X
        XVXi        <- solve(XVX)
      }
      XtZ         <- crossprod(X,Zc - mean.cur)
      mub         <- (1/sig2.me) * XVXi%*%XtZ
      betas       <- c(rmvn(1, mu = mub, sigma = XVXi))
      xi[b,]      <- c(betas)
    }
    
    ### sample lambda, mu ###
    if(pen == TRUE){
      btb         <- t(betas)%*%betas
      if(hc == TRUE){
        lambda[b]   <- rinvgamma(1, K + 1/2, 1/mu[b-1] + (1/2)*btb)
        mu[b]       <- rinvgamma(1, 1, 1/(Al^2) + 1/lambda[b-1])
      } else{
        lambda[b]   <- rgamma(1, K + al, btb/2 + bl)
      }
    }
 
    ### Cholesky projection ###
    if(cholesky){
      Zmat        <- matrix(Zc, nrow = n, byrow = TRUE)
      S           <- (t(Zmat)-betas)%*%t(t(Zmat)-betas)
      Cf          <- cov2cor(S)
      Cfc         <- chol(Cf)
      Cfk         <- kronecker(diag(n), Cf)
      Ccf         <- chol(Cfk)
      Cci         <- solve(Ccf)
    } else{
      Zmat        <- matrix(Zc, nrow = n, byrow = TRUE)
      mean.cur    <- as.vector(X%*%betas)
      sigma       <- solve((1/sig2.me) * kronecker(t(c.mat) %*% c.mat, 
                                                   t(Theta) %*% Theta) + kronecker(diag(1/lambda.psi), 
                                                                                   P.mat))
      mup         <- (1/sig2.me) * sigma %*% t(kronecker(c.mat, Theta)) %*% (Zc - mean.cur)
      bpsi        <- matrix(rmvn(1, mu = mup, sigma = sigma), nrow = Kt, ncol = Kp)
      psi.cur     <- t(bpsi) %*% t(Theta)
      ppT         <- psi.cur %*% t(psi.cur)
      for (c in 1:n) {
        sigma       <- solve((1/sig2.me) * ppT + diag(1, Kp, Kp))
        muc         <- (1/sig2.me) * sigma %*% psi.cur %*% (Zmat[c,] - betas)
        c.mat[c, ]  <- rmvn(1, mu = muc, sigma = sigma)
      }
      pcaef.cur   <- c.mat %*% psi.cur
      Z.cur       <- mean.cur + pcaef.cur
      a.post      <- As + n * 2*K/2
      b.post      <- Bs + 1/2 * crossprod(as.vector(Zc - Z.cur))
      sig2.me     <- 1/rgamma(1, a.post, b.post)
      for (kp in 1:Kp) {
        a.post          <- Kt/2 + Kt/2
        b.post          <- Kt/2 + 1/2 * bpsi[, kp] %*% P.mat %*% bpsi[, kp]
        lambda.psi[kp]  <- 1/rgamma(1, a.post, b.post)
      }
    }


    if(verbose){
      if(mod(b, dots) == 0){
        cat('.')
      }
      
      if(mod(b, up) == 0){
        cat(paste("\n",b,"samples completed\n"))
      }
    }
  }
  
  xip               <- as.matrix(xi[-c(1:burnin),])
  colnames(xip)     <- colnames(X)
  theta             <- pnorm(t(xip))
  L	                <- cbind(diag(K), -diag(K))
  rhop              <- t(L%*%theta)
  colnames(rhop)    <- paste('rho', 1:K)
  geweke  <- c(geweke.diag(xip)$z)
  out <- list(rho = rhop, xis = xip,
              Zc = Zc, geweke = geweke, model = model, mcmc = mcmcs)
  
  ### set class ###
  class(out) <- 'bmvp'
  
  return(out)
}

print.bmvp		<- function(mod){
  K		<- ncol(mod$model$X)/2
  n		<- nrow(mod$model$X)/(2*K)
  cat("Bayesian Multivariate Probit Model for Multivariate Matched Proportions\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probabilities ##
  dmp				<- matrix(apply(mod$rho, 2, median), ncol = K, nrow = 1)
  colnames(dmp)	<- 1:K
  rownames(dmp)	<- ""
  cat("\nMedian difference in marginal probability\n")
  print(round(dmp, digits = 3))
  
}

summary.bmvp	<- function(mod, ci.alpha = 0.05, chains = 4){
  K		<- ncol(mod$model$X)/2
  n		<- nrow(mod$model$X)/(2*K)
  B		<- mod$mcmc$B
  burnin	<- mod$mcmc$burnin
  cat("Bayesian Multivariate Probit Model for Multivariate Matched Proportions\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probability ##
  rho		<- mod$rho
  low		<- ci.alpha/2
  high	<- 1-(ci.alpha/2)
  tab		<- t(apply(rho, 2, quantile, probs = c(0.5, low, high)))
  prob	<- apply(rho > 0, 2, mean)
  GR		<- getGelman(rho, chains)
  tab1	<- cbind(tab, prob, GR)
  colnames(tab1)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""), "P(d > 0)", 'GR Est.', 'Upper GR')
  rownames(tab1)	<- 1:K
  cat("\nDifference in marginal probability\n\n")
  print(round(tab1, 3))
  cat('Notes: P(r > 0) is posterior probability difference > 0.\nGR denotes Gelman-Rubin, < 1.1 suggests convergence.\n')
  
  cat(paste("\nPosterior samples based on", B, "samples after a burn-in of\n"))
  cat(paste(burnin, "(discarded). "))
}

coef.bmvp		<- function(mod, CI = FALSE, ci.alpha = 0.05){
  dmpTable	<- apply(mod$rho, 2, quantile, probs = c(0.5, ci.alpha/2, (1 - ci.alpha/2)))
  
  if(CI){
    dmpt	<- dmpTable
  } else {
    dmpt	<- dmpTable[1, ]
  }
  
  return(dmpt)
}

pprob.bmvp	<- function(mod, rho0 = 0){
  probs	<- apply(mod$rho > rho0, 2, mean)
  
  return(probs)
}


#### ERM, Liu and Chang 2016 ####
erm	<- function(X1, X2, exact = TRUE, alpha = 0.05, direction = "two.sided"){
  K		<- ncol(X1)
  counts	<- getCounts(X1, X2, K)
  n     <- nrow(X1)
  
  n21s  <- counts$n21 + 0.5
  n12s  <- counts$n12 + 0.5
  
  est		<- n21s/n - n12s/n
  sdv		<- sqrt(1/n21s + 1/n12s)
  
  pval	<- vector('numeric', length = K)
  int		<- matrix(0, nrow = K, ncol = 2)
  for(k in 1:K){
    if(exact){
      n21k	<- counts$n21[k]
      n12k	<- counts$n12[k]
      temp	<- binom.test(n21k, n21k + n12k, p = 0.5, alternative = direction, conf.level = 1 - alpha)
      pval[k]	<- temp$p.value
      int[k,]	<- log(temp$conf.int/(1-temp$conf.int))
    } else {
      int[k,]	<- c(est[k] - qnorm(1-alpha/2)*sdv[k], est[k] + qnorm(1-alpha/2)*sdv[k])
      if(direction == "two.sided"){
        pval[k]	<- 1 - pchisq((est[k]/sdv[k])^2, df = 1)
      } else if(direction == "less"){
        pval[k]	<- pnorm(est[k]/sdv[k])
      } else {
        pval[k]	<- 1 - pnorm(est[k]/sdv[k])				
      }
    }	
  }
  
  l		<- list(est = est, sdv = sdv, int = int, pval = pval, data = list(X1 = X1, X2 = X2), args = list(exact = exact, alpha = alpha, direction = direction))
  
  class(l)	<- 'erm'
  
  return(l)
}

print.erm	<- function(mod){
  K		<- ncol(mod$data$X1)
  n		<- nrow(mod$data$X1)
  cat("CMH-based multivariate matched proportions\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in log odds ratios probabilities ##
  dmp				<- matrix(c(mod$est), ncol = K, nrow = 1)
  colnames(dmp)	<- 1:K
  rownames(dmp)	<- ""
  cat("\nLog odds ratio\n")
  print(round(dmp, digits = 3))
  
}

summary.erm	<- function(mod){
  K		<- ncol(mod$data$X1)
  n		<- nrow(mod$data$X1)
  exact	<- mod$args$exact
  cat("CMH-based multivariate matched proportions\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## odds ratio ##
  if(exact){
    rho		<- mod$est
    sdv		<- mod$sdv
    pval	<- mod$pval
    tab1	<- cbind(rho, sdv, pval)
    colnames(tab1)	<- c('Estimate', 'Std.err', 'Exact P')
    rownames(tab1)	<- 1:K
    cat("\nLog odds ratio\n\n")
    print(round(tab1, 3))
    cat('Notes: Estimates based on Lui and Chang (2016).\nStd.err estimated using aysmptotic estimator.\n')	
  } else {
    rho		<- mod$est
    sdv		<- mod$sdv
    z		<- rho/sdv
    pval	<- mod$pval
    tab1	<- cbind(rho, sdv, z, pval)
    colnames(tab1)	<- c('Estimate', 'Std.err', 'Z', 'Pr(>|Z|)')
    rownames(tab1)	<- 1:K
    cat("\nLog odds ratio\n\n")
    print(round(tab1, 3))
    cat('Notes: Estimates based on Lui and Chang (2016).\nStd.err estimated using aysmptotic estimator.\n')	
  }
  
  
}

coef.erm		<- function(mod){
  rho	<- c(mod$est)
  
  return(rho)
}



#### GEE, Klingenberg and Agresti (2006) ####
geemmp	<- function(X1, X2, family = gaussian(link = 'identity'), corstr = 'independence'){
  n		<- nrow(X1)
  K		<- ncol(X1)
  Y		<- as.matrix(cbind(X1, X2))
  gDat	<- data.frame(y	= c(t(Y)), s = rep(1:(2*K), n), id = rep(1:n, each = 2*K))
  
  model	<- geeglm(y ~ factor(s) - 1, family = family, data = gDat, id = id, corstr = corstr)
  
  beta	<- coef(model)
  covb	<- summary(model)$cov.scaled
  # covb	<- (t(Y - beta)%*%(Y - beta))/(n^2)
  
  
  # tests of differences #
  L	      <- cbind(diag(K), -diag(K))
  rho		  <- L%*%beta
  seRho   <- sqrt(diag(L%*%covb%*%t(L)))
  lower	  <- rho - qnorm(0.975)*seRho
  upper	  <- rho + qnorm(0.975)*seRho
  wald	  <- (rho/seRho)^2
  pval	  <- 1 - pchisq(wald, 1)
  
  ## Simultaneous Marginal Homogeneity ##
  # KA Generalized Score Test #
  zse     <- which(diag(covb) == 0)
  # 	if(length(zse) == 0){
  #   	Tgs		<- (n^2)*t(rho)%*%solve(L%*%(t(Y)%*%Y)%*%t(L))%*%rho
  #   	pTgs	<- 1 - pchisq(Tgs, nrow(L))
  # 	} else {
  Tgs		<- NaN
  pTgs	<- NaN
  # }
  
  ## if estimating 
  if(corstr == 'unstructured'){
    cory	<- matrix(0, 2*K, 2*K)
    cory[lower.tri(cory, diag = FALSE)]	<- summary(model)$corr[,1]
    cory[upper.tri(cory, diag = FALSE)]	<- t(cory)[upper.tri(t(cory), diag = FALSE)]		
    diag(cory)	<- rep(1, 2*K)
  } else {
    cory	<- NULL
  }
  
  l		<- list(rho = rho, beta = beta, sdv = seRho, cov = list(covb = covb, cory = cory, zse = zse), lower = lower, upper = upper, wald = wald, pval = pval, model = model, smh = list(Tgs = Tgs, pTgs = pTgs), data = list(X1 = X1, X2 = X2))
  
  class(l)	<- 'gbmp'
  
  # 	if(zse > 0){
  #   	warning('Sparse columns detected, variance estimates may be singular')
  # 	}
  
  return(l)
}

print.gbmp	<- function(mod){
  K		<- ncol(mod$data$X1)
  n		<- nrow(mod$data$X1)
  cat("GEE-based multivariate matched proportions model\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probabilities ##
  dmp				<- matrix(c(mod$rho), ncol = K, nrow = 1)
  colnames(dmp)	<- 1:K
  rownames(dmp)	<- ""
  cat("\nDifference in marginal probability\n")
  print(round(dmp, digits = 3))
  
}

summary.gbmp	<- function(mod){
  K		<- ncol(mod$data$X1)
  n		<- nrow(mod$data$X1)
  cat("GEE-based multivariate matched proportions model\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## odds ratio ##
  rho		<- mod$rho
  sdv		<- mod$sdv
  wald	<- mod$wald
  pval	<- mod$pval
  tab1	<- cbind(rho, sdv, wald, pval)
  colnames(tab1)	<- c('Estimate', 'Std.err', 'Wald', 'Pr(>|W|)')
  rownames(tab1)	<- 1:K
  cat("\nDifference in marginal probability\n\n")
  print(round(tab1, 3))
  cat(paste('Notes: Estimates based on GEE with family: ', mod$model$family$family, ', link: ', mod$model$family$link, '.\nStd.err estimated using robust estimator.\n', sep = ''))
  
  cat("\nTest of simultaneous marginal homogeneity\n")
  cat(paste("Score statistic:", round(mod$smh$Tgs, 3), "on", K, "df, p-value =", round(mod$smh$pTgs, 3)))
  
}

coef.gbmp		<- function(mod){
  dmpt	<- c(mod$rho)
  
  return(dmpt)
}



#### Boot, Westfall, Troendle, and Pennello (2010) ####
bootmmp	<- function(X1, X2, B, int.hyp = FALSE){
  X		<- cbind(X1, X2)
  K		<- ncol(X1)
  n		<- nrow(X1)
  
  L		<- cbind(diag(K), -diag(K))
  
  D		<- t(L%*%t(X))
  Dbar	<- apply(D, 2, mean)
  S		<- diag(apply((t(t(D) - Dbar))^2, 2, mean))
  if(length(which(diag(S) == 0)) > 0){
    Z <- NaN
  } else{
    Z		<- sqrt(n)*solve(sqrt(S))%*%Dbar
  }
  
  Dboot	<- matrix(0, nrow = B, ncol = K)
  if(int.hyp){
    Zboot	<- matrix(0, nrow = B, ncol = K)
    Zmax	<- vector('numeric', length = B)
  }
  
  for(b in 1:B){
    idb			<- sample(1:n, n, replace = TRUE)
    Ds			<- D[idb,]
    Dbb			<- Dboot[b,] <- apply(Ds, 2, mean)
    if(int.hyp){
      Sb			<- sqrt(apply((t(t(Ds) - Dbb))^2, 2, mean))
      Zboot[b,]	<- c(sqrt(n)*ifelse(Sb > 0, Dbb/Sb, 0))
      Zmax[b]		<- max(abs(Zboot[b,]))
    }	
  }
  
  resMat	<- cbind(t(apply(Dboot, 2, quantile, probs = c(0.5, 0.025, 0.975))), apply(Dboot > 0, 2, mean))
  
  if(int.hyp){
    l	<- list(resMat = resMat, Dboot = Dboot, Zboot = Zboot, Zmax = Zmax, Z = Z, S = S, Dbar = Dbar, data = list(X1 = X1, X2 = X2))	
  } else {
    l	<- list(resMat = resMat, Dboot = Dboot, data = list(X1 = X1, X2 = X2))
  }
  
  class(l)	<- 'bootmmp'
  
  return(l)
}

print.bootmmp	<- function(mod){
  K		<- ncol(mod$data$X1)
  n		<- nrow(mod$data$X1)
  cat("Bootstraped multivariate matched proportions model\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probabilities ##
  dmp				<- matrix(c(mod$resMat[,1]), ncol = K, nrow = 1)
  colnames(dmp)	<- 1:K
  rownames(dmp)	<- ""
  cat("\nDifference in marginal probability\n")
  print(round(dmp, digits = 3))
  
}

summary.bootmmp	<- function(mod, ci.alpha = 0.05){
  K		<- ncol(mod$data$X1)
  n		<- nrow(mod$data$X1)
  cat("Bootstraped multivariate matched proportions model\n")
  cat(paste("K =", K, "outcomes for n =", n, "pairs \n"))
  
  ## difference in marginal probabilities ##
  low				<- ci.alpha/2
  high			<- 1-(ci.alpha/2)
  tab1			<- mod$resMat
  colnames(tab1)	<- c('Median', paste(100*low,"%", sep = ""), paste(100*high,"%", sep = ""), "P(d > 0)")
  rownames(tab1)	<- 1:K
  cat("\nDifference in marginal probability\n\n")
  print(round(tab1, 3))
  
}

coef.bootmmp	<- function(mod){
  dmpt	<- c(mod$resMat[,1])
  
  return(dmpt)
}


#### data ####

soc <- rbind(matrix(rep(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1), times = 17), nrow = 17, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0), times = 16), nrow = 16, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0), times = 5), nrow = 5, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 1), times = 5), nrow = 5, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 1), times = 3), nrow = 3, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0), times = 3), nrow = 3, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 0, 0, 1, 0, 0, 1), times = 3), nrow = 3, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 1, 0, 0, 0, 0, 0, 0), times = 2), nrow = 2, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0), times = 2), nrow = 2, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1), times = 2), nrow = 2, byrow = TRUE),
             matrix(rep(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 1, 0, 0, 0, 1, 1, 1, 0), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 0, 1, 0, 1, 0, 1, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 1, 0, 0, 1, 0, 1, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 1, 0, 1, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 1, 0, 1, 0, 1, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 1, 0, 0, 1, 1, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(1, 0, 0, 0, 1, 0, 0, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 1, 1, 0, 1, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(1, 0, 0, 0, 0, 0, 1, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 1, 0, 0, 1, 0, 0, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 0, 0, 1, 0, 0, 1, 0, 1, 1), times = 1), nrow = 1, byrow = TRUE),
             matrix(rep(c(0, 1, 0, 1, 0, 0, 1, 1, 1, 1), times = 1), nrow = 1, byrow = TRUE))
