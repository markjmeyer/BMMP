#### Source Code ####
source('bmmp.R')

#### Manuscript Simulation ####

##### Settings #####
n       <- 75 # 150
K       <- 2
M       <- 500
sCov    <- (0.01^2)*matrix(c(1, 0.5, 0, 0,
                             0.5, 1, 0, 0,
                             0, 0, 1, 0.5,
                             0, 0, 0.5, 1), nrow = 2*K, byrow = TRUE)
theta12   <- seq(0.05, 0.2, by = 0.01)

##### BSpAM specs #####
B         <- 10000
burnin    <- B
dots      <- 50
up        <- 1000

##### storage #####
br1i75s  <- array(0, dim = c(M, 4, ifelse(sparse, length(theta12), length(theta21))))
br2i75s  <- array(0, dim = c(M, 4, ifelse(sparse, length(theta12), length(theta21))))

pw1i75s  <- array(0, dim = c(M, 4, ifelse(sparse, length(theta12), length(theta21))))
pw2i75s  <- array(0, dim = c(M, 4, ifelse(sparse, length(theta12), length(theta21))))

wd1i75s  <- array(0, dim = c(M, 4, ifelse(sparse, length(theta12), length(theta21))))
wd2i75s  <- array(0, dim = c(M, 4, ifelse(sparse, length(theta12), length(theta21))))

cr1i75s  <- array(0, dim = c(M, 4, ifelse(sparse, length(theta12), length(theta21))))
cr2i75s  <- array(0, dim = c(M, 4, ifelse(sparse, length(theta12), length(theta21))))

spci75s  <- array(0, dim = c(M, 2, ifelse(sparse, length(theta12), length(theta21))))

##### simulation run #####
iter <- 0

for(d in 1:ifelse(sparse, length(theta12), length(theta21))){
  tb    <- c(0.05, theta12[d],
              0.25, 0.005)
  for(m in 1:M){
    set.seed(m)
    Ymc   <- rmvnorm(n, mean = tb, sigma = sCov)
    Ymcl  <- c(t(Ymc))
    Yb    <- rbinom(length(Ymcl), 1, prob = Ymcl*(Ymcl > 0)) # abs or set to zero?
    Ys    <- matrix(Yb, nrow = n, byrow = TRUE)
    X1s   <- Ys[,1:K]
    X2s   <- Ys[,(K+1):(2*K)]
    X 	  <- matrix(unlist(getCounts(X1s, X2s, K)), ncol = K, byrow = TRUE)
    
    if(apply(X1s, 2, sum)[2] == 0 & apply(X2s, 2, sum)[2] == 0){
      next
    }
    Xs1Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1s, X2s, K), function(x) x == 0)), 
                                    nrow = 4, byrow = TRUE), 2, sum) > 0)
    Xs2Cols   <- which(apply(matrix(unlist(lapply(getCounts(X1s, X2s, K), function(x) x == 0)), 
                                    nrow = 4, byrow = TRUE), 2, sum) > 0) + K
    XsCols	<- c(Xs1Cols, Xs2Cols)
    
    ## sparsity checks ##
    spci75s[m,1,d]  <- length(XsCols) > 0
    spci75s[m,2,d]  <- sum(apply(Ys, 2, sum) == 0) > 0
    if(sum(apply(Ys, 2, sum) == 0) == 0){
      next
    }
    
    ## run models ##
    bmvpp   <- bmvp(X1s, X2s, B = B, burnin = burnin, pen = TRUE, cholesky = FALSE,
                    Al = 1, Ae = 1, hc = TRUE, al = 0.01, bl = 0.01, 
                    Kp = 2, alpha = 0.001, As = 1, Bs = 1,
                    verbose = FALSE)
    geek    <- geemmp(X1s, X2s)
    boot    <- bootmmp(X1s, X2s, B = 10000)
    lcr121	<- ifelse(X[1,1] == 0, (X[1,1] + 0.5)/(sum(X[,1] + 0.5)), X[1,1]/sum(X[,1]))
    lcr112	<- ifelse(X[2,1] == 0, (X[2,1] + 0.5)/(sum(X[,1] + 0.5)), X[2,1]/sum(X[,1]))
    lcr221	<- ifelse(X[1,2] == 0, (X[1,2] + 0.5)/(sum(X[,2] + 0.5)), X[1,2]/sum(X[,2]))
    lcr212	<- ifelse(X[2,2] == 0, (X[2,2] + 0.5)/(sum(X[,2] + 0.5)), X[2,2]/sum(X[,2]))
    
    lcRho1 	<- lcr121 - lcr112
    lcRho2 	<- lcr221 - lcr212
    
    ## extract results ###
    ## bias ##
    # rho #
    br1i75s[m,,d]   <- (tb[1] - tb[3]) - c(median(bmvpp$rho[,1]), lcRho1, geek$rho[1], boot$resMat[1,1])
    br2i75s[m,,d]   <- (tb[2] - tb[4]) - c(median(bmvpp$rho[,2]), lcRho2, geek$rho[2], boot$resMat[2,1])
    
    ## intervals ##
    # estimates #
    bmvppi1     <- quantile(bmvpp$rho[,1], probs = c(0.025, 0.975))
    geeki1      <- c(geek$lower[1], geek$upper[1])
    booti1      <- boot$resMat[1,2:3]
    
    bmvppi2     <- quantile(bmvpp$rho[,2], probs = c(0.025, 0.975))
    geeki2      <- c(geek$lower[2], geek$upper[2])
    booti2      <- boot$resMat[2,2:3]
    
    oldw 			<- getOption("warn")
    options(warn = -1)
    lcri1 		<- prop.test(n*c(lcr121, lcr112), c(n,n), correct = FALSE)
    
    lcri2 		<- prop.test(n*c(lcr221, lcr212), c(n,n), correct = FALSE)
    options(warn = oldw)
    
    # power #
    pw1i75s[m,,d]   <- 1*c((prod(bmvppi1) > 0), (lcri1$p.val < 0.05), (geek$pval[1] < 0.05), (prod(booti1) > 0))
    pw2i75s[m,,d]   <- 1*c((prod(bmvppi2) > 0), (lcri2$p.val < 0.05), (geek$pval[2] < 0.05), (prod(booti2) > 0))
    
    # width #
    wd1i75s[m,,d]   <- c(diff(bmvppi1), diff(lcri1$conf), diff(geeki1), diff(booti1))
    wd2i75s[m,,d]   <- c(diff(bmvppi2), diff(lcri2$conf), diff(geeki2), diff(booti2))
    
    # coverage #
    cr1i75s[m,,d]   <- 1*c((bmvppi1[1] < tb[1] - tb[3] & bmvppi1[2] > tb[1] - tb[3]),
                           ((lcri1$conf[1] < (tb[1] - tb[3])) & (lcri1$conf[2] > (tb[1] - tb[3]))),
                           (geeki1[1] < tb[1] - tb[3] & geeki1[2] > tb[1] - tb[3]),
                           (booti1[1] < tb[1] - tb[3] & booti1[2] > tb[1] - tb[3]))
    cr2i75s[m,,d]   <- 1*c((bmvppi2[1] < tb[2] - tb[4] & bmvppi2[2] > tb[2] - tb[4]),
                           ((lcri2$conf[1] < (tb[2] - tb[4])) & (lcri2$conf[2] > (tb[2] - tb[4]))),
                           (geeki2[1] < tb[2] - tb[4] & geeki2[2] > tb[2] - tb[4]),
                           (booti2[1] < tb[2] - tb[4] & booti2[2] > tb[2] - tb[4]))
    
    ## simulation controls ##
    iter <- iter + 1
    
    if(mod(iter, dots) == 0){
      cat('.')
    }
    
    if(mod(iter, up) == 0){
      cat(paste("\n",iter,"datasets completed\n"))
    }
  }
}

##### post-process #####
bias1   <- matrix(0, nrow = 4, ncol = ifelse(sparse, length(theta12), length(theta21)))
bias2   <- matrix(0, nrow = 4, ncol = ifelse(sparse, length(theta12), length(theta21)))

mse1    <- matrix(0, nrow = 4, ncol = ifelse(sparse, length(theta12), length(theta21)))
mse2    <- matrix(0, nrow = 4, ncol = ifelse(sparse, length(theta12), length(theta21)))

power1  <- matrix(0, nrow = 4, ncol = ifelse(sparse, length(theta12), length(theta21)))
power2  <- matrix(0, nrow = 4, ncol = ifelse(sparse, length(theta12), length(theta21)))

width1  <- matrix(0, nrow = 4, ncol = ifelse(sparse, length(theta12), length(theta21)))
width2  <- matrix(0, nrow = 4, ncol = ifelse(sparse, length(theta12), length(theta21)))

cover1  <- matrix(0, nrow = 4, ncol = ifelse(sparse, length(theta12), length(theta21)))
cover2  <- matrix(0, nrow = 4, ncol = ifelse(sparse, length(theta12), length(theta21)))

apply(spci75s[,2,], 2, sum) # make sure these are larger than value in line 159

for(i in 1:length(theta12)){
  set.seed(2022)
  sids      <- which(spci75s[,2,i] == 1)
  ids       <- sample(sids, size = 200) # make sure '200' is smaller than min of line 154
  
  b1i75s   <- apply(abs(br1i75s[ids,,i]), 2, mean)
  b2i75s   <- apply(abs(br2i75s[ids,,i]), 2, mean)
  
  m1i75s    <- apply((br1i75s[ids,,i])^2, 2, mean)
  m2i75s    <- apply((br2i75s[ids,,i])^2, 2, mean)
  
  p1i75s   <- apply(pw1i75s[ids,,i], 2, mean, na.rm = TRUE)
  p2i75s   <- apply(pw2i75s[ids,,i], 2, mean)
  
  w1i75s   <- apply(wd1i75s[ids,,i], 2, mean, na.rm = TRUE)
  w2i75s   <- apply(wd2i75s[ids,,i], 2, mean)
  
  c1i75s   <- apply(cr1i75s[ids,,i], 2, mean)
  c2i75s   <- apply(cr2i75s[ids,,i], 2, mean)
  
  bias1[,i]   <- b1i75s
  bias2[,i]   <- b2i75s
  
  mse1[,i]    <- m1i75s
  mse2[,i]    <- m2i75s
  
  power1[,i]  <- p1i75s
  power2[,i]  <- p2i75s
  
  width1[,i]  <- w1i75s
  width2[,i]  <- w2i75s
  
  cover1[,i]  <- c1i75s
  cover2[,i]  <- c2i75s
}

#### save output ####
# save.image(file = 'sim_K2_N75_mvp_bc.RData')

