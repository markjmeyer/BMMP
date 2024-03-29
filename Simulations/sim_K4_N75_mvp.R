#### Source Code ####
source('bmmp.R')

#### Manuscript Simulation ####

##### Settings #####
n       <- 75
K       <- 4
M       <- 500 # 500
mu      <- 0.5
sCov    <- (0.01^2)*diag(2*K)
theta12   <- seq(0.05, 0.2, by = 0.01)

##### BSpAM specs #####
prt    <- list(at = 1, bt = 1, A = 5, ae = 0.001, be = 0.001)
B         <- 10000
burnin    <- B
dots      <- 50
up        <- 1000

##### storage #####
br1i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
br2i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
br3i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
br4i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))

pw1i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
pw2i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
pw3i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
pw4i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))

wd1i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
wd2i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
wd3i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
wd4i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))

cr1i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
cr2i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
cr3i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))
cr4i75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))

spci75s  <- array(0, dim = c(M, 5, ifelse(sparse, length(theta12), length(theta21))))

##### simulation run #####
iter <- 0

for(d in 1:ifelse(sparse, length(theta12), length(theta21))){
  tb    <- c(0.05, theta12[d], 0.005, 0.1,
             0.25, 0.005, 0.05, 0.15)
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
    
    ## run models ##
    bmvpn   <- bmvp(X1s, X2s, B = B, burnin = burnin, pen = TRUE, cholesky = FALSE,
                    Al = 1, Ae = 1, hc = FALSE, al = 0.01, bl = 0.01,
                    Kp = 2, alpha = 0.001, As = 1, Bs = 1,
                    verbose = FALSE)
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
    lcr321	<- ifelse(X[1,3] == 0, (X[1,3] + 0.5)/(sum(X[,3] + 0.5)), X[1,3]/sum(X[,3]))
    lcr312	<- ifelse(X[2,3] == 0, (X[2,3] + 0.5)/(sum(X[,3] + 0.5)), X[2,3]/sum(X[,3]))
    lcr421	<- ifelse(X[1,4] == 0, (X[1,4] + 0.5)/(sum(X[,4] + 0.5)), X[1,4]/sum(X[,4]))
    lcr412	<- ifelse(X[2,4] == 0, (X[2,4] + 0.5)/(sum(X[,4] + 0.5)), X[2,4]/sum(X[,4]))
    
    lcRho1 	<- lcr121 - lcr112
    lcRho2 	<- lcr221 - lcr212
    lcRho3 	<- lcr321 - lcr312
    lcRho4 	<- lcr421 - lcr412
    
    ## extract results ###
    ## bias ##
    # rho #
    br1i75s[m,,d]   <- (tb[1] - tb[1 + K]) - c(median(bmvpn$rho[,1]), median(bmvpp$rho[,1]), lcRho1, geek$rho[1], boot$resMat[1,1])
    br2i75s[m,,d]   <- (tb[2] - tb[2 + K]) - c(median(bmvpn$rho[,2]), median(bmvpp$rho[,2]), lcRho2, geek$rho[2], boot$resMat[2,1])
    br3i75s[m,,d]   <- (tb[3] - tb[3 + K]) - c(median(bmvpn$rho[,3]), median(bmvpp$rho[,3]), lcRho3, geek$rho[3], boot$resMat[3,1])
    br4i75s[m,,d]   <- (tb[4] - tb[4 + K]) - c(median(bmvpn$rho[,4]), median(bmvpp$rho[,4]), lcRho4, geek$rho[4], boot$resMat[4,1])
    
    ## intervals ##
    # estimates #
    bmvpni1    <- quantile(bmvpn$rho[,1], probs = c(0.025, 0.975))
    bmvppi1     <- quantile(bmvpp$rho[,1], probs = c(0.025, 0.975))
    geeki1      <- c(geek$lower[1], geek$upper[1])
    booti1      <- boot$resMat[1,2:3]
    
    bmvpni2    <- quantile(bmvpn$rho[,2], probs = c(0.025, 0.975))
    bmvppi2     <- quantile(bmvpp$rho[,2], probs = c(0.025, 0.975))
    geeki2      <- c(geek$lower[2], geek$upper[2])
    booti2      <- boot$resMat[2,2:3]
    
    bmvpni3    <- quantile(bmvpn$rho[,3], probs = c(0.025, 0.975))
    bmvppi3     <- quantile(bmvpp$rho[,3], probs = c(0.025, 0.975))
    geeki3      <- c(geek$lower[3], geek$upper[3])
    booti3      <- boot$resMat[3,2:3]
    
    bmvpni4    <- quantile(bmvpn$rho[,4], probs = c(0.025, 0.975))
    bmvppi4     <- quantile(bmvpp$rho[,4], probs = c(0.025, 0.975))
    geeki4      <- c(geek$lower[4], geek$upper[4])
    booti4      <- boot$resMat[4,2:3]
    
    oldw 			<- getOption("warn")
    options(warn = -1)
    lcri1 		<- prop.test(n*c(lcr121, lcr112), c(n,n), correct = FALSE)
    
    lcri2 		<- prop.test(n*c(lcr221, lcr212), c(n,n), correct = FALSE)
    
    lcri3 		<- prop.test(n*c(lcr321, lcr312), c(n,n), correct = FALSE)
    
    lcri4 		<- prop.test(n*c(lcr421, lcr412), c(n,n), correct = FALSE)
    options(warn = oldw)
    
    # power #
    pw1i75s[m,,d]   <- 1*c((prod(bmvpni1) > 0), (prod(bmvppi1) > 0), (lcri1$p.val < 0.05), (geek$pval[1] < 0.05), (prod(booti1) > 0))
    pw2i75s[m,,d]   <- 1*c((prod(bmvpni2) > 0), (prod(bmvppi2) > 0), (lcri2$p.val < 0.05), (geek$pval[2] < 0.05), (prod(booti2) > 0))
    pw3i75s[m,,d]   <- 1*c((prod(bmvpni3) > 0), (prod(bmvppi3) > 0), (lcri3$p.val < 0.05), (geek$pval[3] < 0.05), (prod(booti3) > 0))
    pw4i75s[m,,d]   <- 1*c((prod(bmvpni4) > 0), (prod(bmvppi4) > 0), (lcri4$p.val < 0.05), (geek$pval[4] < 0.05), (prod(booti4) > 0))
    
    # width #
    wd1i75s[m,,d]   <- c(diff(bmvpni1), diff(bmvppi1), diff(lcri1$conf), diff(geeki1), diff(booti1))
    wd2i75s[m,,d]   <- c(diff(bmvpni2), diff(bmvppi2), diff(lcri2$conf), diff(geeki2), diff(booti2))
    wd3i75s[m,,d]   <- c(diff(bmvpni3), diff(bmvppi3), diff(lcri3$conf), diff(geeki3), diff(booti3))
    wd4i75s[m,,d]   <- c(diff(bmvpni4), diff(bmvppi4), diff(lcri4$conf), diff(geeki4), diff(booti4))
    
    # coverage #
    cr1i75s[m,,d]   <- 1*c((bmvpni1[1] < tb[1] - tb[1 + K] & bmvpni1[2] > tb[1] - tb[1 + K]),
                           (bmvppi1[1] < tb[1] - tb[1 + K] & bmvppi1[2] > tb[1] - tb[1 + K]),
                           ((lcri1$conf[1] < (tb[1] - tb[1 + K])) & (lcri1$conf[2] > (tb[1] - tb[1 + K]))),
                           (geeki1[1] < tb[1] - tb[1 + K] & geeki1[2] > tb[1] - tb[1 + K]),
                           (booti1[1] < tb[1] - tb[1 + K] & booti1[2] > tb[1] - tb[1 + K]))
    
    cr2i75s[m,,d]   <- 1*c((bmvpni2[1] < tb[2] - tb[2 + K] & bmvpni2[2] > tb[2] - tb[2 + K]),
                           (bmvppi2[1] < tb[2] - tb[2 + K] & bmvppi2[2] > tb[2] - tb[2 + K]),
                           ((lcri2$conf[1] < (tb[2] - tb[2 + K])) & (lcri2$conf[2] > (tb[2] - tb[2 + K]))),
                           (geeki2[1] < tb[2] - tb[2 + K] & geeki2[2] > tb[2] - tb[2 + K]),
                           (booti2[1] < tb[2] - tb[2 + K] & booti2[2] > tb[2] - tb[2 + K]))
    
    cr3i75s[m,,d]   <- 1*c((bmvpni3[1] < tb[3] - tb[3 + K] & bmvpni3[2] > tb[3] - tb[3 + K]),
                           (bmvppi3[1] < tb[3] - tb[3 + K] & bmvppi3[2] > tb[3] - tb[3 + K]),
                           ((lcri3$conf[1] < (tb[3] - tb[3 + K])) & (lcri3$conf[2] > (tb[3] - tb[3 + K]))),
                           (geeki3[1] < tb[3] - tb[3 + K] & geeki3[2] > tb[3] - tb[3 + K]),
                           (booti3[1] < tb[3] - tb[3 + K] & booti3[2] > tb[3] - tb[3 + K]))
    
    cr4i75s[m,,d]   <- 1*c((bmvpni4[1] < tb[4] - tb[4 + K] & bmvpni4[2] > tb[4] - tb[4 + K]),
                           (bmvppi4[1] < tb[4] - tb[4 + K] & bmvppi4[2] > tb[4] - tb[4 + K]),
                           ((lcri4$conf[1] < (tb[4] - tb[4 + K])) & (lcri4$conf[2] > (tb[4] - tb[4 + K]))),
                           (geeki4[1] < tb[4] - tb[4 + K] & geeki4[2] > tb[4] - tb[4 + K]),
                           (booti4[1] < tb[4] - tb[4 + K] & booti4[2] > tb[4] - tb[4 + K]))
    
    ## sparsity checks ##
    spci75s[m,1,d]  <- length(XsCols) > 0
    spci75s[m,2,d]  <- sum(apply(Ys, 2, sum) == 0) > 0
    
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
bias1   <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
bias2   <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
bias3   <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
bias4   <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))

mse1    <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
mse2    <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
mse3    <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
mse4    <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))

power1  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
power2  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
power3  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
power4  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))

width1  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
width2  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
width3  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
width4  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))

cover1  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
cover2  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
cover3  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))
cover4  <- matrix(0, nrow = 5, ncol = ifelse(sparse, length(theta12), length(theta21)))

apply(spci75s[,2,], 2, sum) # make sure these are more than value in line 220

for(i in 1:length(theta12)){
  set.seed(2022)
  sids      <- which(spci75s[,2,i] == 1)
  ids       <- sample(sids, size = 200) # make sure '200' is smaller than min of line 215
  
  b1i75s   <- apply(abs(br1i75s[ids,,i]), 2, mean)
  b2i75s   <- apply(abs(br2i75s[ids,,i]), 2, mean)
  b3i75s   <- apply(abs(br3i75s[ids,,i]), 2, mean)
  b4i75s   <- apply(abs(br4i75s[ids,,i]), 2, mean)
  
  m1i75s    <- apply((br1i75s[ids,,i])^2, 2, mean)
  m2i75s    <- apply((br2i75s[ids,,i])^2, 2, mean)
  m3i75s    <- apply((br3i75s[ids,,i])^2, 2, mean)
  m4i75s    <- apply((br4i75s[ids,,i])^2, 2, mean)
  
  p1i75s   <- abs(apply(pw1i75s[ids,,i], 2, mean, na.rm = TRUE))
  p2i75s   <- abs(apply(pw2i75s[ids,,i], 2, mean))
  p3i75s   <- abs(apply(pw3i75s[ids,,i], 2, mean, na.rm = TRUE))
  p4i75s   <- abs(apply(pw4i75s[ids,,i], 2, mean))
  
  w1i75s   <- apply(wd1i75s[ids,,i], 2, mean, na.rm = TRUE)
  w2i75s   <- apply(wd2i75s[ids,,i], 2, mean)
  w3i75s   <- apply(wd3i75s[ids,,i], 2, mean)
  w4i75s   <- apply(wd4i75s[ids,,i], 2, mean)
  
  c1i75s   <- abs(apply(cr1i75s[ids,,i], 2, mean))
  c2i75s   <- abs(apply(cr2i75s[ids,,i], 2, mean))
  c3i75s   <- abs(apply(cr3i75s[ids,,i], 2, mean))
  c4i75s   <- abs(apply(cr4i75s[ids,,i], 2, mean))
  
  bias1[,i]   <- b1i75s
  bias2[,i]   <- b2i75s
  bias3[,i]   <- b3i75s
  bias4[,i]   <- b4i75s
  
  mse1[,i]    <- m1i75s
  mse2[,i]    <- m2i75s
  mse3[,i]    <- m3i75s
  mse4[,i]    <- m4i75s
  
  power1[,i]  <- p1i75s
  power2[,i]  <- p2i75s
  power3[,i]  <- p3i75s
  power4[,i]  <- p4i75s
  
  width1[,i]  <- w1i75s
  width2[,i]  <- w2i75s
  width3[,i]  <- w3i75s
  width4[,i]  <- w4i75s
  
  cover1[,i]  <- c1i75s
  cover2[,i]  <- c2i75s
  cover3[,i]  <- c3i75s
  cover4[,i]  <- c4i75s
}

#### save output ####
# save.image(file = 'sim_K4_N75_mvp.RData')



