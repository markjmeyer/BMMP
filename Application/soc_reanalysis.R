#### Source Code ####
source('bmmp.R')

#### load data ####
pc		<- soc[,1:5] # primary care contact
colnames(pc)	<- c('PC_DD', 'PC_MH', 'PC_JJ', 'PC_CPS', 'PC_ED')

sc		<- soc[,6:10] # specialty care contact
colnames(sc)	<- c('SC_DD', 'SC_MH', 'SC_JJ', 'SC_CPS', 'SC_ED')

X1		<- pc
X2		<- sc
K     <- ncol(X1)

#### model specs ####
B         <- 10000
burnin    <- B

#### Bayesian MMP ####

model_bmvp    <- bmvp(X1, X2, B = B, burnin = burnin, pen = TRUE,
                      Al = 5, hc = TRUE, al = 0.01, bl = 0.01, 
                      Kp = 2, alpha = 0.001, As = 1, Bs = 1,
                      verbose = TRUE)
summary(model_bmvp)

#### trace plots ####
labsr    <- c(expression(rho[1]), expression(rho[2]), expression(rho[3]), expression(rho[4]), expression(rho[5]))
par(ask = TRUE)
for(i in 1:K){
  plot(1:B, model_bmvp$rho[,i], type = 'n', ylab = labsr[i], xlab = 'b', main = '')
  abline(h = axTicks(2), v = axTicks(1), lty = 3, col = 'lightgray')
  lines(1:B, model_bmvp$rho[,i], type = 'l', lwd = 2, col = 'forestgreen')
}

#### GEE and bootstrap-based approach ####
geek    <- geemmp(X1, X2)
boot    <- bootmmp(X1, X2, B = 10000)


