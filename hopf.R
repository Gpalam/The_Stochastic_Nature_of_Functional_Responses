Jacob <- function(alpha, beta, deltaA, nu){
  M <- matrix(0, 3, 3)
  M[1,1] = -beta*deltaA/alpha
  M[1,2] = -deltaA
  M[1,3] = 0
  M[2,1] = -beta*(1-deltaA/alpha)
  M[2,2] = -2*deltaA
  M[2,3] = 2*nu
  M[3,1] = beta*(1-deltaA/alpha)
  M[3,2] = deltaA
  M[3,3] = -nu
  return(M)
}

toSolve <- function(alpha, pars){
  J <- Jacob(alpha, pars[1], pars[2], pars[3])
  ev <- Re(eigen(J)$values)
  idx <- which.min(abs(ev))
  return(ev[idx])
}
pars <- c(0.001,1.5,0.1)
alphac <- 3*pars[2]+pars[3]^2/2/(pars[2]+pars[3])
toSolve(alphac,pars)

library(rootSolve)
uniroot(toSolve, c(alphac*0.9,alphac*1.1), pars)$root

# Hopf
nu <- 1
deltaA <- 5
pars <- c(0.001,deltaA,nu)
alphac <- 3*pars[2]+pars[3]^2/2/(pars[2]+pars[3])
end <- 50
dbeta <- 0.1
stab <- matrix(0, end/dbeta, 3)
i <- 1
for (beta in seq(0.001,end,dbeta)){
  pars <- c(beta,deltaA,nu)
  rt <- uniroot(toSolve, c(alphac*0.9,alphac*1.1), pars, tol = 1e-10)
  stab[i,1] <- rt$root
  stab[i,2] <- beta
  stab[i,3] <- rt$f.root
  alphac <- stab[i,1]
  i <- i+1
}
plot(stab[,1],stab[,2],xlim = c(0,20), type = "l", lwd = 3, lty = 2, xlab = "alpha", ylab = "beta", 
     main = "deltaR = 0, delta A = 5, lambdaR = 0, lambdaA = 0, nu = 1")
lines(rep(pars[2],length(stab[,1])),stab[,2], lwd = 3, lty = 2)
text(2,25,"Stable resource")
text(9,25,"Stable coexistence")

# Transcritical
pars <- c(0.001,deltaA,nu)
alphac <- 3*pars[2]+pars[3]^2/2/(pars[2]+pars[3])
J <- Jacob(5.1, pars[1], pars[2], pars[3])
eigen(J)$values

# Oscillations
toSolve1 <- function(alpha, pars){
  J <- Jacob(alpha, pars[1], pars[2], pars[3])
  ev <- Im(eigen(J)$values)
  idx <- which(ev>0)
  if (length(idx)==0){
    ret = -1
  } else {
    ret = abs(ev[idx])
  }
  return(ret)
}

i <- 20
pars <- c(stab[i,2],deltaA,nu)
uniroot(toSolve1, c(deltaA,stab[i,1]), pars)$root

stab1 <- matrix(0, end/dbeta, 3)
i <- 1
for (beta in seq(0.001,end,dbeta)){
  pars <- c(beta,deltaA,nu)
  rt <- uniroot(toSolve1, c(deltaA,stab[i,1]), pars, tol = 1e-10)
  stab1[i,1] <- rt$root
  stab1[i,2] <- beta
  stab1[i,3] <- rt$f.root
  i <- i+1
}

X11()

pp <- 5
cex <- 2
par(mar=c(pp,pp,pp,pp))
plot(stab[,1],stab[,2],xlim = c(0,20), type = "l", lwd = 3, lty = 2, 
     cex.lab=2,cex.axis=2,xlab = expression(alpha), ylab = expression(beta), 
     main = "")
lines(rep(pars[2],length(stab[,1])),stab[,2], lwd = 3, lty = 2)
lines(stab1[,1],stab1[,2], lwd = 3, lty = 2)
text(2,25,"Stable \n resource",cex=cex)
text(7.5,27,"Stable",cex=cex)
text(7.5,23,"coexistence",cex=cex)
text(12.5,27,"Damped",cex=cex)
text(12.5,23,"oscillations",cex=cex)
text(18,27,"Limit",cex=cex)
text(18,23,"cycles",cex=cex)

#
pars
J <- Jacob(12, pars[1], pars[2], pars[3])
eigen(J)$values

