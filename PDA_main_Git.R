#******************************************************
#******************************************************
#******************************************************
# Paper: A Simple Procedure for Analyzing Reliability Data from Double-Stage Accelerated Life Tests

# ------------------------------------------------------------------------------------------
# Description:
# This code is used to perform the PDA procedure described in the paper
# Device-A data is used as an example

# Authors: Xiao Liu
# Last revision: 06/24/2019
# Contact: liuxiao314923@gmail.com
# Google page: https://sites.google.com/site/liuxiaosite1/
#******************************************************
#******************************************************
#******************************************************

# Sub-rountine: MLE for ALT data at a stress level (Weibull lifetime distribution with log-linear life-stress relationship)
like = function(para,data){  # negative log-likelihood
  y = data[,1]
  c = data[,2]
  mu = (para[1])  # para[1]:mu
  sigma = (para[2])  #para[2]:sigma
  yy = (y-mu)/sigma
  n = length(y)
  l = 0
  for (i in 1:n){
    if (c[i]==1){
      l = l + (-log(sigma)+yy[i]-exp(yy[i]))
    }else{
      l = l - exp(yy[i])
    }
  }
  l = -l
  return(l)
}

#******************************************************
#******************************************************
#******************************************************
# Data: Device-A data is used as an example
tmp1 = c(283, 361, 515, 638, 854, 1024, 1030, 1045, 1767, 1777, 1856, 1951, 1964, 2884,5000)
tmp2 = c( rep(1,14), 0)
data.h = cbind(log(tmp1), tmp2)

tmp1 = c(581, 95, 1432, 1586, 2452, 2734, 2772, 4106, 4674, rep(5000,11))
tmp2 = c( rep(1,9), rep(0,11))
data.m = cbind(log(tmp1), tmp2)

tmp1 = c(1298, 1390, 3187, 3241, 3261, 3313, 4501, 4568,4841, 4982, rep(5000,90))
tmp2 = c( rep(1,10), rep(0,90))
data.l = cbind(log(tmp1), tmp2)

# Stress levels:
S = c(37.0767, 34.8498, 32.8754)
gamma1 = 0.7


#******************************************************
#******************************************************
#******************************************************
# PDA

# Step 1:
output =optim(par=c(7,0.1),like,data=data.h,hessian=TRUE,method="L-BFGS-B",
      lower=c(-Inf,0.001),upper=c(Inf,Inf))
para.est.h = (output$par)
var.mat.h = solve(output$hessian)


# Step 3:
# mid stress
data.m.2 = cbind( data.h[,1] + gamma1 * (S[2]-S[3]), data.h[,2])
data.mm = rbind(data.m, data.m.2)
output =optim(par=c(7,0.1),like,data=data.mm,hessian=TRUE,method="L-BFGS-B",
              lower=c(-Inf,0.001),upper=c(Inf,Inf))
para.est.m = (output$par)
var.mat.m = solve(output$hessian)

# low stress
data.l.2 = cbind( data.h[,1] + gamma1 * (S[1]-S[3]), data.h[,2])
data.ll = rbind(data.l, data.l.2)
output =optim(par=c(10,0.1),like,data=data.ll,hessian=TRUE,method="L-BFGS-B",
              lower=c(-Inf,0.001),upper=c(Inf,Inf))
para.est.l = (output$par)
var.mat.l = solve(output$hessian)

# Step 4:
mu = c(para.est.l[1],para.est.m[1],para.est.h[1])
V = diag(c(var.mat.l[1,1], var.mat.m[1,1], var.mat.h[1,1]))
X = cbind( rep(1,length(S)), S)

# Equation (10):
output = solve( t(X)%*%solve(V)%*%X ) %*% ( t(X)%*%solve(V)%*% mu )
gamma0.est = output[1]
gamma1.est = output[2]

# Equation (11):
sigma.est = solve(  1/var.mat.l[2,2]+1/var.mat.m[2,2]+1/var.mat.h[2,2]  ) * 
  (  para.est.l[2]/var.mat.l[2,2]+para.est.m[2]/var.mat.m[2,2]+para.est.h[2]/var.mat.h[2,2]  ) 

# Equation (12):
var.mat =  rbind( cbind( solve( t(X)%*%solve(V)%*%X ), c(0,0) ), rep(0,3) )
var.mat[3,3] = solve(  1/var.mat.l[2,2]+1/var.mat.m[2,2]+1/var.mat.h[2,2]  )

# Equation (13)
y.est = gamma0.est + 11605/293 * gamma1.est + log(-log(1-0.1))*sigma.est


# output:
print(c(gamma0.est,gamma1.est,sigma.est))
print(y.est)







