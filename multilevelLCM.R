# MLCM
setwd("/Volumes/GoogleDrive/My Drive/BiostatMS/CIVIC")

library(rjags)
library(R2jags)
library(coda)
library(gtools)

# 1.a define variable
N = nrow(data)
G = max(data[,1])
K = max(data[,2])
J = ncol(data) - 2


# 2. Model Specification 
######--------
multilevel_LCM <- function(){
  # likelihood 
  for(i in 1:N){
    U[i] ~ dcat(PI[Grp[i],])
    for (j in 1:J){
      A[i,j] ~ dbern(P[j,U[i]])
    }
  }
  
  # Priors
  ## PI - GxK
  for (g in 1:G){ 
    PI[g,1:K] ~ ddirch(alpha[1:K]) #specify vector/mat 
  } 
  ## P - A|U
  for (j in 1:J){ 
    for(k in 1:K){
      P[j,k] ~ dbeta(1,1)
    }
  }
}

dat <- list(
  'N' = N, 
  'J' = J,
  'K' = K,
  'G' = G,
  'alpha' = rep(1,K),
  'A' = as.matrix(data[,3:(2+J)]), # matrix
  'Grp' = data[,1]
)

# Sincce previously the estimation of PI which is latentclass U | cluster G is around 0.5
# set initial PI (0.2, 0.8) deviating from (0.5,0.5).
init.PI = matrix(0.2, nrow = G, ncol = K)
init.PI[,2] = 1 - init.PI[,1]
model.inits <- function(){
  set.seed(1080)
  list(
    P = matrix(runif(J*K, 0,1), nrow = J, ncol = K), # A|U JxK
    U = sample(1:K, size = N, replace = T,  prob = rep(1/K, K)),    # Latent Parameter 1/2/1/2
    PI =  init.PI
  )  # row sum to 1, randomly
}

params <- c("P","PI", "U")
res <- jags(data=dat,inits=model.inits,parameters.to.save = params,
            model.file=multilevel_LCM,n.chains=2,
            # since n.iter 2e4 doesn't converge, set 4e4
            n.iter=40000,n.burnin=30000)
res.mcmc = as.mcmc(res)
real.P = P
real.PI = PI
####--------------


# 3. Evalation 
######----------
library(lattice)
library(ggplot2)
library(dplyr)
library(reshape)

# 3.1 traceplot of some estimation
xyplot(
  res.mcmc[[1]][,2:27],
  layout = c(6,5))
xyplot(
  res.mcmc[[1]][,28:31])

# 3.2 Violinplot for evaluating estimation 

## 3.2.1 For P: J*K
df.P = as.data.frame(res.mcmc[[1]][,2:(K*J + 1)]) 
colMeans(df.par)
store.P = real.P
store.P[2,] = real.P[10,]
store.P[10,] = real.P[2,]
df.P = melt(df.par, variable_name = "parameter")

df.P$TrueVal = rep(t(store.P), each = N)
df.P %>% ggplot(aes(x = parameter, y = value)) + geom_violin(scale = "width") + 
  stat_summary(fun.data = median_hilow, geom="pointrange") + 
  geom_point(aes(x = parameter, y = TrueVal), shape = 4, color = "#56B4E9") + 
  ggtitle(paste("Multi-level latent class model: G = ",G, ", K = ", K, ", N = ", N, ", J = ", J))

## 3.2.2 For PI: G*K
df.PI = as.data.frame(res.mcmc[[1]][,seq(K*J+2, K*J+G*K,length.out = G)])   
store.PI = real.PI
store.PI[2] = real.PI[10]
store.PI[10] = real.PI[2]
df.PI = melt(df.PI, variable_name = "parameter")
df.PI$TrueVal = rep(store.PI, each = N)
df.PI %>% ggplot(aes(x = parameter, y = value)) + geom_violin(scale = "width") + 
  stat_summary(fun.data = median_hilow, geom="pointrange") + 
  geom_point(aes(x = parameter, y = TrueVal), shape = 4, color = "#56B4E9") + 
  ggtitle(paste("Multi-level latent class model: G = ",G, ", K = ", K, ", N = ", N, ", J = ", J))
