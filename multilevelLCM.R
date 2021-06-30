# MLCM
#setwd("/Volumes/GoogleDrive/My Drive/BiostatMS/CIVIC")


rm(list=ls())
source("simData_cluster.R")

library(rjags)
library(R2jags)
library(coda)
library(gtools)

# 1.a define variable
N = nrow(data)
G = max(data[,1]) # this might have problem if some clusters are not sampled; minor issue.
K = max(data[,2]) # same comment as above.
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
      P[j,k] ~ dbeta(1,1) # here we are just using fixed hyperparameters (1,1); theoretically, we can also let them be random (I set them to be 5,5 in the simData_cluster.R)
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
#init.PI = matrix(1/K, nrow = G, ncol = K)
#init.PI[,2] = 1 - init.PI[,1]
model.inits <- function(){
  #set.seed(1080) # don't set seed here.
  list(
    P = matrix(runif(J*K, 0,1), nrow = J, ncol = K), # A|U JxK
    U = sample(1:K, size = N, replace = T,  prob = rep(1/K, K)),    # Latent Parameter 1/2/1/2
    PI =  matrix(1/K, nrow = G, ncol = K)
  )  # row sum to 1, randomly
}

params <- c("P","PI", "U")
res <- jags(data=dat,inits=model.inits,parameters.to.save = params,
            model.file=multilevel_LCM,n.chains=2,
            # since n.iter 2e4 doesn't converge, set 4e4
            n.iter=1000,n.burnin=500)
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
colMeans(df.P)
store.P = real.P
# store.P[2,] = real.P[10,]
# store.P[10,] = real.P[2,]
df.P = melt(df.P, variable_name = "parameter")

df.P$parameter <- factor(df.P$parameter,levels=c("P[1,1]",  "P[1,2]",  "P[2,1]",  "P[2,2]",  "P[3,1]",  "P[3,2]",  "P[4,1]",  "P[4,2]",  "P[5,1]", 
                                                 "P[5,2]","P[6,1]","P[6,2]","P[7,1]","P[7,2]","P[8,1]","P[8,2]","P[9,1]","P[9,2]","P[10,1]","P[10,2]"))

df.P$TrueVal = rep(c(t(store.P)), each = nrow(res.mcmc[[1]])) ## this needs to be modified to this.

df.P %>% ggplot(aes(x = parameter, y = value)) + geom_violin(scale = "width") + 
  stat_summary(fun.data = median_hilow, geom="pointrange") + 
  geom_point(aes(x = parameter, y = TrueVal), shape = 4, color = "#56B4E9") + 
  ggtitle(paste("Multi-level latent class model: G = ",G, ", K = ", K, ", N = ", N, ", J = ", J))

## 3.2.2 For PI: G*K
df.PI = as.data.frame(res.mcmc[[1]][,seq(K*J+2, K*J+G*K,length.out = G)])   
store.PI = real.PI
# store.PI[2,] = real.PI[10,]
# store.PI[10,] = real.PI[2,]
df.PI = melt(df.PI, variable_name = "parameter")
df.PI$parameter <- factor(df.PI$parameter,levels=c("PI[1,1]","PI[2,1]",  "PI[3,1]",  "PI[4,1]",  "PI[5,1]",  "PI[6,1]",  "PI[7,1]",  "PI[8,1]",  "PI[9,1]","PI[10,1]" ))

df.PI$TrueVal = rep(store.PI[,1], each = nrow(res.mcmc[[1]])) ## modified.

df.PI %>% ggplot(aes(x = parameter, y = value)) + geom_violin(scale = "width") + 
  stat_summary(fun.data = median_hilow, geom="pointrange") + 
  geom_point(aes(x = parameter, y = TrueVal), shape = 4, color = "#56B4E9") + # this is PI[g,1]
  #geom_point(aes(x = parameter, y = 1-TrueVal), shape = 4, color = "red") + # this is for 1-PI[g,1]
  ggtitle(paste("Multi-level latent class model: G = ",G, ", K = ", K, ", N = ", N, ", J = ", J))
